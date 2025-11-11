## ============================================================
## helpers
## ============================================================

build_prediction_grid <- function(data,
                                  strat_vars,
                                  covariate_var,
                                  age_from_to,
                                  age_int) {
  
  base_grid <- list(age = seq(age_from_to[1], age_from_to[2], by = age_int))
  
  cov_vals <- if (!is.null(covariate_var)) {
    purrr::map(covariate_var, ~ unique(data[[.x]])) |>
      rlang::set_names(covariate_var)
  } else list()
  
  strat_vals <- if (!is.null(strat_vars)) {
    purrr::map(strat_vars, ~ unique(data[[.x]])) |>
      rlang::set_names(strat_vars)
  } else list()
  
  rlang::exec(tidyr::crossing, !!!c(base_grid, cov_vals, strat_vals))
}

add_splines <- function(data,
                        age_var    = "age",
                        spline_type = "ns",
                        spline_df   = 4) {
  
  age_vals <- data[[age_var]]
  
  if (identical(spline_type, "ns")) {
    basis_fit <- splines::ns(age_vals, df = spline_df)
  } else if (identical(spline_type, "bs")) {
    basis_fit <- splines::bs(age_vals, df = spline_df)
  } else {
    stop("spline_type must be 'ns' or 'bs'")
  }
  
  basis_mat <- predict(basis_fit, newx = age_vals)
  basis_tbl <- tibble::as_tibble(basis_mat)
  colnames(basis_tbl) <- paste0("age_spline", seq_len(ncol(basis_tbl)))
  
  data_out <- dplyr::bind_cols(data, basis_tbl)
  
  list(
    data      = data_out,
    basis_fit = basis_fit
  )
}

add_splines_to_grid <- function(pred_grid, basis_fit) {
  basis_mat <- predict(basis_fit, newx = pred_grid$age)
  basis_tbl <- tibble::as_tibble(basis_mat)
  colnames(basis_tbl) <- paste0("age_spline", seq_len(ncol(basis_tbl)))
  dplyr::bind_cols(pred_grid, basis_tbl)
}

msm_model <- function(data, Q, covariate_var, extra_covars = NULL) {
  all_covs <- c(covariate_var, extra_covars)
  all_covs <- all_covs[!is.null(all_covs)]
  
  if (length(all_covs) == 0) {
    cov_form <- NULL
  } else {
    cov_form <- stats::reformulate(all_covs)
  }
  
  suppressWarnings(
    msm::msm(
      formula    = state_msm ~ age,
      subject    = hhidpn,
      data       = data,
      qmatrix    = Q,
      obstype    = obstype,
      deathexact = 3,
      covariates = cov_form,
      control    = list(fnscale = 5000, maxit = 25000),
      gen.inits  = TRUE,
      method     = "BFGS"
    )
  )
}

safe_covariate_list <- function(x, wanted = NULL) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) {
    x <- as.list(x[1, , drop = FALSE])
  }
  if (!is.null(wanted)) {
    x <- x[wanted]
  }
  x
}

tidy_q <- function(qmat) {
  if (is.null(qmat)) {
    return(tibble::tibble(from = NA_integer_, to = NA_integer_, estimate = NA_real_))
  }
  if (!is.matrix(qmat)) {
    stop("Expected matrix in tidy_q(), got: ", class(qmat)[1])
  }
  
  df <- as.data.frame(as.table(qmat)) %>%
    dplyr::rename(from = Var1, to = Var2, estimate = Freq)
  
  # first attempt: parse numbers out of the labels
  from_num <- suppressWarnings(
    as.integer(readr::parse_number(as.character(df$from)))
  )
  to_num <- suppressWarnings(
    as.integer(readr::parse_number(as.character(df$to)))
  )
  
  # if ALL are NA, it means labels were non-numeric ("A","B","C") or missing
  if (all(is.na(from_num))) {
    # figure out matrix size from qmat directly
    n <- nrow(qmat)
    # repeat 1:n for each column
    from_num <- rep(seq_len(n), times = n)
  }
  if (all(is.na(to_num))) {
    n <- ncol(qmat)
    # repeat each 1:n for each row
    to_num <- rep(seq_len(n), each = n)
  }
  
  df$from <- from_num
  df$to   <- to_num
  
  df
}

## ============================================================
## main fitter (non-boot)
## ============================================================
fit_msm <- function(
    data,
    strat_vars    = NULL,
    covariate_var = NULL,
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    spline_df     = NULL,
    spline_type   = "ns",
    Q             = default_q
) {
  
  ## 0) drop NAs in covariates/strata to avoid "Covariate NA unknown"
  drop_vars <- c(strat_vars, covariate_var)
  if (length(drop_vars)) {
    data <- data %>%
      dplyr::filter(
        dplyr::if_all(dplyr::all_of(drop_vars), ~ !is.na(.x))
      )
  }
  
  ## 1) maybe add splines to the DATA
  spline_names <- character(0)
  basis_fit    <- NULL
  if (!is.null(spline_df)) {
    spl_out      <- add_splines(
      data,
      age_var     = "age",
      spline_type = spline_type,
      spline_df   = spline_df
    )
    data         <- spl_out$data
    basis_fit    <- spl_out$basis_fit
    spline_names <- grep("^age_spline", names(data), value = TRUE)
  }
  
  ## 2) build flat prediction grid from this data
  pred_grid <- build_prediction_grid(
    data         = data,
    strat_vars   = strat_vars,
    covariate_var= covariate_var,
    age_from_to  = age_from_to,
    age_int      = age_int
  )
  
  ## 3) add splines to the prediction grid too
  if (!is.null(spline_df)) {
    pred_grid <- add_splines_to_grid(pred_grid, basis_fit)
  }
  
  ## 4) now nest so every row has all covariates (including splines)
  nest_vars <- c(strat_vars, covariate_var, "age")
  pred_grid <- pred_grid %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(nest_vars))) %>%
    dplyr::group_nest(.key = "data_fit") %>%
    dplyr::mutate(
      data_fit = purrr::map(.data$data_fit, ~ as.list(.x[1, , drop = FALSE]))
    )
  
  ## 5) fit per stratum in the original data
  out <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(strat_vars))) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(group_data) {
      
      # current stratum (e.g. female = 0/1)
      strat_vals <- group_data[1, strat_vars, drop = FALSE]
      
      # fit msm with ALL covariates (period + spline cols)
      mod <- msm_model(
        data         = group_data,
        Q            = Q,
        covariate_var= covariate_var,
        extra_covars = spline_names
      )
      
      # subset prediction grid to this stratum
      match_rows <- rep(TRUE, nrow(pred_grid))
      for (v in strat_vars) {
        match_rows <- match_rows & (pred_grid[[v]] == strat_vals[[v]])
      }
      pred_sub <- pred_grid[match_rows, ]
      
      # evaluate Q for every row in this (stratum-specific) grid
      pred_sub <- pred_sub %>%
        dplyr::mutate(
          q = purrr::map(
            .data$data_fit,
            ~{
              covs <- safe_covariate_list(
                .x,
                wanted = c(covariate_var, spline_names)
              )
              
              # robust check for missing/empty covariates
              bad_cov <- FALSE
              if (length(covs) == 0L) {
                bad_cov <- TRUE
              } else {
                bad_cov <- any(vapply(
                  covs,
                  function(v) {
                    # empty
                    if (length(v) == 0L) return(TRUE)
                    # plain NA
                    if (is.na(v)) return(TRUE)
                    # character NA / "NA"
                    if (is.character(v) && (length(v) == 0L ||
                                            is.na(v) ||
                                            v == "NA")) return(TRUE)
                    FALSE
                  },
                  logical(1)
                ))
              }
              
              if (bad_cov) {
                return(matrix(NA_real_, nrow = nrow(Q), ncol = ncol(Q)))
              }
              
              # get qmatrix; suppress benign msm warnings
              qmat <- tryCatch(
                suppressWarnings(
                  msm::qmatrix.msm(mod, covariates = covs, ci = "none")
                ),
                error = function(e)
                  matrix(NA_real_, nrow = nrow(Q), ncol = ncol(Q))
              )
              
              if (is.matrix(qmat)) qmat else qmat$estimate
            }
          )
        ) %>%
        dplyr::select(-data_fit)
      
      pred_sub
    })
  
  ## 6) tidy to long with hazards
  out %>%
    dplyr::mutate(haz = purrr::map(q, tidy_q)) %>%
    tidyr::unnest(haz) %>%
    dplyr::select(
      dplyr::all_of(c(strat_vars, covariate_var,
                      "age", "from", "to", "estimate"))
    ) %>%
    dplyr::rename(haz = "estimate") %>%
    dplyr::mutate(
      from = readr::parse_number(as.character(.data$from)) |> as.integer(),
      to   = readr::parse_number(as.character(.data$to))   |> as.integer()
    )
}


## ============================================================
## boot wrapper (this is the one your run is probably picking up
## from zzz_boot_test_functions.R, so overwrite it with this)
## ============================================================
fit_msm_boot <- function(
    data,
    strat_vars    = NULL,
    covariate_var = NULL,
    times         = 8,
    weight        = "pwt",
    group         = "hhidpn",
    n_cores       = 4,
    Q,
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    spline_type   = "ns",
    spline_df     = NULL
) {
  
  boot_rset <- group_bootstraps2(
    data   = data,
    group  = group,
    weight = weight,
    times  = times
  )
  
  if (Sys.getenv("RSTUDIO") != 1 & .Platform$OS.type == "unix") {
    future::plan("multicore", workers = n_cores)
  } else {
    future::plan("multisession", workers = n_cores)
  }
  
  cat(" Running in parallel with", future::nbrOfWorkers(), "worker(s)\n")
  
  # --- prebuild reference grid from FULL dataset ---
  reference_grid <- build_prediction_grid(
    data         = data,
    strat_vars   = strat_vars,
    covariate_var= covariate_var,
    age_from_to  = age_from_to,
    age_int      = age_int
  )
  # include spline structure if needed
  if (!is.null(spline_df)) {
    spl_out <- add_splines(data,
                           age_var     = "age",
                           spline_type = spline_type,
                           spline_df   = spline_df)
    reference_grid <- add_splines_to_grid(reference_grid, spl_out$basis_fit)
  }
  
  results <- furrr::future_map_dfr(
    boot_rset$splits,
    ~{
      d_resample <- rsample::analysis(.x) |> dplyr::arrange(hhidpn, age)
      fit_msm(
        data         = d_resample,
        strat_vars   = strat_vars,
        covariate_var= covariate_var,
        age_from_to  = age_from_to,
        age_int      = age_int,
        spline_type  = spline_type,
        spline_df    = spline_df,
        Q            = Q
      )
    },
    .options = furrr::furrr_options(
      seed     = TRUE,
      packages = c("dplyr", "tidyr", "purrr", "msm", "readr", "tibble"),
      globals  = list(
        fit_msm             = fit_msm,
        msm_model           = msm_model,
        tidy_q              = tidy_q,
        build_prediction_grid = build_prediction_grid,
        add_splines         = add_splines,
        add_splines_to_grid = add_splines_to_grid,
        safe_covariate_list = safe_covariate_list,
        group_bootstraps2   = group_bootstraps2,
        reference_grid      = reference_grid,
        Q                   = Q,
        age_from_to         = age_from_to,
        age_int             = age_int,
        spline_type         = spline_type,
        spline_df           = spline_df,
        strat_vars          = strat_vars,
        covariate_var       = covariate_var
      )
    ),
    .progress = TRUE,
    .id       = "replicate"
  )
  
  future::plan("sequential")
  gc()
  
  results %>%
    dplyr::group_by(dplyr::across(-c(.data$replicate, .data$haz))) %>%
    dplyr::summarise(
      haz   = mean(.data$haz, na.rm = TRUE),
      haz_l = stats::quantile(.data$haz, 0.025, na.rm = TRUE),
      haz_u = stats::quantile(.data$haz, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

