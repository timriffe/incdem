## ============================================================
## helpers
## ============================================================

build_prediction_grid <- function(data,
                                  strat_vars,
                                  covariate_var,
                                  age_from_to,
                                  age_int) {
  
  # base age grid
  base_grid <- list(age = seq(age_from_to[1], age_from_to[2], by = age_int))
  
  # don't re-add age if user put "age" in covariate_var
  covariate_var_nage <- setdiff(covariate_var, "age")
  
  # values for covariates (except age)
  cov_vals <- if (!is.null(covariate_var_nage) && length(covariate_var_nage)) {
    purrr::map(covariate_var_nage, ~ unique(data[[.x]])) |>
      rlang::set_names(covariate_var_nage)
  } else {
    list()
  }
  
  # values for strat vars
  strat_vals <- if (!is.null(strat_vars) && length(strat_vars)) {
    purrr::map(strat_vars, ~ unique(data[[.x]])) |>
      rlang::set_names(strat_vars)
  } else {
    list()
  }
  
  # full cartesian grid
  rlang::exec(
    tidyr::crossing,
    !!!c(base_grid, cov_vals, strat_vals)
  )
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
  # 0) drop NAs in covariates/strata
  drop_vars <- c(strat_vars, covariate_var)
  if (length(drop_vars)) {
    data <- data %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(drop_vars), ~ !is.na(.x)))
  }
  
  # 1) maybe add splines to DATA
  spline_names <- character(0)
  basis_fit    <- NULL
  if (!is.null(spline_df)) {
    spl_out <- add_splines(
      data,
      age_var     = "age",
      spline_type = spline_type,
      spline_df   = spline_df
    )
    data         <- spl_out$data
    basis_fit    <- spl_out$basis_fit
    spline_names <- grep("^age_spline", names(data), value = TRUE)
  }
  
  # 2) build prediction grid
  pred_grid <- build_prediction_grid(
    data         = data,
    strat_vars   = strat_vars,
    covariate_var= covariate_var,
    age_from_to  = age_from_to,
    age_int      = age_int
  )
  
  # 3) add splines to GRID
  if (!is.null(spline_df)) {
    pred_grid <- add_splines_to_grid(pred_grid, basis_fit)
  }
  
  # 4) fit per stratum, evaluate row-by-row
  out <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(strat_vars))) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(group_data) {
      
      # current stratum values
      strat_vals <- if (length(strat_vars)) group_data[1, strat_vars, drop = FALSE] else NULL
      
      # fit msm for this stratum
      mod <- msm_model(
        data         = group_data,
        Q            = Q,
        covariate_var= covariate_var,
        extra_covars = spline_names
      )
      
      # subset prediction grid to this stratum
      pred_sub <- pred_grid
      if (length(strat_vars)) {
        for (v in strat_vars) {
          pred_sub <- pred_sub[pred_sub[[v]] == strat_vals[[v]], , drop = FALSE]
        }
      }
      
      # align factor covariates in grid to model data
      if (!is.null(covariate_var)) {
        for (cv in covariate_var) {
          if (cv %in% names(group_data) && is.factor(group_data[[cv]])) {
            pred_sub[[cv]] <- factor(pred_sub[[cv]], levels = levels(group_data[[cv]]))
          }
        }
      }
      
      # evaluate q for each row
      pred_sub$q <- purrr::map(seq_len(nrow(pred_sub)), function(i) {
        row_i <- pred_sub[i, , drop = FALSE]
        
        cov_cols <- c(covariate_var, spline_names)
        covs <- NULL
        if (length(cov_cols)) {
          covs <- lapply(cov_cols, function(cc) row_i[[cc]][[1]])
          names(covs) <- cov_cols
        }
        
        # check for bad/missing covariates
        bad_cov <- FALSE
        if (!is.null(covs) && length(covs)) {
          for (cc in names(covs)) {
            v <- covs[[cc]]
            if (length(v) == 0L ||
                is.na(v) ||
                (is.character(v) && (is.na(v) || v == "" || v == "NA"))) {
              bad_cov <- TRUE
              break
            }
          }
        }
        
        if (bad_cov) {
          return(matrix(NA_real_, nrow = nrow(Q), ncol = ncol(Q)))
        }
        
        qmat <- tryCatch(
          suppressWarnings(
            msm::qmatrix.msm(mod, covariates = covs, ci = "none")
          ),
          error = function(e)
            matrix(NA_real_, nrow = nrow(Q), ncol = ncol(Q))
        )
        
        if (is.matrix(qmat)) qmat else qmat$estimate
      })
      
      pred_sub
    })
  
  # 5) tidy
  out %>%
    dplyr::mutate(haz = purrr::map(q, tidy_q)) %>%
    tidyr::unnest(haz) %>%
    dplyr::select(dplyr::all_of(c(
      strat_vars,
      covariate_var,
      "age",
      "from",
      "to",
      "estimate"
    ))) %>%
    dplyr::rename(haz = "estimate")
}


# helper: make sure a resample has all levels of the categorical covariates
ensure_levels_in_resample <- function(resample,
                                      full_data,
                                      cat_vars,
                                      n_per_level = 20) {
  # resample: the bootstrap sample we just drew
  # full_data: the original full dataset (has all levels)
  # cat_vars: character vector of vars we want all levels for
  # n_per_level: how many rows to borrow per missing level
  out <- resample
  
  for (v in cat_vars) {
    full_lvls <- unique(full_data[[v]])
    rs_lvls   <- unique(out[[v]])
    
    missing_lvls <- setdiff(full_lvls, rs_lvls)
    
    if (length(missing_lvls)) {
      # for each missing level, borrow some rows from the full data
      for (ml in missing_lvls) {
        pool <- full_data[full_data[[v]] == ml, , drop = FALSE]
        # if pool is small, just take all
        if (nrow(pool) <= n_per_level) {
          add <- pool
        } else {
          add <- pool[sample.int(nrow(pool), n_per_level, replace = TRUE), , drop = FALSE]
        }
        out <- dplyr::bind_rows(out, add)
      }
    }
  }
  
  out
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
  
  # detect categorical covariates from full data
  cat_covs   <- character(0)
  cov_levels <- list()
  if (!is.null(covariate_var)) {
    for (cv in covariate_var) {
      if (is.factor(data[[cv]]) || is.character(data[[cv]])) {
        cat_covs <- c(cat_covs, cv)
        cov_levels[[cv]] <- unique(data[[cv]])
      }
    }
  }
  
  # 1) master prediction grid (common to all reps)
  master_grid <- build_prediction_grid(
    data         = data,
    strat_vars   = strat_vars,
    covariate_var= covariate_var,
    age_from_to  = age_from_to,
    age_int      = age_int
  ) %>%
    tidyr::crossing(
      from = seq_len(nrow(Q)),
      to   = seq_len(ncol(Q))
    ) %>%
    dplyr::mutate(
      from = as.integer(from),
      to   = as.integer(to),
      dplyr::across(dplyr::all_of(strat_vars), as.character),
      dplyr::across(dplyr::all_of(covariate_var), as.character),
      age  = round(as.numeric(age), 3)
    )
  
  # make join keys unique even if "age" is in covariate_var
  join_keys <- unique(c(strat_vars, covariate_var, "age", "from", "to"))
  
  # 2) bootstrap splits
  boot_rset <- group_bootstraps2(
    data   = data,
    group  = group,
    weight = weight,
    times  = times
  )
  
  # 3) parallel plan
  if (Sys.getenv("RSTUDIO") != 1 & .Platform$OS.type == "unix") {
    future::plan("multicore", workers = n_cores)
  } else {
    future::plan("multisession", workers = n_cores)
  }
  
  cat(" Running in parallel with", future::nbrOfWorkers(), "worker(s)\n")
  
  results <- furrr::future_map_dfr(
    seq_along(boot_rset$splits),
    ~{
      sp <- boot_rset$splits[[.x]]
      d_resample <- rsample::analysis(sp)
      
      # 3a) top up missing categorical levels so msm can build contrasts
      if (length(cat_covs)) {
        d_resample <- ensure_levels_in_resample(
          resample     = d_resample,
          full_data    = data,
          cat_vars     = cat_covs,
          n_per_level  = 20
        )
      }
      
      # 3b) robust relabel: every repeated id becomes a separate subject
      d_resample <- d_resample %>%
        dplyr::group_by(.data[[group]]) %>%
        dplyr::mutate(.dup_id = dplyr::row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          !!group := dplyr::if_else(
            .dup_id == 1L,
            .data[[group]],
            paste0(.data[[group]], "_", .dup_id)
          )
        ) %>%
        dplyr::select(-.dup_id)
      
      # 3c) order within subject
      d_resample <- d_resample %>%
        dplyr::arrange(.data[[group]], age) %>%
        dplyr::group_by(.data[[group]]) %>%
        dplyr::arrange(age, .by_group = TRUE) %>%
        dplyr::ungroup()
      
      # 3d) re-impose factor levels from full data
      if (length(cat_covs)) {
        for (cv in cat_covs) {
          if (cv %in% names(d_resample)) {
            d_resample[[cv]] <- factor(d_resample[[cv]],
                                       levels = cov_levels[[cv]])
          }
        }
      }
      
      # 3e) safe fit
      out_i <- tryCatch(
        {
          fit_msm(
            data         = d_resample,
            strat_vars   = strat_vars,
            covariate_var= covariate_var,
            age_from_to  = age_from_to,
            age_int      = age_int,
            spline_type  = spline_type,
            spline_df    = spline_df,
            Q            = Q
          ) %>%
            dplyr::mutate(
              dplyr::across(dplyr::all_of(strat_vars), as.character),
              dplyr::across(dplyr::all_of(covariate_var), as.character),
              age  = round(as.numeric(age), 3),
              from = as.integer(from),
              to   = as.integer(to)
            )
        },
        error = function(e) {
          # fallback: empty hazards for this replicate
          master_grid %>%
            dplyr::mutate(haz = NA_real_)
        }
      )
      
      # 3f) align to master grid
      master_grid %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(join_keys)), .keep_all = TRUE) %>%
        dplyr::left_join(out_i, by = join_keys) %>%
        dplyr::mutate(replicate = .x)
    },
    .options = furrr::furrr_options(
      seed     = TRUE,
      packages = c("dplyr","tidyr","purrr","msm","readr","tibble"),
      globals  = list(
        fit_msm               = fit_msm,
        msm_model             = msm_model,
        tidy_q                = tidy_q,
        build_prediction_grid = build_prediction_grid,
        add_splines           = add_splines,
        add_splines_to_grid   = add_splines_to_grid,
        safe_covariate_list   = safe_covariate_list,
        group_bootstraps2     = group_bootstraps2,
        ensure_levels_in_resample = ensure_levels_in_resample,
        Q                     = Q,
        age_from_to           = age_from_to,
        age_int               = age_int,
        spline_type           = spline_type,
        spline_df             = spline_df,
        strat_vars            = strat_vars,
        covariate_var         = covariate_var,
        cat_covs              = cat_covs,
        cov_levels            = cov_levels,
        group                 = group,
        master_grid           = master_grid,
        join_keys             = join_keys
      )
    ),
    .progress = TRUE,
    .id       = "replicate"
  )
  
  future::plan("sequential")
  gc()
  
  # 4) summarise
  results %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(join_keys))) %>%
    dplyr::summarise(
      haz     = mean(haz, na.rm = TRUE),
      haz_l   = stats::quantile(haz, 0.025, na.rm = TRUE),
      haz_u   = stats::quantile(haz, 0.975, na.rm = TRUE),
      n_reps  = sum(!is.na(haz)),
      .groups = "drop"
    )
}











