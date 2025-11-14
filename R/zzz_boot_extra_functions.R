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
# -------------------------------------------------------------------
# Bootstrap wrapper around fit_msm()
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Bootstrap wrapper around fit_msm()  — id_col & weights first-class
# -------------------------------------------------------------------
fit_msm_boot <- function(
    data,
    id_col            = "hhidpn",
    strat_vars        = NULL,
    covariate_var     = NULL,
    age_from_to       = c(50, 100),
    age_int           = 0.25,
    spline_df         = NULL,
    spline_type       = "ns",
    Q,
    times             = 8,
    weight_col        = NULL,
    n_cores           = 1,
    ci_level          = 0.95,
    return_replicates = FALSE,
    parallel          = c("auto","none","future","mclapply")
) {
  
  parallel <- match.arg(parallel)
  
  # ---- 0) Master prediction grid so all replicates align ----
  master_grid <- build_prediction_grid(
    data         = data,
    strat_vars   = strat_vars,
    covariate_var= covariate_var,
    age_from_to  = age_from_to,
    age_int      = age_int
  ) |>
    tidyr::crossing(
      from = seq_len(nrow(Q)),
      to   = seq_len(ncol(Q))
    ) |>
    dplyr::mutate(
      from = as.integer(from),
      to   = as.integer(to)
    ) |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(strat_vars), as.character),
      dplyr::across(dplyr::all_of(covariate_var), as.character),
      age = round(as.numeric(age), 6)
    )
  
  join_keys <- unique(c(strat_vars, covariate_var, "age", "from", "to"))
  master_grid <- master_grid |>
    dplyr::distinct(dplyr::across(dplyr::all_of(join_keys)), .keep_all = TRUE)
  
  # ---- 1) Stratified, group bootstrap replicates (weighted if provided) ----
  # Prefer your rsample-based infra if present; fallback to a weighted builder.
  boot_tbl <- if (exists("group_bootstraps2", mode = "function")) {
    rset <- group_bootstraps2(
      data   = data,
      group  = id_col,
      weight = weight_col,
      times  = times
    )
    tibble::tibble(
      replicate = seq_along(rset$splits),
      data      = lapply(rset$splits, rsample::analysis)
    )
  } else {
    # Fallback: stratified, weighted group bootstrap
    #  - sample IDs with replacement within each stratum
    #  - size equals # unique IDs in that stratum
    #  - prob ∝ subject-level mean(weight_col) or equal if NULL
    if (!length(id_col) || !id_col %in% names(data)) {
      stop("id_col must be a column in `data`.")
    }
    strata_vars <- strat_vars %||% character(0)
    
    id_strata <- data |>
      dplyr::distinct(dplyr::across(dplyr::all_of(c(id_col, strata_vars))))
    
    if (!is.null(weight_col) && weight_col %in% names(data)) {
      id_w <- data |>
        dplyr::group_by(.data[[id_col]]) |>
        dplyr::summarise(.w = mean(.data[[weight_col]], na.rm = TRUE), .groups="drop")
      id_strata <- dplyr::left_join(id_strata, id_w, by = id_col)
      id_strata$.w[is.na(id_strata$.w)] <- 1
    } else {
      id_strata$.w <- 1
    }
    
    split_ids <- if (length(strata_vars)) {
      id_strata |>
        dplyr::group_split(dplyr::across(dplyr::all_of(strata_vars)), .keep = TRUE)
    } else list(id_strata)
    
    sample_ids_once <- function() {
      parts <- lapply(split_ids, function(tbl) {
        n  <- nrow(tbl)
        p  <- tbl$.w / sum(tbl$.w)
        ix <- sample.int(n, n, replace = TRUE, prob = p)
        tbl[ix, , drop = FALSE]
      })
      dplyr::bind_rows(parts)
    }
    
    reps <- vector("list", times)
    for (t in seq_len(times)) {
      samp_ids <- sample_ids_once()[[id_col]]
      # bind rows for each selected id occurrence, relabeling duplicates on the fly
      out <- lapply(seq_along(samp_ids), function(k) {
        this_id <- samp_ids[k]
        rows <- data[data[[id_col]] == this_id, , drop = FALSE]
        occ_k <- sum(samp_ids[seq_len(k)] == this_id)
        rows[[id_col]] <- if (occ_k == 1L) this_id else paste0(this_id, "_", occ_k)
        rows
      })
      reps[[t]] <- dplyr::bind_rows(out)
    }
    tibble::tibble(replicate = seq_len(times), data = reps)
  }
  
  # ---- 2) Per-replicate safe fit (drops constant covariates if needed) ----
  safe_fit_one <- function(df) {
    # Order within subject by time
    df <- df |>
      dplyr::arrange(.data[[id_col]], age) |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::arrange(age, .by_group = TRUE) |>
      dplyr::ungroup()
    
    # Attempt fit; if contrasts fail (single-level factor), drop constant covariates
    try_once <- function(covars) {
      fit_msm(
        data          = df,
        strat_vars    = strat_vars,
        covariate_var = covars,
        age_from_to   = age_from_to,
        age_int       = age_int,
        spline_df     = spline_df,
        spline_type   = spline_type,
        Q             = Q
      )
    }
    
    out <- tryCatch(try_once(covariate_var), error = identity)
    
    if (inherits(out, "error")) {
      cov_keep <- covariate_var
      if (length(covariate_var)) {
        vary <- vapply(covariate_var, function(v) {
          if (v %in% names(df)) dplyr::n_distinct(df[[v]]) > 1 else FALSE
        }, logical(1))
        cov_keep <- covariate_var[vary]
      }
      out <- tryCatch(try_once(cov_keep), error = identity)
    }
    
    if (inherits(out, "error") || !nrow(out)) {
      master_grid |>
        dplyr::mutate(haz = NA_real_)
    } else {
      out |>
        dplyr::mutate(
          dplyr::across(dplyr::all_of(strat_vars), as.character),
          dplyr::across(dplyr::all_of(covariate_var), as.character),
          age  = round(as.numeric(age), 6),
          from = as.integer(from),
          to   = as.integer(to)
        ) |>
        dplyr::right_join(master_grid, by = join_keys) |>
        dplyr::mutate(haz = as.numeric(haz))
    }
  }
  
  # ---- 3) Run fits over replicates (optional parallel) ----
  run_indices <- seq_len(nrow(boot_tbl))
  run_one <- function(i) {
    out <- safe_fit_one(boot_tbl$data[[i]])
    out$replicate <- i
    out
  }
  
  os_unix    <- .Platform$OS.type == "unix"
  in_rstudio <- identical(Sys.getenv("RSTUDIO"), "1")
  
  # Determine backend
  if (parallel == "none") {
    
    message("▶ Running sequentially on 1 core.")
    replicates_df <- purrr::map_dfr(run_indices, run_one, .id = NULL)
    
  } else if (parallel == "mclapply") {
    
    if (!os_unix) {
      stop("mclapply is only available on Unix (Linux/macOS).")
    }
    message("▶ Running with parallel::mclapply on ", n_cores, " cores.")
    parts <- parallel::mclapply(run_indices, run_one, mc.cores = n_cores)
    replicates_df <- dplyr::bind_rows(parts)
    
  } else if (parallel %in% c("future","auto")) {
    
    # Decide which future backend we want
    if (parallel == "auto") {
      # auto mode chooses multicore *only* when safe
      if (os_unix && !in_rstudio) {
        desired_plan <- "multicore"
      } else {
        desired_plan <- "multisession"
      }
    } else {
      # parallel == "future": if user already set a plan, keep it
      cur_plan <- future::plan()
      plan_name <- class(cur_plan)[1L]
      
      if (!identical(plan_name, "sequential")) {
        # use user’s existing plan
        desired_plan <- plan_name
      } else {
        # decide based on environment
        if (os_unix && !in_rstudio) {
          desired_plan <- "multicore"
        } else {
          desired_plan <- "multisession"
        }
      }
    }
    
    # Apply the plan if needed
    if (desired_plan == "multicore") {
      message("▶ Setting future plan: multicore (", n_cores, " workers).")
      future::plan(future::multicore, workers = n_cores)
    } else if (desired_plan == "multisession") {
      message("▶ Setting future plan: multisession (", n_cores, " workers).")
      future::plan(future::multisession, workers = n_cores)
    } else {
      message("▶ Using existing user-defined future plan: ", desired_plan)
    }
    
    # guarantee cleanup
    on.exit(try(future::plan("sequential"), silent = TRUE), add = TRUE)
    
    # run in parallel
    replicates_df <- furrr::future_map_dfr(
      run_indices,
      run_one,
      .options = furrr::furrr_options(seed = TRUE)
    )
    
  } else {
    
    stop("Unknown parallel mode: ", parallel)
    
  }
  
  
  # ---- 4) Sanity: do we have variation across replicates? ----
  # var_check <- replicates_df |>
  #   dplyr::group_by(dplyr::across(dplyr::all_of(join_keys))) |>
  #   dplyr::summarise(unique_haz = dplyr::n_distinct(.data$haz[!is.na(.data$haz)]),
  #                    .groups = "drop")
  # 
  # message("✔ Replicate variation summary (unique hazard values per cell):")
  # print(summary(var_check$unique_haz))
  
  # ---- 5) Summarise to bootstrap CIs ----
  # ---- 5) Summarise to bootstrap CIs (renamed columns) ----
  alpha <- (1 - ci_level) / 2
  lo <- alpha
  hi <- 1 - alpha
  
  summary_df <- replicates_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(join_keys))) |>
    dplyr::summarise(
      hazard = mean(.data$haz, na.rm = TRUE),
      lower  = stats::quantile(.data$haz, probs = lo, na.rm = TRUE, type = 7),
      upper  = stats::quantile(.data$haz, probs = hi, na.rm = TRUE, type = 7),
      n_reps = sum(!is.na(.data$haz)),
      .groups = "drop"
    )
  
  if (isTRUE(return_replicates)) {
    return(list(summary = summary_df, replicates = replicates_df))
  } else {
    return(summary_df)
  }
  
}













