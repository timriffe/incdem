
# whenever we add packages, after installing locally, run renv::snapshot()
# renv::snapshot(force=TRUE)


# renv::restore()
# ------------------------------------------------------------------- #
library(doParallel)
library(tidyverse)
library(parallel)
library(splines)
library(slider)
library(haven)
library(expm)
library(msm)
library(zoo)
library(rsample)
# these 4 can be commented out !!!!!!
library(multidplyr)
library(future)
library(furrr)

# impute age using some simple logic
impute_age <- function(age, wave){
  if (all(is.na(age))){
    return(age)
  }
  age <-  zoo::na.approx(age, x = wave, na.rm = FALSE)
  if (any(is.na(age))){
    if (is.na(age[1])){
      # handle leading NAs
      diffs            <- diff(wave) * 2
      nNAs             <- rle(is.na(age))$lengths[1]
      first_non_NA_age <- age[!is.na(age)][1]
      subtract_this    <- diffs[1:nNAs] |> rev() |> cumsum() |> rev()
      age[1:nNAs]      <- first_non_NA_age - subtract_this
    }
    n <- length(age)
    if (is.na(age[n])){
      diffs            <- diff(wave) * 2
      nNAs             <- rle(is.na(rev(age)))$lengths[1]
      last_non_NA_age  <- age[!is.na(age)] |> rev() %>% '['(1)
      add_this         <- diffs[(n-nNAs):(n-1)] |> cumsum()
      age[(n-nNAs+1):n]  <- last_non_NA_age + add_this
    }
  }
  age
}


# ----------------------------------------------------
# custom bootstrap function built atop rsample package:
# ----------------------------------------------------

# Clone IDs inside a bootstrap sample so each resampled copy becomes its own subject
clone_ids_in_boot_sample <- function(data,
                                     id_col   = "hhidpn",
                                     time_var = "age") {
  id_sym   <- rlang::sym(id_col)
  time_sym <- rlang::sym(time_var)
  
  data %>%
    # within each (id, time) cell, give each row its own "copy index"
    dplyr::group_by(!!id_sym, !!time_sym) %>%
    dplyr::mutate(
      .copy_idx = dplyr::row_number(),
      !!id_sym := paste0(!!id_sym, "_", .copy_idx)
    ) %>%
    dplyr::ungroup() %>%
    # nice for msm: sorted by new id and time
    dplyr::arrange(!!id_sym, !!time_sym)
}


group_bootstraps2 <- function(data,
                              group,          # <- expects a character string!
                              times = 25,
                              apparent = FALSE,
                              ...,
                              strata = NULL,
                              pool = 0.1,
                              weight = NULL) {
  # Validate the group column exists
  if (!group %in% names(data)) {
    stop(glue::glue("The group column '{group}' does not exist in the data."))
  }
  
  # If strata is passed as character, check it exists
  if (!is.null(strata) && !strata %in% names(data)) {
    stop(glue::glue("The strata column '{strata}' does not exist in the data."))
  }
  
  # Use validate_group with explicit all_of() to avoid tidyselect deprecation warning
  group_var <- rsample:::validate_group(tidyselect::all_of(group), data)
  
  # Validate strata if used
  if (!is.null(strata)) {
    strata <- rsample:::check_grouped_strata(tidyselect::all_of(group), strata, pool, data)
  }
  
  # Generate splits
  split_objs <- group_boot_splits2(
    data = data,
    group = group_var,
    times = times,
    strata = strata,
    pool = pool,
    weight = weight
  )
  
  split_objs$splits <- purrr::map(split_objs$splits, rsample:::rm_out)
  
  if (apparent) {
    split_objs <- dplyr::bind_rows(split_objs, rsample::apparent(data))
  }
  
  if (is.null(strata)) {
    strata <- FALSE
  }
  
  boot_att <- list(
    times = times,
    apparent = apparent,
    strata = strata,
    pool = pool,
    group = group_var
  )
  
  rsample::new_rset(
    splits = split_objs$splits,
    ids = split_objs$id,
    attrib = boot_att,
    subclass = c("group_bootstraps", "bootstraps", "group_rset", "rset")
  )
}




group_boot_splits2 <- function(data, group, times = 25, strata = NULL, pool = 0.1, weight = NULL) {
  group_vals <- getElement(data, group)
  
  if (!is.null(strata)) {
    strata_vals <- getElement(data, strata)
    strata_vals <- as.character(strata_vals)
    strata_vals <- rsample:::make_strata(strata_vals, pool = pool)
  } else {
    strata_vals <- NULL
  }
  
  n <- nrow(data)
  
  if (!is.null(weight)) {
    first_rows <- data[!duplicated(group_vals), , drop = FALSE]
    wt <- first_rows[[weight]]
    names(wt) <- as.character(first_rows[[group]])
    wt <- ifelse(is.na(wt), 0, wt)
    wt <- wt / sum(wt)
  } else {
    wt <- NULL
  }
  
  indices <- make_groups2(
    data = data,
    group = group_vals,
    v = times,
    balance = "prop",
    replace = TRUE,
    strata = strata_vals,
    weights = wt
  )
  
  group_vals <- as.character(group_vals)
  group_to_rows <- split(seq_along(group_vals), group_vals)
  
  indices <- lapply(indices, function(fold_df) {
    group_ids <- rep(fold_df$..group, fold_df$..replications)
    matched <- group_to_rows[as.character(group_ids)]
    
    if (any(sapply(matched, is.null))) {
      missing_ids <- group_ids[sapply(matched, is.null)]
      stop(glue::glue("One or more sampled group IDs did not match any rows: {paste(missing_ids, collapse = ', ')}"))
    }
    
    unlist(matched, use.names = FALSE)
  })
  
  split_objs <- purrr::map(indices, function(idx) {
    if (length(idx) == 0 || any(is.na(idx)) || any(idx <= 0)) {
      stop("Invalid analysis index: empty or contains NA / non-positive value.")
    }
    
    rsample::make_splits(
      list(analysis = idx, assessment = integer(0)),
      data = data,
      class = c("group_boot_split", "boot_split")
    )
  })
  
  if (getOption("group_boot.verbose", TRUE)) {
    all_assessable <- purrr::map_int(split_objs, ~nrow(rsample::assessment(.x)))
    if (any(all_assessable == 0)) {
      cli::cli_warn(
        c(
          "Some assessment sets contained zero rows.",
          i = "Consider using a non-grouped resampling method."
        ),
        call = rlang::caller_env()
      )
    }
  }
  
  list(
    splits = split_objs,
    id = rsample:::names0(length(split_objs), "Bootstrap")
  )
}

make_groups2 <- function(data, group, v, balance = c("groups", "observations", "prop"), strata = NULL, ...) {
  balance <- rlang::arg_match(balance, error_call = rlang::caller_env())
  data_ind <- tibble(..index = seq_len(nrow(data)), ..group = group)
  data_ind$..group <- as.character(data_ind$..group)
  
  res <- switch(balance, 
                groups = rsample:::balance_groups(data_ind = data_ind, v = v, strata = strata, ...), 
                observations = rsample:::balance_observations(data_ind = data_ind, v = v, strata = strata, ...),
                prop = {
                  keys = balance_prop2(data_ind = data_ind, v = v, strata = strata, ...)
                  list(data_ind = data_ind, keys = keys)
                }
  )
  
  data_ind <- res$data_ind
  keys <- res$keys
  if (nrow(keys) == 0) {
    stop("No valid bootstrap groups were generated. Check that group column is not empty.")
  }
  
  data_ind$..group <- as.character(data_ind$..group)
  keys$..group <- as.character(keys$..group)
  
  split(keys, keys$..folds)
}


balance_prop2 <- function(data_ind, v, replace = TRUE, strata = NULL, weights = NULL, ...) {
  # Split data_ind by strata if given
  data_splits <- if (is.null(strata)) {
    list(data_ind)
  } else {
    rsample:::split_unnamed(data_ind, strata)
  }
  
  # Sample groups per stratum using helper
  keys <- purrr::map_dfr(
    data_splits,
    balance_prop_helper2,
    v = v,
    replace = replace,
    weights = weights
  )
  
  keys
}

balance_prop_helper2 <- function(data_ind, v, replace = TRUE, weights = NULL) {
  group_ids <- sort(unique(as.character(data_ind$..group)))
  
  if (length(group_ids) == 0) {
    warning("No group IDs found; skipping bootstrap resampling for this stratum.")
    return(tibble(..group = character(0), ..replications = integer(0), ..folds = integer(0)))
  }
  
  # Normalize weights
  if (!is.null(weights)) {
    weights_vec <- weights[group_ids]
    weights_vec[is.na(weights_vec)] <- 0
    if (sum(weights_vec) == 0) {
      warning("All weights are zero or NA for the selected groups; using uniform weights.")
      weights_vec <- rep(1, length(group_ids))
    }
    weights_vec <- weights_vec / sum(weights_vec)
  } else {
    weights_vec <- rep(1, length(group_ids)) / length(group_ids)
  }
  # message("Debug: v = ", v)
  purrr::map_dfr(seq_len(v), function(x) {
    if (length(group_ids) == 0 || anyNA(group_ids)) {
      return(tibble(..group = character(0), ..replications = integer(0), ..folds = integer(0)))
    }
    
    sampled_groups <- sample(
      group_ids,
      size = length(group_ids),
      replace = TRUE,
      prob = weights_vec
    )
    
    rep_counts <- table(sampled_groups)
    tibble(
      ..group = names(rep_counts),
      ..replications = as.integer(rep_counts),
      ..folds = x
    )
  })
}


# ------------------------------------
# fit_msm() functions
# ------------------------

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
  if (length(drop_vars)>0) {
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
      strat_vals <- if (length(strat_vars)>0) group_data[1, strat_vars, drop = FALSE] else NULL
      
      # fit msm for this stratum
      mod <- msm_model(
        data         = group_data,
        Q            = Q,
        covariate_var= covariate_var,
        extra_covars = spline_names
      )
      if (is.null(mod)) {
        # whole stratum fails: return pred_sub with NA matrices
        pred_sub <- pred_grid
        if (length(strat_vars)) {
          for (v in strat_vars) {
            pred_sub <- pred_sub[pred_sub[[v]] == strat_vals[[v]], , drop = FALSE]
          }
        }
        pred_sub$q <- replicate(
          nrow(pred_sub),
          matrix(NA_real_, nrow = nrow(Q), ncol = ncol(Q)),
          simplify = FALSE
        )
        return(pred_sub)
      }
      # subset prediction grid to this stratum
      pred_sub <- pred_grid
      if (length(strat_vars)>0) {
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
        if (length(cov_cols)>0) {
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
    
    if (length(missing_lvls)>0) {
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

# Build a shared bootstrap design (id-level)
make_boot_design <- function(
    data,
    id_col     = "hhidpn",
    times      = 1000,
    weight_col = NULL
) {
  group_bootstraps2(
    data   = data,
    group  = id_col,
    times  = times,
    weight = weight_col
  )
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
# -------------------------------------------------------------------
# Bootstrap wrapper around fit_msm() — id_col & weights first-class,
# with optional shared boot_rset and seed.
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
    parallel          = c("auto","none","future","mclapply"),
    boot_rset         = NULL,  # 
    seed              = NULL   # 
) {
  
  parallel <- match.arg(parallel)
  # if parallel == "future" but the user's current plan is sequential ----
  if (parallel == "future") {
    current_strategy <- class(future::plan())[1]
    
    # future::sequential --> SequentialFuture / multicore disabled
    if (grepl("sequential", tolower(current_strategy))) {
      stop(
        paste0(
          "parallel = 'future' was requested, but the current future plan is sequential.\n",
          "Please run this in your console BEFORE calling fit_msm_boot():\n\n",
          "    future::plan(future::multisession, workers = ", n_cores, ")\n\n",
          "Or choose a different parallel method: parallel = 'none' or 'mclapply'.\n"
        ),
        call. = FALSE
      )
    }
  }
  
  # ---- helper: standardise join cols for safe joins ----
  coerce_join_cols <- function(df) {
    if (!is.null(strat_vars) && length(strat_vars)) {
      df <- df %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(strat_vars), as.character))
    }
    if (!is.null(covariate_var) && length(covariate_var)) {
      df <- df %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(covariate_var), as.character))
    }
    df %>%
      dplyr::mutate(
        age  = round(as.numeric(.data$age), 6),
        from = as.integer(.data$from),
        to   = as.integer(.data$to)
      )
  }
  
  # ---- 0) Master prediction grid so all replicates align ----
  
  master_grid <- build_prediction_grid(
    data          = data,
    strat_vars    = strat_vars,
    covariate_var = covariate_var,
    age_from_to   = age_from_to,
    age_int       = age_int
  ) %>%
    tidyr::crossing(
      from = seq_len(nrow(Q)),
      to   = seq_len(ncol(Q))
    ) %>%
    dplyr::mutate(
      from = as.integer(from),
      to   = as.integer(to)
    ) %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(strat_vars), as.character),
      dplyr::across(dplyr::all_of(covariate_var), as.character),
      age = round(as.numeric(age), 6)
    )
  
  join_keys <- unique(c(strat_vars, covariate_var, "age", "from", "to"))
  
  master_grid <- master_grid %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(join_keys)), .keep_all = TRUE)
  
  master_grid <- coerce_join_cols(master_grid)
  

  
    fit_msm_on_split <- function(split, rep_idx) {
      tryCatch(
        {
          dat_b <- rsample::analysis(split) %>%
            clone_ids_in_boot_sample(id_col = id_col, time_var = "age")
          
          haz_b <- fit_msm(
            data          = dat_b,
            strat_vars    = strat_vars,
            covariate_var = covariate_var,
            age_from_to   = age_from_to,
            age_int       = age_int,
            spline_df     = spline_df,
            spline_type   = spline_type,
            Q             = Q
          )
          
          if (!nrow(haz_b)) {
            out <- master_grid
            out$haz <- NA_real_
          } else {
            haz_b2 <- haz_b %>%
              dplyr::mutate(
                dplyr::across(dplyr::all_of(strat_vars), as.character),
                dplyr::across(dplyr::all_of(covariate_var), as.character),
            age  = round(as.numeric(.data$age), 6),
            from = as.integer(.data$from),
            to   = as.integer(.data$to)
            )
      
      out <- master_grid %>%
        dplyr::left_join(
          haz_b2,
          by     = join_keys,
          suffix = c("", ".rep")
        )
        }
          
          out$replicate <- rep_idx
          out
          },
      error = function(e) {
        # This is where your "numerical overflow" lands
        message(
          "Warning: fit_msm() failed for replicate ", rep_idx,
          " with error: ", conditionMessage(e)
        )
        # Fallback: full grid with NA hazards for this replicate
        out <- master_grid
        out$haz <- NA_real_
        out$replicate <- rep_idx
        out
      }
      )
  }
    
  
  # ---- 1) Stratified, group bootstrap replicates (weighted if provided) ----
  # Use provided boot_rset if given; otherwise build it here.
  
  if (is.null(boot_rset)) {
    if (!is.null(seed)) set.seed(seed)
    boot_rset <- group_bootstraps2(
      data   = data,
      group  = id_col,       # character, as expected by group_bootstraps2
      times  = times,
      weight = weight_col
    )
  }
  
  n_boot <- length(boot_rset$splits)
  
  boot_tbl <- tibble::tibble(
    replicate = seq_len(n_boot),
    split     = boot_rset$splits
  )
  
  # ---- 2) Evaluate across replicates (serial or parallel) ----
  
  eval_serial <- function(boot_tbl) {
    boot_tbl %>%
      dplyr::group_by(replicate) %>%
      dplyr::group_modify(
        ~ fit_msm_on_split(
          split   = .x$split[[1]],
          rep_idx = .y$replicate[[1]]
        )
      ) %>%
      dplyr::ungroup()
  }
  
  eval_future <- function(boot_tbl) {
    furrr::future_map_dfr(
      seq_len(nrow(boot_tbl)),
      function(i) {
        safe_fit_msm_on_split(
          split   = boot_tbl$split[[i]],
          rep_idx = i
        )
      },
      .options = furrr::furrr_options(seed = FALSE)
    )
  }
  
  eval_mclapply <- function(boot_tbl) {
    out_list <- parallel::mclapply(
      X = seq_len(nrow(boot_tbl)),
      FUN = function(i) {
        fit_msm_on_split(
          split   = boot_tbl$split[[i]],
          rep_idx = i
        )
      },
      mc.cores = n_cores
    )
    dplyr::bind_rows(out_list)
  }
  
  
  if (parallel == "auto") {
    if (n_cores <= 1) {
      parallel <- "none"
    } else {
      parallel <- "future"
    }
  }
  
  boot_haz <- switch(
    parallel,
    "none"     = eval_serial(boot_tbl),
    "future"   = eval_future(boot_tbl),
    "mclapply" = eval_mclapply(boot_tbl),
    # fallback
    eval_serial(boot_tbl)
  )
  
  boot_haz <- boot_haz %>%
    coerce_join_cols() %>%
    dplyr::mutate(
      replicate = as.integer(.data$replicate),
      haz       = as.numeric(.data$haz)
    ) %>%
    dplyr::select(-dplyr::any_of(".groups"))
  
  # ---- 3) Summarise across replicates ----
  
  alpha  <- 1 - ci_level
  q_low  <- alpha / 2
  q_high <- 1 - alpha / 2
  
  haz_summary <- boot_haz %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(join_keys))) %>%
    dplyr::summarise(
      haz_mean = mean(haz, na.rm = TRUE),
      haz_low  = as.numeric(stats::quantile(
        haz, probs = q_low,  na.rm = TRUE,
        names = FALSE, type = 6
      )),
      haz_high = as.numeric(stats::quantile(
        haz, probs = q_high, na.rm = TRUE,
        names = FALSE, type = 6
      )),
      .groups  = "drop"
    )
  
  if (return_replicates) {
    return(list(
      summary    = haz_summary,
      replicates = boot_haz
    ))
  } else {
    return(haz_summary)
  }
}


# -----------------------------------------------------
# prevalence functions
# -----------------------------------------------------


# 14-11-2025 The below is still untested, as I've got hazard bootstrapping 
# running in the background. To continue in a later sitting.
#' Fit age-specific prevalence with GAM
#' Fit age-specific prevalence with GAM
#'
#' @param data          Long-format individual data (one bootstrap sample or full data).
#' @param strat_vars    Character vector of stratification vars (e.g. "female").
#' @param covariate_var Character vector of covariates (e.g. "period").
#' @param age_from_to   Numeric length-2, lower/upper bounds of age grid.
#' @param age_int       Numeric age step for prediction grid.
#' @param state_var     Name of state variable (default "state_msm").
#' @param id_col        Subject ID column (default "hhidpn").
#' @param condition_state Integer code of the "condition" state (default 2).
#' @param spline_df     Optional degrees of freedom for age spline (NULL = linear age).
#' @param spline_type   Type of spline for age ("ns" or "bs"), as in fit_msm().
#' @param weight_col           Optional name of person-weight column. If NULL, uses equal weights.
#'
#' @return Tibble with strat_vars, covariate_var, age, prevalence.
fit_prev <- function(
    data,
    strat_vars      = NULL,
    covariate_var   = NULL,
    age_from_to     = c(50, 100),
    age_int         = 0.25,
    state_var       = "state_msm",
    id_col          = "hhidpn",
    condition_state = 2L,
    spline_df       = NULL,
    spline_type     = "ns",
    weight_col      = NULL,
    exclude_state   = 3L
) {
  # --- 0) sanity / NA handling ---------------------------------------------

  # drop rows with NAs
  drop_vars <- unique(c(strat_vars, covariate_var, "age", state_var, id_col, weight_col))
  if (length(drop_vars)>0) {
    data <- data %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(drop_vars), ~ !is.na(.x)))
  }
  
  # drop excluded states (e.g. dead) BEFORE any grouping
  if (!is.null(exclude_state)) {
    data <- data %>%
      dplyr::filter(!.data[[state_var]] %in% exclude_state)
  }
  state_sym <- rlang::sym(state_var)
  id_sym    <- rlang::sym(id_col)
  
  # keep only subjects with >1 row (as in student's code)
  dat_use <- data %>%
    dplyr::group_by(!!id_sym) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()
  
  # define weights (person weights if supplied, else 1)
  if (is.null(weight_col)) {
    dat_use <- dat_use %>%
      dplyr::mutate(.w = 1)
  } else {
    w_sym <- rlang::sym(weight_col)
    dat_use <- dat_use %>%
      dplyr::mutate(.w = !!w_sym)
  }
  
  # --- 1) empirical (weighted) prevalence by whole year --------------------
  group_vars  <- c(strat_vars, covariate_var, "age")
  
  prev_counts <- dat_use %>%
    dplyr::mutate(
      age  = floor(.data$age),
      .cond = (!!state_sym) == condition_state
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      N = sum(.w * .cond, na.rm = TRUE),   # weighted "cases"
      n = sum(.w,           na.rm = TRUE), # weighted total
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      prev = dplyr::if_else(n > 0, N / n, NA_real_)
    ) %>%
    dplyr::filter(dplyr::between(age, age_from_to[1], age_from_to[2]))
  
  if (!nrow(prev_counts)) {
    # nothing to fit
    return(tibble::tibble())
  }
  
  # --- 2) prediction grid (same helper as hazard code) ---------------------
  new_data <- build_prediction_grid(
    data          = prev_counts,
    strat_vars    = strat_vars,
    covariate_var = covariate_var,
    age_from_to   = age_from_to,
    age_int       = age_int
  )
  
  # --- 3) GAM formula with optional spline on age --------------------------
  # construct age term (linear or spline)
  age_term <- if (is.null(spline_df)) {
    "age"
  } else {
    if (identical(spline_type, "ns")) {
      sprintf("splines::ns(age, df = %d)", as.integer(spline_df))
    } else if (identical(spline_type, "bs")) {
      sprintf("splines::bs(age, df = %d)", as.integer(spline_df))
    } else {
      stop("spline_type must be 'ns' or 'bs'")
    }
  }
  
  cov_terms <- covariate_var
  cov_terms <- cov_terms[!is.na(cov_terms) & nzchar(cov_terms)]
  
  rhs_terms <- c(age_term, cov_terms)
  rhs_str   <- paste(rhs_terms, collapse = " + ")
  
  form <- stats::as.formula(
    paste("cbind(N, n - N) ~", rhs_str)
  )
  
  # --- 4) fit per stratum & predict ---------------------------------------
  
  split_prev <- if (length(strat_vars)>0) {
    prev_counts %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(strat_vars))) %>%
      dplyr::group_split()
  } else {
    list(prev_counts)
  }
  
  out_list <- purrr::map(split_prev, function(df_stratum) {
    # current stratum values
    strat_vals <- if (length(strat_vars)>0) df_stratum[1, strat_vars, drop = FALSE] else NULL
    
    # fit gam for this stratum on weighted counts
    mod <- mgcv::gam(
      formula = form,
      family  = stats::binomial(link = "logit"),
      data    = df_stratum
    )
    
    # prediction grid subset for this stratum
    nd <- new_data
    if (length(strat_vars)>0) {
      for (v in strat_vars) {
        nd <- nd[nd[[v]] == strat_vals[[v]], , drop = FALSE]
      }
    }
    
    if (!nrow(nd)) {
      return(NULL)
    }
    
    nd$prevalence <- stats::predict(mod, newdata = nd, type = "response")
    nd
  })
  
  out <- out_list %>%
    purrr::compact() %>%
    dplyr::bind_rows()
  
  out %>%
    dplyr::select(dplyr::all_of(c(strat_vars, covariate_var, "age", "prevalence")))
}



#' Fit age-specific prevalence with GAM
#'
#' @param data          Long-format individual data (one bootstrap sample or full data).
#' @param strat_vars    Character vector of stratification vars (e.g. "female").
#' @param covariate_var Character vector of covariates (e.g. "period").
#' @param age_from_to   Numeric length-2, lower/upper bounds of age grid.
#' @param age_int       Numeric age step for prediction grid.
#' @param state_var     Name of state variable (default "state_msm").
#' @param id_col        Subject ID column (default "hhidpn").
#' @param condition_state Integer code of the "condition" state (default 2).
#' @param spline_df     Optional degrees of freedom for age spline (NULL = linear age).
#' @param spline_type   Type of spline for age ("ns" or "bs"), as in fit_msm().
#' @param weight_col    Optional name of person-weight column. If NULL, uses equal weights.
#' @param boot_rset optional external design created by group_bootstraps2()
#' @param seed integer optional seed
#' @param exclude_state remove this state before calculating prevalence
fit_prev_boot <- function(
    data,
    id_col            = "hhidpn",
    strat_vars        = NULL,
    covariate_var     = NULL,
    age_from_to       = c(50, 100),
    age_int           = 0.25,
    state_var         = "state_msm",
    condition_state   = 2L,
    spline_df         = NULL,
    spline_type       = "ns",
    weight_col        = NULL,
    times             = 8,
    ci_level          = 0.95,
    return_replicates = FALSE,
    boot_rset         = NULL,  # << new: optional external design
    seed              = NULL,  # << new: optional seed
    exclude_state     = 3L     # << pass through to fit_prev()
) {
  
  # ---------------- helper: standardise join cols ---------------------
  
  coerce_join_cols <- function(df) {
    if (!is.null(strat_vars) && length(strat_vars)) {
      df <- df %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(strat_vars), as.character))
    }
    if (!is.null(covariate_var) && length(covariate_var)) {
      df <- df %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(covariate_var), as.character))
    }
    df %>%
      dplyr::mutate(age = round(as.numeric(.data$age), 6))
  }
  
  # ---------------- 0) master prediction grid -------------------------
  
  master_grid <- build_prediction_grid(
    data          = data,
    strat_vars    = strat_vars,
    covariate_var = covariate_var,
    age_from_to   = age_from_to,
    age_int       = age_int
  )
  
  master_grid <- coerce_join_cols(master_grid)
  join_keys   <- unique(c(strat_vars, covariate_var, "age"))
  
  # ---------------- helper: fit on one bootstrap split ----------------
  
  fit_prev_on_split <- function(split) {
    
    dat_b <- rsample::analysis(split)
    
    prev_b <- fit_prev(
      data            = dat_b,
      strat_vars      = strat_vars,
      covariate_var   = covariate_var,
      age_from_to     = age_from_to,
      age_int         = age_int,
      state_var       = state_var,
      id_col          = id_col,
      condition_state = condition_state,
      spline_df       = spline_df,
      spline_type     = spline_type,
      weight_col      = weight_col,
      exclude_state   = exclude_state  # ensure we drop state 3 here
    )
    
    template <- master_grid
    
    if (!nrow(prev_b)) {
      template$prevalence <- NA_real_
      return(template)
    }
    
    prev_b2 <- coerce_join_cols(prev_b)
    
    out <- template %>%
      dplyr::left_join(
        prev_b2,
        by     = join_keys,
        suffix = c("", ".rep")
      )
    
    out
  }
  
  # ---------------- 1) bootstrap replicates (grouped by id) -----------
  
  # use provided design or build our own
  if (is.null(boot_rset)) {
    if (!is.null(seed)) set.seed(seed)
    boot_rset <- group_bootstraps2(
      data   = data,
      group  = id_col,
      times  = times,
      weight = weight_col
    )
  }
  
  # enforce replicate IDs as 1:n_boot
  n_boot   <- length(boot_rset$splits)
  boot_tbl <- tibble::tibble(
    replicate = seq_len(n_boot),
    split     = boot_rset$splits
  )
  
  # ---------------- 2) serial evaluation over replicates --------------
  
  boot_prev <- boot_tbl %>%
    dplyr::group_by(replicate) %>%
    dplyr::group_modify(
      ~ fit_prev_on_split(.x$split[[1]])
    ) %>%
    dplyr::ungroup()
  
  boot_prev <- boot_prev %>%
    coerce_join_cols() %>%
    dplyr::mutate(
      replicate  = as.integer(.data$replicate),
      prevalence = as.numeric(.data$prevalence)
    ) %>%
    dplyr::select(-dplyr::any_of(".groups"))
  
  # ---------------- 3) summarise over replicates ----------------------
  
  alpha  <- 1 - ci_level
  q_low  <- alpha / 2
  q_high <- 1 - alpha / 2
  
  prev_summary <- boot_prev %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(join_keys))) %>%
    dplyr::summarise(
      prev_mean = mean(prevalence, na.rm = TRUE),
      prev_low  = as.numeric(stats::quantile(
        prevalence, probs = q_low,  na.rm = TRUE,
        names = FALSE, type = 6
      )),
      prev_high = as.numeric(stats::quantile(
        prevalence, probs = q_high, na.rm = TRUE,
        names = FALSE, type = 6
      )),
      .groups   = "drop"
    )
  
  if (return_replicates) {
    return(list(
      summary    = prev_summary,
      replicates = boot_prev
    ))
  } else {
    return(prev_summary)
  }
}









