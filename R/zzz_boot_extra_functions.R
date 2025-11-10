fit_msm <- function(data,
                    strat_vars = NULL,
                    covariate_var = NULL,
                    age_from_to = c(50, 100),
                    age_int = 0.25,
                    spline_type = "ns",
                    spline_df = NULL,
                    Q = default_q) {
  
  # Step 1: Add splines (if any)
  if (!is.null(spline_df)) {
    splines_out      <- add_splines(data)
    data             <- splines_out$data
    spline_basis_fit <- splines_out$basis_fit
  }
  
  # Step 2: Build prediction grid
  prediction_grid <- build_prediction_grid(data, strat_vars, covariate_var, age_from_to, age_int)
  if (!is.null(spline_df)) {
    prediction_grid <- add_splines_to_grid(prediction_grid, spline_basis_fit)
  }
  
  # Step 3: Group data and compute Q matrices directly
  result <- data %>%
    group_by(across(all_of(strat_vars))) %>%
    group_split() %>%
    purrr::map_dfr(function(group_data) {
      strat_vals <- group_data[1, strat_vars, drop = FALSE]
      mod <- msm_model(group_data, Q, covariate_var)
      
      # Subset prediction grid rows that match current stratification
      match_rows <- rep(TRUE, nrow(prediction_grid))
      for (var in strat_vars) {
        match_rows <- match_rows & (prediction_grid[[var]] == strat_vals[[var]])
      }
      pred_sub <- prediction_grid[match_rows, ]
      
      # Compute Q matrices
      pred_sub <- pred_sub %>%
        mutate(q = map(data_fit, function(covars) {
          covs <- safe_covariate_list(covars, covariate_var)
          if (anyNA(covs)) return(matrix(NA, nrow = nrow(Q), ncol = ncol(Q)))
          qmat <- tryCatch(qmatrix.msm(mod, covariates = covs, ci = "none"),
                           error = function(e) matrix(NA, nrow = nrow(Q), ncol = ncol(Q)))
          if (is.matrix(qmat)) qmat else qmat$estimate
        }))
      
      # Drop bulky nested column
      select(pred_sub, -data_fit)
    })
  
  # Step 4: Tidy final output
  output <- result %>%
    mutate(haz = map(q, tidy_q)) %>%
    unnest(haz) %>%
    select(all_of(c(strat_vars, covariate_var, "age", "from", "to", "estimate"))) %>%
    rename(haz = estimate) %>%
    mutate(
      across(any_of(c(strat_vars, covariate_var)), ~ if (is.character(.x)) factor(.x) else .x),
      from = readr::parse_number(as.character(from)) |> as.integer(),
      to   = readr::parse_number(as.character(to)) |> as.integer()
    )
  
  return(output)
}



# Bootstrap wrapper
fit_msm_boot <- function(data,
                         strat_vars = NULL,
                         covariate_var = NULL,
                         times = 8,
                         weight = "pwt",
                         group = "hhidpn",
                         n_cores = 4,
                         Q,
                         age_from_to = c(50, 100),
                         age_int = 0.25,
                         spline_type = "ns",
                         spline_df = NULL) {
  
  boot_rset <- group_bootstraps2(
    data = data,
    group = group,
    weight = weight,
    times = times
  )
  
  if (Sys.getenv("RSTUDIO") != 1 & .Platform$OS.type == "unix"){
    future::plan("multicore", workers = n_cores)
  } else {
    future::plan("multisession", workers = n_cores)
  }
 

  cat("🚀 Running in parallel with", future::nbrOfWorkers(), "worker(s)\n")
  
  results <- furrr::future_map_dfr(
    boot_rset$splits,
    ~ fit_msm(
      data = rsample::analysis(.x) |> arrange(hhidpn, age),
      strat_vars = strat_vars,
      covariate_var = covariate_var,
      age_from_to = age_from_to,
      age_int = age_int,
      spline_type = spline_type,
      spline_df = spline_df,
      Q = Q
    ),
    .options = furrr::furrr_options(
      seed = TRUE,
      packages = c("dplyr", "tidyr", "msm", "readr"),
      globals = list(
        fit_msm = fit_msm,
        msm_model = msm_model,
        tidy_q = tidy_q,
        build_prediction_grid = build_prediction_grid,
        add_splines = add_splines,
        add_splines_to_grid = add_splines_to_grid,
        safe_covariate_list = safe_covariate_list,
        crossing = tidyr::crossing,
        Q = Q,
        age_from_to = age_from_to,
        age_int = age_int,
        spline_type = spline_type,
        spline_df = spline_df
      )
    ),
    .progress = TRUE,
    .id = "replicate",
    .keep = "none"
  )
  
  rm(boot_rset)
  future::plan("sequential")  # Clean up workers
  gc()
  
  ci_result <- results %>%
    group_by(across(-replicate)) %>%
    summarise(
      haz = mean(haz, na.rm = TRUE),
      haz_l = quantile(haz, 0.025, na.rm = TRUE),
      haz_u = quantile(haz, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(ci_result)
}


# Build the prediction grid from strata/covariates
build_prediction_grid <- function(data, strat_vars, covariate_var, age_from_to, age_int) {
  base_grid <- list(age = seq(age_from_to[1], age_from_to[2], by = age_int))
  
  
  if (!is.null(covariate_var)) {
    cov_vals <- map(covariate_var, ~ unique(data[[.x]]))
    names(cov_vals) <- covariate_var
  } else {
    cov_vals <- list()
  }
  
  if (!is.null(strat_vars)) {
    strat_vals <- map(strat_vars, ~ unique(data[[.x]]))
    names(strat_vals) <- strat_vars
  } else {
    strat_vals <- list()
  }
  
  full_grid <- exec(crossing, !!!c(base_grid, cov_vals, strat_vals))
  full_grid <- full_grid %>%
    group_by(across(all_of(c(strat_vars, covariate_var, "age")))) %>%
    group_nest(.key = "data_fit") %>%
    mutate(data_fit = map(data_fit, ~ as.list(.x[1, , drop = FALSE])))
  return(full_grid)
}

# Add spline basis to data
add_splines <- function(data, age_var = "age", spline_type = "ns", spline_df = 4) {
  age_values <- data[[age_var]]
  
  # Build spline basis once using the full age range from the data
  basis_constructor <- match.fun(spline_type)
  spline_basis_fit <- basis_constructor(age_values, df = spline_df)
  
  # Apply fitted basis to data
  spline_cols <- predict(spline_basis_fit, newx = age_values) %>%
    as_tibble() %>%
    set_names(paste0("age_spline", seq_len(ncol(.))))
  
  # Return the updated data and the basis fit for reuse
  data_out <- bind_cols(data, spline_cols)
  list(data = data_out, basis_fit = spline_basis_fit)
}

add_splines_to_grid <- function(pred_grid, basis_fit) {
  # Use same spline fit to ensure basis alignment
  spline_cols <- predict(basis_fit, newx = pred_grid$age) %>%
    as_tibble() %>%
    set_names(paste0("age_spline", seq_len(ncol(.))))
  bind_cols(pred_grid, spline_cols)
}

# MSM model fitting helper
msm_model <- function(data, Q, covariate_var) {
  suppressWarnings(
    msm(
      formula = state_msm ~ age,
      subject = hhidpn,
      data = data,
      qmatrix = Q,
      obstype = obstype,
      deathexact = 3,
      covariates = reformulate(covariate_var),
      control = list(fnscale = 5000, maxit = 25000),
      gen.inits = TRUE,
      method = "BFGS"
    )
  )
}

# Tidy Q matrix to long form
tidy_q <- function(qmat) {
  if (is.null(qmat)) return(tibble(from = NA, to = NA, estimate = NA))
  if (!is.matrix(qmat)) stop("Expected a matrix in `tidy_q()` but got: ", class(qmat))
  as.data.frame(as.table(qmat)) %>%
    rename(from = Var1, to = Var2, estimate = Freq)
}

safe_covariate_list <- function(x, covariate_var = NULL) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) {
    x_out <- as.list(x[1, , drop = FALSE])
    if (!is.null(covariate_var)) {
      x_out <- x_out[covariate_var]
    }
    return(x_out)
  }
  if (is.list(x)) return(x)
  stop("Unexpected type for covariate input")
}
