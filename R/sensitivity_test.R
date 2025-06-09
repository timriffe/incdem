# ------------------------------------------------------------------- #
# load data
source("R/01b_data_preparation.R")
# ------------------------------------------------------------------- #
# Q tells msm what the valid transitions are
Q <- rbind(
  c(0, 0.1, 0.1),  # healthy can go to dementia or death
  c(0, 0,   0.1),  # dementia can go to death
  c(0, 0,   0)     # death is absorbing
)
diag(Q) <- -rowSums(Q)
# ------------------------------------------------------------------- #
#' Fit Multistate Model with Sensitivity to Covariates
#'
#' Fits a multistate model by including one or more covariates and generating predicted transition rates and probabilities.
#' The function nests the data by stratification and covariate levels, fits the MSM model for each group, 
#' then produces a prediction grid for each combination of covariates and age, and computes transition rates and probabilities.
#'
#' @param .data A dataframe with individual data. in our case hrs_msm
#' @param strat_vars A character vector of variable names to stratify the MSM (group and nest)
#' @param covariate_var A character vector of covariate variable names to include in the model.
#' @param age_int A numeric value specifying the interval for prediction ages (default is 0.25 years).
#' @param cont_grid This parameter is only used if one of the stratification or vovariate vars is continuous. Say we have continuous time, in this case this variable should take the values like 2000, 2010, for data fitting
#' @return A data frame with estimated transition rates and probabilities 
#' 

# ------------------------------------------------------------------- #
# check the comments and example code below
# this function has 4 (four) parts:
# 1 creating the grid, 
# 2 fitting the msm
# 3 predicting from the created grid and model
# 4 reshaping and returning the outcome
# can take any variables as covariate_var or strat_vars
# if covariate or strat var is continuous cont_grid should be provided
# it should be something ordinal or nominal like c(2000, 2010) 
# a assume that age is always provided 
# and 3 splines exist as covariates by default for now 
# NOTE: Dont pay attention to warning. IT DOES NOT AFFECT RESULTS
# I suggest using newly created obs_date as a time variable
fit_msm_sensitivity <- function(.data,
                                strat_vars    = NULL,
                                covariate_var = NULL,
                                age_int       = 0.25,
                                cont_grid     = NULL) {
  # --------------------------------------------------------------#
  # PT1
  # We start by constructing the grid
  # the ifelse part takes arbitrary specified columns and 
  # makes a list with the correspodning unique values
  # This first part is only used in case we have continuous covariate 
  if(!is.null(cont_grid)) {
    
    vars        <- map(c(strat_vars), 
                       ~ unique(na.omit(.data[[.x]])))
    
    vars[[(length(vars) + 1)]] <- cont_grid 
    
    names(vars) <- c(strat_vars, covariate_var)
    
    # othervise we use this part
  } else { 
    
    vars        <- map(c(covariate_var, strat_vars), 
                       ~ unique(na.omit(.data[[.x]]))) |>
      set_names(c(covariate_var, strat_vars))
    
  }
  
  # add base predictors. in our model it is age
  base_grid <- list(
    age    = seq(50, 100, by = age_int)
  )
  
  # construct prediction grid dynamically
  prediction_grid <- exec(crossing, !!!c(base_grid, vars))
  
  # predict the spline basis for the age specification
  spline_basis <- predict(spline_basis_fit,
                          newx = prediction_grid$age) |>
    as.data.frame() |>
    setNames(paste0("age_spline", 1:3))
  
  # bind data together
  prediction_grid <- bind_cols(prediction_grid, spline_basis) |>
    # mutate(female = ifelse(female == 1, "Female", "Male")) |>
    # nest data by correspodnding covariates
    group_by(across(all_of(c(
      strat_vars, covariate_var, "age"
    )))) |>
    group_nest(keep = TRUE, .key = "data_fit") |>
    # create a nested list of fitting values
    mutate(data_fit = map(data_fit, ~ .x |>
                            dplyr::select(-c(
                              all_of(strat_vars), age
                            )) |>
                            as.list()))
  
  # --------------------------------------------------------------#
  # PT2
  # create fit model with arbitrary covariates
  # always assume that age aplines are used and age is provided
  # NOTE: Takes time
  base_covs   <- c("age_spline1", "age_spline2", "age_spline3")
  all_covs    <- c(base_covs, covariate_var)
  cov_formula <- reformulate(all_covs)
  
  # Stratify by one or more variables
  result <- .data |>
    group_by(across(all_of(strat_vars))) |>
    group_nest() |>
    na.omit() |>
    mutate(model = map(
      data,
      ~ msm(
        formula    = state_msm ~ age,
        subject    = hhidpn,
        data       = .x,
        qmatrix    = Q,
        obstype    = obstype,
        deathexact = 3,
        covariates = cov_formula,
        control    = list(fnscale = 5000, maxit = 25000),
        gen.inits  = TRUE,
        method     = "BFGS"
      )
    )) |>
    dplyr::select(-data)
  
  # --------------------------------------------------------------#
  # PT3
  # here we predict the model by the corresponding created grid and 
  # calculate the rate and prob. matrices
  pred <- prediction_grid |>
    # add our age interval
    mutate(age_interval = age_int) |>
    # join the model
    left_join(result, by = strat_vars) |>
    # calculate rate
    mutate(
      rate = map2(
        .x = model,
        .y = data_fit,
        ~ qmatrix.msm(x          = .x, 
                      covariates = .y)$estimates |>
          as.table() |>
          as.data.frame() |>
          rename(from = Var1, 
                 to   = Var2, 
                 rate = Freq)
      ),
      # calculate probability
      prob = pmap(
        list(model, data_fit, age_interval),
        ~ pmatrix.msm(x          = ..1, 
                      t          = ..3, 
                      covariates = ..2) |>
          as.table() |>
          as.data.frame() |>
          rename(from = Var1, 
                 to   = Var2, 
                 prob = Freq)
      )
    )
  # ------------------------------------------------------------------- #
  # PT4 
  # return resulting data with rates and probabilities 
  # unnest rate and chosen variables
  rate <- pred |>
    dplyr::select(all_of(c(strat_vars, covariate_var)), age, rate) |>
    unnest(rate)
  
  # unnest probability
  prob <- pred |>
    dplyr::select(all_of(c(strat_vars, covariate_var)), age, prob) |>
    unnest(prob)
  
  # this is the final result
  result_df <- rate |>
    # join together
    full_join(prob) |>
    # remove D-D transitions (empty)
    filter(from != "State 3") |>
    # remove recovery possibility (empty)
    filter(!(from == "State 2" & to == "State 1"))
  
  return(result_df)
  
}

# ------------------------------------------------------------------- #
# create obs_date code by TR
hrs_msm <- hrs_msm |>
  # a test condition to check pandemic effect on mort trend
  mutate(
    birth_date         = as_date(birth_date, origin = "1960-01-01"),
    birth_date_decimal = decimal_date(birth_date),
    obs_date           = birth_date_decimal + age
  ) |>
  filter(obs_date < (as_date("2019-jan-01") |> decimal_date())) |>
  # factor variables
  mutate(across(c(female, race, hispanic, education,
                  hypertension, diabetes, heart_disease,
                  stroke, ever_dementia), ~ as.factor(.)))

# ------------------------------------------------------------------- #
# lets take a look at some results
# example 1: stratify by sex and use obs_date as a covariate
# NOTE: since covariate is continuous, I provide cont_grid
# that is used in the prediction part to predict data on these particular point
# in this example it will return predictions for 2000 2010 and 2020
results1 <- hrs_msm |>
  fit_msm_sensitivity(
    strat_vars    = c("female"),
    covariate_var = c("obs_date"),
    age_int       = 0.25,
    cont_grid     = c(2000, 2010, 2020)
  )

# here is the corresponding plot for example 1
results1 |>
  mutate(
    to = to |> as.character() |> parse_number(),
    from = from |> as.character() |> parse_number()
  ) |>
  filter(to > from) |>
  mutate(Period = as.factor(obs_date),
         transition = paste0("haz", from, to)) |>
  ggplot(aes(
    x        = age,
    y        = rate,
    color    = transition,
    linetype = Period
  )) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    x        = "Age",
    y        = "Hazard Rate",
    color    = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap( ~ female) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------- #
# Example 2
# stratify by sex and use stroke and education as covariate
# NOTE: since there is no continuous covariates  cont_grid is set to NULL
results2 <- hrs_msm |>
  # filter(!is.na(Date), !is.na(hypertension)) |>
  fit_msm_sensitivity(
    strat_vars    = c("female"),
    covariate_var = c("education", "stroke"),
    age_int       = 0.25,
    cont_grid     = NULL # !!!
  )

# here is the figure for example 2
results2 |>
  mutate(
    to = to |> as.character() |> parse_number(),
    from = from |> as.character() |> parse_number()
  ) |>
  filter(to > from) |>
  mutate(transition = paste0("haz", from, to)) |>
  ggplot(aes(
    x        = age,
    y        = rate,
    color    = transition,
    linetype = stroke
  )) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    x        = "Age",
    y        = "Hazard Rate",
    color    = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap(education ~ female) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------- #
# Example 3
# Lets separate data by time in 2 periods
hrs_msm <- hrs_msm |>
  # before and after 2010 just as an example
  mutate(Period = ifelse(obs_date < 2010, "<2010", "2011+")) |>
  # factor it
  mutate(Period = as.factor(Period))

# Here we use stroke as covariate and sex and Period as stratas
results3 <- hrs_msm |>
  # filter(!is.na(Date), !is.na(hypertension)) |>
  fit_msm_sensitivity(
    strat_vars    = c("female", "Period"),
    covariate_var = c("stroke"),
    age_int       = 0.25,
    cont_grid     = NULL
  )

# here is the plot for example 3
results3 |>
  mutate(to = to |> as.character() |> parse_number(),
         from = from |> as.character()|> parse_number()) |> 
  filter(to > from) |> 
  mutate(transition = paste0("haz", from, to)) |>
  ggplot(aes(x        = age, 
             y        = rate, 
             color    = transition, 
             linetype = Period)) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    x        = "Age",
    y        = "Hazard Rate",
    color    = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap(stroke ~ female) + 
  theme(legend.position = "bottom")

# ------------------------------------------------------------------- #
# Example 4
# same as 3 but now Period is in covariates
results4 <- hrs_msm |>
  # filter(!is.na(Date), !is.na(hypertension)) |>
  fit_msm_sensitivity(
    strat_vars    = c("female"),
    covariate_var = c("stroke", "Period"),
    age_int       = 0.25,
    cont_grid     = NULL
  )

# here is the plot for example 4
results4 |>
  mutate(to = to |> as.character() |> parse_number(),
         from = from |> as.character()|> parse_number()) |> 
  filter(to > from) |> 
  mutate(transition = paste0("haz", from, to)) |>
  ggplot(aes(x        = age, 
             y        = rate, 
             color    = transition, 
             linetype = Period)) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    x        = "Age",
    y        = "Hazard Rate",
    color    = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap(stroke ~ female) + 
  theme(legend.position = "bottom")

# ------------------------------------------------------------------- #
# end