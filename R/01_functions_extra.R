# Build P at a given age for a single group, using standardised cols
build_P_at_age <- function(probs_grp) {
  ages <- sort(unique(probs_grp$age))
  if (length(ages) != 1L) {
    stop("build_P_at_age() expects a single age in probs_grp.")
  }
  
  states <- sort(unique(c(probs_grp$from, probs_grp$to)))
  n      <- length(states)
  
  P <- matrix(0, nrow = n, ncol = n,
              dimnames = list(states, states))
  
  i_from <- match(probs_grp$from, states)
  i_to   <- match(probs_grp$to,   states)
  P[cbind(i_from, i_to)] <- probs_grp$p
  
  P
}

# Absorbing state = row that sends out (almost) nothing but receives something
detect_absorbing_state <- function(P, tol = 1e-12) {
  n <- nrow(P)
  out_to_others   <- numeric(n)
  in_from_others  <- numeric(n)
  
  for (j in seq_len(n)) {
    out_to_others[j]  <- sum(P[j, -j, drop = FALSE])
    in_from_others[j] <- sum(P[-j, j, drop = FALSE])
  }
  
  idx <- which(out_to_others < tol & in_from_others > tol)
  if (length(idx) != 1L) {
    stop("Expected exactly one absorbing state, found ", length(idx), ".")
  }
  
  rownames(P)[idx]
}

init_constant_stationary <- function(probs_grp, from_age) {
  # probs_grp is already standardised (age, from, to, p) and single group
  P_full <- build_P_at_age(probs_grp[probs_grp$age == from_age, , drop = FALSE])
  
  abs_state <- detect_absorbing_state(P_full)
  living    <- setdiff(rownames(P_full), abs_state)
  
  if (length(living) < 2L) {
    stop("Need at least two non-absorbing states for constant_stationary init.")
  }
  
  U <- P_full[living, living, drop = FALSE]
  
  ev  <- eigen(t(U))
  idx <- which.min(Mod(ev$values - 1))
  v   <- Re(ev$vectors[, idx])
  
  v <- pmax(v, 0)
  if (sum(v) <= 0) {
    stop("Computed stationary eigenvector has non-positive sum.")
  }
  v <- v / sum(v)
  names(v) <- living
  
  max_p <- max(v)
  if (max_p > 0.99) {
    warning(
      "init_constant_stationary(): stationary distribution is highly concentrated (max = ",
      round(max_p, 3),
      "). In one-way/no-return models this typically converges to everyone in ",
      "the worst state. Consider `init_method = 'rev_stationary'` or explicit `init`."
    )
  }
  
  v
}

# Simple Markov backward projection helpers (state-generic)

rpi2u <- function(rpivec,
                  from_state,
                  to_state,
                  start_age,
                  age_interval) {
  
  n   <- length(rpivec)
  out <- cbind(0, rbind(diag(rpivec), 0))
  
  ages       <- ((0:n) * age_interval) + start_age - age_interval
  from_names <- paste(from_state, ages, sep = "::")
  to_names   <- paste(to_state,   ages, sep = "::")
  
  dimnames(out) <- list(to_names, from_names)
  out
}

proj_pop_back <- function(rU, origin_state_index, n_age) {
  # rU is (3*n_age) x (3*n_age)
  N <- nrow(rU)
  if (N != 3 * n_age) {
    stop("proj_pop_back(): dimension mismatch; expected 3*n_age rows.")
  }
  
  ones  <- rep(1, n_age)
  zeros <- rep(0, n_age)
  
  init <- c(zeros, zeros, zeros)
  idx_start <- (origin_state_index - 1L) * n_age + 1L
  idx_end   <- origin_state_index * n_age
  init[idx_start:idx_end] <- ones
  
  pop <- matrix(0, nrow = N, ncol = n_age + 1L)
  pop[, 1L] <- init
  
  for (i in seq_len(n_age)) {
    pop[, i + 1L] <- rU %*% pop[, i]
  }
  
  pop
}

proj_prev_back <- function(rU, origin_state_index, n_age, healthy_block_index = 1L) {
  pop <- proj_pop_back(rU, origin_state_index, n_age)
  
  row_start <- (healthy_block_index - 1L) * n_age + 1L
  row_end   <- healthy_block_index * n_age
  
  PH <- pop[row_start:row_end, 1:n_age, drop = FALSE]
  
  # assume exactly one other living state shares the living mass;
  # here we'll only use the *proportion* in "healthy" within living.
  # With 2 living states, this is PH / (PH + PU); we reconstruct denominator:
  living_rows <- seq_len(2L * n_age)   # two living blocks stacked
  living_pop  <- pop[living_rows, 1:n_age, drop = FALSE]
  living_tot  <- colSums(living_pop)
  
  pi <- PH / (matrix(living_tot, nrow = n_age, ncol = n_age, byrow = TRUE))
  
  pi
}

# Detect states (2 living, 1 absorbing) from standardised probs
detect_three_state_structure <- function(probs_grp) {
  states <- sort(unique(c(probs_grp$from, probs_grp$to)))
  
  if (length(states) != 3L) {
    stop("rev_stationary currently assumes exactly 3 states (2 living + 1 absorbing).")
  }
  
  ages <- sort(unique(probs_grp$age))
  age0 <- ages[1L]
  
  P0 <- build_P_at_age(probs_grp[probs_grp$age == age0, , drop = FALSE])
  abs_state <- detect_absorbing_state(P0)
  living    <- setdiff(states, abs_state)
  
  list(
    living  = living,
    dead    = abs_state,
    states  = states
  )
}
build_TRsub_from_probs <- function(probs_grp,
                                   from_age,
                                   age_interval,
                                   init_for_surv,
                                   trans_col = "transition") {
  
  # Ensure age grid subset
  probs_sub <- probs_grp[probs_grp$age >= from_age, , drop = FALSE]
  
  # Transition label, compatible with mscalc
  probs_sub[[trans_col]] <- paste(probs_sub$from, probs_sub$to, sep = "-")
  
  # State structure: 2 living + 1 dead
  st <- detect_three_state_structure(probs_sub)
  s1 <- st$living[1L]
  s2 <- st$living[2L]
  sd <- st$dead
  
  # Survivorship (lx-scale) from mscalc
  surv <- mscalc::calc_survivorship(
    p_list       = probs_sub,
    init         = init_for_surv,
    delim        = "-",
    age_interval = age_interval,
    trans_col    = trans_col
  )
  # surv must have at least: age, <s1>, <s2>, <sd>
  needed_cols <- c("age", s1, s2, sd)
  if (!all(needed_cols %in% names(surv))) {
    stop("calc_survivorship output must contain columns: ",
         paste(needed_cols, collapse = ", "))
  }
  
  # Merge survivorship into transition rows by age and from-state
  # (we want lx_i at each age, for every i->j row)
  surv_long <- surv[, c("age", s1, s2, sd)]
  surv_long <- tidyr::pivot_longer(
    surv_long,
    cols      = c(s1, s2, sd),
    names_to  = "state",
    values_to = "lx"
  )
  
  probs_with_lx <- dplyr::left_join(
    probs_sub,
    surv_long,
    by = c("age" = "age", "from" = "state")
  )
  
  if (any(is.na(probs_with_lx$lx))) {
    stop("Missing lx after join; check calc_survivorship / state labels.")
  }
  
  # Transfers: expected number moving i->j in each step, per unit initial mass
  probs_with_lx$tr <- probs_with_lx$lx * probs_with_lx$p
  
  # Pivot transfers wide by from/to pairs for readability
  tr_wide <- probs_with_lx[, c("age", "from", "to", "tr")]
  tr_wide <- tidyr::pivot_wider(
    tr_wide,
    id_cols     = "age",
    names_from  = c("from", "to"),
    names_sep   = "_",
    values_from = "tr",
    values_fill = list(tr = 0)
  )
  
  # Attach survivorship (lx) to that same age grid
  TRsub <- dplyr::left_join(
    surv[, c("age", s1, s2, sd)],
    tr_wide,
    by = "age"
  )
  
  # For readability, rename lx columns to H/U/D style *internally only*
  # but preserve original state names in attributes so we can map back.
  colnames(TRsub)[colnames(TRsub) == s1] <- "Hx"
  colnames(TRsub)[colnames(TRsub) == s2] <- "Ux"
  colnames(TRsub)[colnames(TRsub) == sd] <- "Dx"
  
  attr(TRsub, "state_labels") <- list(
    s1 = s1,
    s2 = s2,
    sd = sd
  )
  
  TRsub
}

init_back_prob_from_TRsub <- function(TRsub,
                                      start_age,
                                      age_interval) {
  
  # Extract original state labels (incdem)
  lab <- attr(TRsub, "state_labels")
  if (is.null(lab)) {
    stop("TRsub must carry 'state_labels' attribute.")
  }
  
  # Expect transfers named like "<from>_<to>" with H/U/D mapped from s1/s2/sd
  required_cols <- c("age", "Hx", "Ux", "Dx",
                     "H_U", "H_D", "H_H",
                     "U_H", "U_D", "U_U")
  if (!all(required_cols %in% colnames(TRsub))) {
    stop("TRsub is missing required transfer columns: ",
         paste(setdiff(required_cols, colnames(TRsub)), collapse = ", "))
  }
  
  ID <- TRsub
  
  # cumulative deaths from living mass (generalising your CDx)
  total0 <- ID$Hx[1L] + ID$Ux[1L]
  CDx    <- total0 - (ID$Hx + ID$Ux)
  Dx_inc <- c(0, diff(CDx))
  
  ID$CDx <- CDx
  ID$Dx_inc <- Dx_inc
  
  # Transfers (already there):
  # H_H, H_U, H_D, U_H, U_U, U_D
  
  # Reverse probabilities (age step direction is upwards)
  n_age <- nrow(ID)
  
  ID$r_hh <- ID$H_H / dplyr::lead(ID$Hx)
  ID$r_hu <- ID$U_H / dplyr::lead(ID$Hx)
  
  ID$r_uu <- ID$U_U / dplyr::lead(ID$Ux)
  ID$r_uh <- ID$H_U / dplyr::lead(ID$Ux)
  
  ID$r_dh <- ID$H_D / dplyr::lead(ID$Dx_inc)
  ID$r_du <- ID$U_D / dplyr::lead(ID$Dx_inc)
  
  # Drop the last step (no lead info there)
  r_hh_vec <- ID$r_hh[seq_len(n_age - 1L)]
  r_hu_vec <- ID$r_hu[seq_len(n_age - 1L)]
  r_uu_vec <- ID$r_uu[seq_len(n_age - 1L)]
  r_uh_vec <- ID$r_uh[seq_len(n_age - 1L)]
  r_dh_vec <- ID$r_dh[seq_len(n_age - 1L)]
  r_du_vec <- ID$r_du[seq_len(n_age - 1L)]
  
  # Construct reverse projection blocks (H/U/D internal labels)
  rHH <- rpi2u(r_hh_vec, "H", "H", start_age, age_interval)
  rUH <- rpi2u(r_uh_vec, "U", "H", start_age, age_interval)
  rDH <- rpi2u(r_dh_vec, "D", "H", start_age, age_interval)
  
  rUU <- rpi2u(r_uu_vec, "U", "U", start_age, age_interval)
  rHU <- rpi2u(r_hu_vec, "H", "U", start_age, age_interval)
  rDU <- rpi2u(r_du_vec, "D", "U", start_age, age_interval)
  
  zeros_block <- matrix(0, nrow = n_age, ncol = 3L * n_age)
  
  rU2 <- rbind(
    cbind(rHH, rUH, rDH),
    cbind(rHU, rUU, rDU),
    zeros_block
  )
  
  rU2[is.na(rU2)] <- 0
  
  # Back-project prevalence of H starting from H, U, D
  rpiH_H <- proj_prev_back(rU2, origin_state_index = 1L, n_age = n_age)
  rpiH_U <- proj_prev_back(rU2, origin_state_index = 2L, n_age = n_age)
  rpiH_D <- proj_prev_back(rU2, origin_state_index = 3L, n_age = n_age)
  
  last_col <- ncol(rpiH_H)
  idx      <- last_col - 1L
  
  pH <- (rpiH_H[1L, idx] +
           rpiH_U[1L, idx] +
           rpiH_D[1L, idx]) / 3
  
  c(H = pH, U = 1 - pH)
}


init_rev_stationary_from_probs <- function(probs_grp,
                                           from_age,
                                           age_interval,
                                           init_for_surv = c(`1` = 1, `2` = 0),
                                           trans_col     = "transition") {
  
  # Build TRsub using mscalc survivorship and incdem transitions
  TRsub <- build_TRsub_from_probs(
    probs_grp    = probs_grp,
    from_age     = from_age,
    age_interval = age_interval,
    init_for_surv = init_for_surv,
    trans_col    = trans_col
  )
  
  # Compute reverse-stationary prevalence in H/U internal space
  init_HU <- init_back_prob_from_TRsub(
    TRsub      = TRsub,
    start_age  = from_age,
    age_interval = age_interval
  )
  
  # Map back to original state labels
  lab <- attr(TRsub, "state_labels")
  s1  <- lab$s1
  s2  <- lab$s2
  
  out <- c(0, 0, 0)
  names(out) <- c(s1, s2, lab$sd)
  
  out[s1] <- init_HU["H"]
  out[s2] <- init_HU["U"]
  
  # Dead starts at 0 by construction
  out
}



align_init_vec <- function(init, probs_grp) {
  states <- sort(unique(c(probs_grp$from, probs_grp$to)))
  
  if (is.null(names(init))) {
    stop("`init` must be a named vector with state labels.")
  }
  
  # sanity: init must refer only to states that exist
  unknown <- setdiff(names(init), states)
  if (length(unknown)) {
    stop("`init` contains states not present in probs: ",
         paste(unknown, collapse = ", "))
  }
  
  # just renormalise on the states given in init
  if (sum(init) <= 0) {
    stop("Sum of `init` must be positive.")
  }
  
  init / sum(init)
}

get_init_vec_from_probs <- function(probs_grp,
                                    from_age,
                                    age_interval,
                                    init,
                                    init_method = c("init",
                                                    "constant_stationary",
                                                    "rev_stationary"),
                                    trans_col = "transition") {
  
  init_method <- match.arg(init_method)
  
  if (init_method == "init") {
    # just check & renormalise, don't add absorbing
    return(align_init_vec(init, probs_grp))
  }
  
  if (init_method == "constant_stationary") {
    return(init_constant_stationary(
      probs_grp = probs_grp,
      from_age  = from_age
    ))
  }
  
  if (init_method == "rev_stationary") {
    return(init_rev_stationary_from_probs(
      probs_grp     = probs_grp,
      from_age      = from_age,
      age_interval  = age_interval,
      init_for_surv = init,
      trans_col     = trans_col
    ))
  }
  
  stop("Unknown init_method: ", init_method)
}


calc_exs <- function(probs,
                     from_age      = 50,
                     age_interval  = 0.25,
                     init          = c(`1` = 1, `2` = 0),
                     init_method   = c("init",
                                       "constant_stationary",
                                       "rev_stationary"),
                     from_col      = "from",
                     to_col        = "to",
                     age_col       = "age",
                     p_col         = "p",
                     group_cols    = c("replicate", "female", "period5"),
                     trans_col     = "transition") {
  
  init_method <- match.arg(init_method)
  
  # standardise cols
  from_sym   <- rlang::sym(from_col)
  to_sym     <- rlang::sym(to_col)
  age_sym    <- rlang::sym(age_col)
  p_sym      <- rlang::sym(p_col)
  group_syms <- rlang::syms(group_cols)
  trans_sym  <- rlang::sym(trans_col)
  
  probs_std <- probs |>
    dplyr::rename(
      from = !!from_sym,
      to   = !!to_sym,
      age  = !!age_sym,
      p    = !!p_sym
    )
  
  probs_std |>
    dplyr::filter(age >= from_age) |>
    dplyr::group_by(!!!group_syms) |>
    dplyr::group_modify(
      ~ {
        # .x contains: from, to, age, p
        # .y contains: replicate, female, period5
        
        grp_full <- dplyr::bind_cols(.y, .x)
        
        # 1. initial distribution (full transitions; includes self + absorbing)
        init_vec <- get_init_vec_from_probs(
          probs_grp    = grp_full,
          from_age     = from_age,
          age_interval = age_interval,
          init         = init,
          init_method  = init_method,
          trans_col    = trans_col
        )
        
        # 2. off-diagonal transitions for calc_occupancies()
        p_list <- grp_full |>
          dplyr::filter(to > from) |>
          dplyr::mutate(!!trans_sym := paste0(from, "-", to)) |>
          dplyr::select(!!!group_syms, age, p, !!trans_sym)
        
        # 3. mscalc::calc_occupancies() (already Lx-like, discounted)
        mscalc::calc_occupancies(
          transitions  = p_list,
          init         = init_vec,
          delim        = "-",
          age_interval = age_interval,
          trans_col    = trans_col
        )
      }
    ) |>
    dplyr::rename(
      Lx1 = `1`,
      Lx2 = `2`
    ) |>
    dplyr::group_by(!!!group_syms) |>
    dplyr::summarise(
      DFLE = sum(Lx1),
      DLE  = sum(Lx2),
      .groups = "drop"
    )
}
















