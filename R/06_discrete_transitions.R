

library(vroom)
library(dplyr)
library(purrr)
library(expm)

# untested
hazards_to_discrete <- function(hazards,
                          age_interval = .25,
                          id_cols = NULL,
                          state_space = NULL,
                          haz_col = "haz",
                          from_col = "from",
                          to_col = "to") {
  # hazards: data.frame / tibble in long format
  
  df <- as_tibble(hazards)
  
  # Which columns define a unique Q-matrix?
  # Default: all non-(from,to,haz) columns (including age, boot id, etc.)
  if (is.null(id_cols)) {
    id_cols <- setdiff(
      names(df),
      c(from_col, to_col, haz_col)
    )
  }
  
  # State labels
  if (is.null(state_space)) {
    state_space <- sort(unique(c(df[[from_col]], df[[to_col]])))
  }
  n_states <- length(state_space)
  
  # Map original state labels to 1..n
  idx <- seq_len(n_states)
  names(idx) <- as.character(state_space)
  
  build_P_one_group <- function(grp_df) {
    # Build Q for this group
    Q <- matrix(0, n_states, n_states)
    
    f <- idx[as.character(grp_df[[from_col]])]
    t <- idx[as.character(grp_df[[to_col]])]
    
    # Fill off-diagonals with hazards
    Q[cbind(f, t)] <- grp_df[[haz_col]]
    
    # Diagonals = -row sum of off-diagonals
    diag(Q) <- -rowSums(Q)
    
    # Matrix exponential for interval length age_interval
    P <- expm(Q * age_interval)
    
    tibble(
      from = rep(state_space, each = n_states),
      to   = rep(state_space, times = n_states),
      p    = as.vector(P)
    )
  }
  
  df %>%
    group_by(across(all_of(id_cols))) %>%
    group_modify(~ build_P_one_group(.x)) %>%
    ungroup()
}

