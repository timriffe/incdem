
source("R/00_package_and_functions.R")



hazards_to_discrete <- function(hazards,
                                age_interval = 0.25,
                                id_cols = NULL,
                                state_space = NULL,
                                haz_col = "haz",
                                from_col = "from",
                                to_col = "to",
                                n_cores = 1,
                                parallel = c("auto", "none", "future", "mclapply")) {
  parallel <- match.arg(parallel)
  
  df <- tibble::as_tibble(hazards)
  
  # Which columns define a unique Q-matrix?
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
    Q <- matrix(0, n_states, n_states)
    
    f <- idx[as.character(grp_df[[from_col]])]
    t <- idx[as.character(grp_df[[to_col]])]
    
    Q[cbind(f, t)] <- grp_df[[haz_col]]
    
    diag(Q) <- -rowSums(Q)
    
    P <- expm::expm(Q * age_interval)
    
    tibble::tibble(
      from = rep(state_space, each = n_states),
      to   = rep(state_space, times = n_states),
      p    = as.vector(P)
    )
  }
  
  # Decide parallel strategy
  if (parallel == "auto") {
    parallel <- if (n_cores <= 1L) "none" else "future"
  }
  
  if (parallel == "future") {
    current_strategy <- class(future::plan())[1]
    if (grepl("sequential", tolower(current_strategy))) {
      stop(
        paste0(
          "parallel = 'future' was requested, but the current future plan is sequential.\n",
          "Please run this in your console BEFORE calling hazards_to_discrete():\n\n",
          "  future::plan(future::multisession, workers = ", n_cores, ")\n\n",
          "Or choose a different parallel method: parallel = 'none' or 'mclapply'.\n"
        ),
        call. = FALSE
      )
    }
  }
  
  # Grouping structure without group_split()
  if (length(id_cols) == 0L) {
    df_grouped <- dplyr::group_by(df, ..dummy = 1L)
    id_cols    <- "..dummy"
  } else {
    df_grouped <- dplyr::group_by(df, dplyr::across(dplyr::all_of(id_cols)))
  }
  
  gindex     <- dplyr::group_indices(df_grouped)  # integer vector, same length as df
  group_keys <- dplyr::group_keys(df_grouped)     # one row per group
  n_groups   <- nrow(group_keys)
  
  compute_one_group <- function(i) {
    grp_df  <- df[gindex == i, , drop = FALSE]
    key_row <- group_keys[i, , drop = FALSE]
    
    base <- build_P_one_group(grp_df)
    
    key_rep <- key_row[rep(1, nrow(base)), , drop = FALSE]
    dplyr::bind_cols(key_rep, base)
  }
  
  eval_serial <- function() {
    purrr::map_dfr(seq_len(n_groups), compute_one_group)
  }
  
  eval_future <- function() {
    furrr::future_map_dfr(
      seq_len(n_groups),
      compute_one_group,
      .options = furrr::furrr_options(seed = FALSE)
    )
  }
  
  # ---- CHUNKED mclapply: one chunk per core ----
  eval_mclapply <- function() {
    cores <- max(1L, n_cores)
    cores <- min(cores, n_groups)  # don't create more cores than groups
    
    # split group indices into ~equal chunks
    all_idx <- seq_len(n_groups)
    # This divides into 'cores' chunks with roughly equal length
    chunk_id <- ceiling(seq_along(all_idx) * cores / n_groups)
    idx_chunks <- split(all_idx, chunk_id)
    
    chunk_worker <- function(idx_vec) {
      # idx_vec is a vector of group indices for this core
      res_list <- vector("list", length(idx_vec))
      for (k in seq_along(idx_vec)) {
        i <- idx_vec[k]
        res_list[[k]] <- compute_one_group(i)
      }
      dplyr::bind_rows(res_list)
    }
    
    out_list <- parallel::mclapply(
      X        = idx_chunks,
      FUN      = chunk_worker,
      mc.cores = cores
    )
    
    dplyr::bind_rows(out_list)
  }
  
  out <- switch(
    parallel,
    "none"     = eval_serial(),
    "future"   = eval_future(),
    "mclapply" = eval_mclapply(),
    eval_serial()
  )
  
  if ("..dummy" %in% names(out)) {
    out <- dplyr::select(out, -"..dummy")
  }
  
  dplyr::ungroup(out)
}

haz <- read_csv("Data/model1/adj_haz_replicates.csv.gz") |> 
  mutate(from = substr(transition,2,2) |> as.integer(),
         to = substr(transition,3,3) |> as.integer()) |> 
  select(-transition)
object.size(haz) |> print(units="Mb")

probs1f <- 
  haz |> 
  filter(female == 1) |> 
  hazards_to_discrete(
    age_interval = 0.25,
    id_cols      = c("age", "replicate", "female", "period5"),
    n_cores      = 7,
    parallel     = "mclapply")

write_csv(probs1f, file = "Data/model1/probs_f.csv.gz")
rm(probs1f);gc()

probs1m <- 
  haz |> 
  filter(female == 0) |> 
  hazards_to_discrete(
    age_interval = 0.25,
    id_cols      = c("age", "replicate", "female", "period5"),
    n_cores      = 6,
    parallel     = "mclapply")

write_csv(probs1m, file = "Data/model1/probs_m.csv.gz")
rm(probs1m);gc()
