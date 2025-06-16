group_bootstraps2 <- function(data,
                              group,
                              times = 25,
                              apparent = FALSE,
                              ...,
                              strata = NULL,
                              pool = 0.1,
                              weight = NULL) {
  # rlang::check_dots_empty()
  
  group <- rsample:::validate_group({{ group }}, data)
  
  if (!missing(strata)) {
    strata <- rsample:::check_grouped_strata({{ group }}, {{ strata }}, pool, data)
  }
  
  split_objs <- group_boot_splits2(
    data = data,
    group = group,
    times = times,
    strata = strata,
    pool = pool,
    weight = weight
  )
  
  split_objs$splits <- purrr::map(split_objs$splits, rsample:::rm_out)
  
  if (apparent) {
    split_objs <- dplyr::bind_rows(split_objs, apparent(data))
  }
  
  if (is.null(strata)) {
    strata <- FALSE
  }
  
  boot_att <- list(
    times = times,
    apparent = apparent,
    strata = strata,
    pool = pool,
    group = group
  )
  
  new_rset(
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
    strata_vals <- make_strata(strata_vals, pool = pool)
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
    prop = 1,
    replace = TRUE,
    strata = strata_vals,
    weights = wt
  )
  
  # Convert group IDs back to row indices
  group_vals <- as.character(group_vals)
  group_to_rows <- split(seq_along(group_vals), group_vals)
  
  # indices <- lapply(indices, function(fold_ids) {
  #   unlist(group_to_rows[as.character(fold_ids)], use.names = FALSE)
  # })
 
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
    # Sanity check: valid positive integer vector?
    if (length(idx) == 0 || any(is.na(idx)) || any(idx <= 0)) {
      stop("Invalid analysis index: empty or contains NA / non-positive value.")
    }
    make_splits(
      list(analysis = idx, assessment = integer(0)),
      data = data,
      class = c("group_boot_split", "boot_split")
    )
  })
  
  all_assessable <- purrr::map_int(split_objs, ~nrow(assessment(.x)))
  if (any(all_assessable == 0)) {
    cli::cli_warn(
      c(
        "Some assessment sets contained zero rows.",
        i = "Consider using a non-grouped resampling method."
      ),
      call = rlang::caller_env()
    )
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
                  keys <- balance_prop2(data_ind = data_ind, v = v, strata = strata, ...)
                  list(data_ind = data_ind, keys = keys)
                }
  )
  
  data_ind <- res$data_ind
  keys <- res$keys
  data_ind$..group <- as.character(data_ind$..group)
  keys$..group <- as.character(keys$..group)
  
  # Ensure each fold has its own data frame of group IDs and replications
  split(keys, keys$..folds)
}


balance_prop2 <- function(prop, data_ind, v, replace = TRUE, strata = NULL, weights = NULL, ...) {
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
    prop = prop,
    v = v,
    replace = replace,
    weights = weights
  )
  
  keys
}

      
balance_prop_helper2 <- function(data_ind, prop, v, replace = TRUE, weights = NULL) {
  group_ids <- sort(unique(as.character(data_ind$..group)))
  
  if (!is.null(weights)) {
    weights_vec <- weights[group_ids]
    weights_vec[is.na(weights_vec)] <- 0
    weights_vec <- weights_vec / sum(weights_vec)
  } else {
    weights_vec <- rep(1, length(group_ids))
    weights_vec <- weights_vec / sum(weights_vec)
  }
  
  purrr::map_dfr(seq_len(v), function(x) {
    sampled_groups <- sample(
      group_ids,
      size = ceiling(length(group_ids) * prop),
      replace = replace,
      prob = weights_vec
    )
    
    # tibble(..group = sampled_groups) %>%
    #   count(..group, name = "..replications") %>%
    #   mutate(..folds = x)
    # Sample n groups (with replacement, weighted)
    sampled_groups <- sample(group_ids, size = length(group_ids), replace = TRUE, prob = weights_vec)
    
    # Record how often each group was sampled
    rep_counts <- table(sampled_groups)
    tibble(
      ..group = names(rep_counts),
      ..replications = as.integer(rep_counts),
      ..folds = x
    )
    
  })
}
