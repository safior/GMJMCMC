#' @export 
#' Summary function that works with gmjmcmc.re.
#' This function summaries the correlation features in addition to the regular features.
summary.gmjmcmc.re <- function (object, pop = "best", tol = 0.0001, labels = FALSE, effects = NULL, data = NULL, ...) {
  transforms.bak <- set.transforms(object$transforms)
  if (pop == "all") {
    results <- list()
    results[[1]] <- object
    merged <- merge_results.ts(results, pop, 2, 0.0000001, data = data)
    
    best <- max(sapply(merged$results, function (y) y$best))
    feats.strings <- sapply(merged$features, FUN = function(x) print.feature.ts(x = x, labels = labels, round = 2))
    
    if (!is.null(effects) & !is.null(labels)) {
      effects <- compute_effects(merged,labels = labels, quantiles = effects)
    }
    
    return(summary_internal(best = merged$crit.best, feats.strings, merged$marg.probs, effects = effects,
                     best.pop = merged$pop.best, thread.best = merged$thread.best,  
                     reported = merged$reported, rep.pop = merged$rep.pop, rep.thread = merged$rep.thread, tol = tol))
  }
  
  if (pop == "last") pop <- length(object$models)
  else if (pop == "best") {
    pop <- which.max(unlist(object$best.margs))
    feats.strings <- sapply(object$populations[[pop]], FUN = function(x) print.feature.ts(x = x, labels = labels, round = 2))
    re.feats.strings <- sapply(object$re_populations[[pop]], FUN = function(x) print.feature.re(x))
    feats.strings <- c(feats.strings, re.feats.strings)
  }
  
  if (!is.null(effects) & !is.null(labels)) {
    effects <- compute_effects(object, labels = labels, quantiles = effects)
  }
  
  obj <- summary_internal(
    best = object$best,
    marg.probs = object$marg.probs[[pop]],
    effects = effects,
    feats.strings = feats.strings,
    best.pop = which.max(unlist(object$best.margs)),
    reported = object$best.margs[[pop]],
    rep.pop = pop,
    tol = tol
  )
  set.transforms(transforms.bak)
  return(obj)
}

# Summary function for gmjmcmc.parallel.re
summary.gmjmcmc_merged.re <- function (object, tol = 0.0001, labels = FALSE, effects = NULL, pop = NULL, data_ts, window_list, add_lagged_response, ...) {
  transforms.bak <- set.transforms(object$transforms)
  transforms.bak.ts <- set.transforms(object$transforms)

  if (!is.null(pop)) {
    lookback_window <- object[[1]]$lookback_window
    object <- merge_results.re2(object$results.raw, populations = pop, complex.measure = 2, tol = 0.0000001, data = data, lw = lookback_window)
  }
  
  best <- max(sapply(object$results, function (y) y$best))
  feats.strings <- sapply(object$features, FUN = function(x) print.feature.ts(x = x, labels = labels, round = 2))
  re.feats.strings <- sapply(object$re.features, FUN = function(x) print.feature.re(x))
  feats.strings <- c(feats.strings, re.feats.strings)
  
  if (!is.null(effects) & !is.null(labels)) {
    effects <- compute_effects(object,labels = labels, quantiles = effects)
  }
  obj <- summary_internal(best = object$crit.best, feats.strings, object$marg.probs, effects = effects,
                   best.pop = object$pop.best, thread.best = object$thread.best,  
                   reported = object$reported, rep.pop = object$rep.pop, rep.thread = object$rep.thread, tol = tol)
  set.transforms(transforms.bak)
  set.transforms.ts(transforms.bak.ts)
  return(obj)
}

#' @export merge_results.re2
#' Naming is due to this being version 2
merge_results.re2 <- function (results, populations = NULL, complex.measure = NULL, tol = NULL, data = NULL, 
                                lw, marg.lik.method) {
  # Default values
  if (is.null(populations))
    populations <-"best"
  if (is.null(complex.measure))
    complex.measure <- 2
  if (is.null(tol))
    tol <- 0.0000001

  add_lagged_response <- NULL
  i <- 1
  # Checks if lagged repsonse is to be added as covariate. While loop in case first run failed.
  while(is.null(add_lagged_response)) {
    if(!is.atomic(results[[i]])) {
      add_lagged_response <- results[[i]]$add_lagged_response
    }
    i <- i + 1
  }

  # Check and filter results that did not run successfully
  results <- filter.results(results)
  raw.results <- results
  res.count <- length(results)

  # Get correlation features for each population and corresponding indices in marginal probability vector.
  re.out <- get_pop_re.features(results, populations)
  re.features <- re.out$re.features
  re.inds.tot <- re.out$re.inds.tot

  # Get renomarlized features
  renorm_feat <- get_renomarlized_features(results, populations)
  features <- renorm_feat$features
  renorms <- renorm_feat$renorms
  results <- renorm_feat$results
  pw <- renorm_feat$pw

  renorms <- unlist(renorms)
  na.feats <- which(is.na(renorms))
  if (length(na.feats) != 0) {
    warning("Underflow occurred,", length(na.feats), "features removed.\n")
    renorms <- renorms[-na.feats]
    features <- features[-na.feats]
  }

  ## Detect equivalent features
  # Mock data not implemented currently
  # Account for possible lagged response covariate
  if (add_lagged_response) {
    data <- add_lag_resp(data) 
  }
  data2 <- data[(lw + 1) : nrow(data), ]
  mock.data <- check.data(data2, FALSE)
  mock.data.ts <- check.data(data, FALSE)

  mock.data.precalc <- precalc.features.ts(mock.data, mock.data.ts, lw, features)[,-(1:2)]

  # Renormalize regular features
  feats.map <- get_feats.map(mock.data.precalc, features, renorms[-re.inds.tot], tol)

  # Select the simplest features based on the specified complexity measure and sort them
  feats.simplest.ids <- unique(feats.map[complex.measure, ])
  feats.simplest.ids <- feats.simplest.ids[order(feats.map[4, feats.simplest.ids])]
  counts <- sapply(feats.simplest.ids, function(x) sum(feats.map[complex.measure,] == x))
  feats.simplest <- features[feats.simplest.ids]
  importance <- feats.map[4, feats.simplest.ids, drop = FALSE]

  # Sum estimated marginal inclusion probabilities of identical correlation features after renormalization
  renorms.re <- renorms[re.inds.tot]

  re.sums <- sum_same_re(re.features, renorms.re, marg.lik.method)
  # Assign unique correlation features
  re.features.new <- re.sums$re.features.new
  # Convert vector to matrix and concatenate with matrix of recalculated marginal inclusion probabilites for the regular features
  re.importance <- matrix(re.sums$re.importance, nrow = 1, ncol = length(re.sums$re.importance))
  importance <- cbind(importance, re.importance)

  # Get best results
  best <- get_best_results(results)

  merged <- list(
    re.features = re.features.new,
    features = feats.simplest,
    marg.probs = importance,
    counts = counts,
    results = results,
    results.raw = raw.results,
    pop.best = best$pop.best,
    thread.best = best$thread.best,
    crit.best = best$crit.best,
    reported = pw$best,
    rep.pop = pw$pop.best,
    best.log.posteriors = best$bests,
    rep.thread = pw$thread.best,
    transforms = results[[1]]$transforms,
    transforms.ts = results[[1]]$transforms.ts,
    add_lagged_response = results[[1]]$add_lagged_respons,
    lookback_window = lw
  )
  attr(merged, "class") <- "gmjmcmc_merged"
  return(merged)
}

# Function for summing identical correlation features
sum_same_re <- function(re.features, renorms.re, marg_lik_method) {
  string.re.features <- get.string.re.feat(re.features, marg_lik_method)
  summed.inds <- c()
  re.features.new <- list()
  re.importance <- c()
  for (i in 1:length(re.features)) {
    s.re.feat <- string.re.features[i]
    # Match current correlation feature with all identical features. 
    same.re.inds <- which(string.re.features == s.re.feat)
    # Account for already summed indices, so now double summing occurs.
    if(i %in% summed.inds) next
    # Append new correlation feature
    re.features.new <- c(re.features.new, re.features[i])
    # Sum over the same correlation features
    same.ind.sum <- sum(unlist(renorms.re[same.re.inds]))
    # Append the sum
    re.importance <- c(re.importance, same.ind.sum)
    # Append summed indices
    summed.inds <- c(summed.inds, same.re.inds)
  }
  return(list(re.features.new = re.features.new, re.importance = re.importance))
}

# Get string representation of correlation features
get.string.re.feat <- function(re.features, marg_lik_method) {
  if (marg_lik_method == 'nlme') {
    string.re.features <- vector("list")
    for(i in 1:length(re.features)){
      re.feat <- re.features[[i]]
      # Concatenate string of cor struct argument and grouping argument
      string.re <- paste(re.feat$cor_arg, re.feat$random_arg)
      string.re.features <- append(string.re.features, string.re)
    }
  }
  # inla is already a uniue string.
  else if(marg_lik_method == 'inla') {
    string.re.features <- re.features
  }
  else {
    stop("Marginal likelihood method is not implemented!")
  }
  return(unlist(string.re.features))
}

# Get a list of all correlation features and their corresponding incices in the marg.probs vector.
get_pop_re.features <- function(results, populations) {
  res.count <- length(results)
  # Select populations to use
  pops.use <- select_pops(results, res.count, populations)
  
  # Collect all features and their renormalized weighted values
  re.features <- vector("list")
  re.inds.tot <- c()
  tot.feat <- 0
  for (i in 1:res.count) {
    for (pop in pops.use[[i]]) {
      # Get all correlation features for all populations
      re.feats <- results[[i]]$re_populations[[pop]]
      re.features <- append(re.features, re.feats)

      # Get indices for the correlation features
      # Find the correct index to delimit the populations
      tot.feat.pop <- length(results[[i]]$marg.probs[[pop]])
      tot.feat <- tot.feat + tot.feat.pop
      # Correlation feature indices always corresponds to the last elements in the marg.probs vector. 
      # Recall that the marg.probs vector is a binary vector over the correlation features as well.
      re.inds <- (tot.feat - length(re.feats) + 1) : tot.feat
      re.inds.tot <- c(re.inds.tot, re.inds)
    }
  }
  return(list(re.features = re.features, re.inds.tot = re.inds.tot))
}