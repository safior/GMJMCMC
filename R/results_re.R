#' @export 
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
merge_results.re2 <- function (results, populations = NULL, complex.measure = NULL, tol = NULL, data = NULL, lw, marg.lik.method) {
  # Default values
  if (is.null(populations))
    populations <-"best"
  if (is.null(complex.measure))
    complex.measure <- 2
  if (is.null(tol))
    tol <- 0.0000001

  add_lagged_response <- results[[1]]$add_lagged_response

  # Check and filter results that did not run successfully
  results <- filter.results(results)
  raw.results <- results
  res.count <- length(results)

  # Get random effects and indices in marginal probability vector corresponding to random effects
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
  # Generate mock data to compare features with
  if (add_lagged_response) {
    data <- add_lag_resp(data) 
  }
  data2 <- data[(lw + 1) : nrow(data), ]
  #print(data)
  mock.data <- check.data(data2, FALSE)
  mock.data.ts <- check.data(data, FALSE)

  mock.data.precalc <- precalc.features.ts(mock.data, mock.data.ts, lw, features)[,-(1:2)]

  # Renorms
  feats.map <- get_feats.map(mock.data.precalc, features, renorms[-re.inds.tot], tol)#renorms[1:(length(renorms)-n.re)], tol)

  # Select the simplest features based on the specified complexity measure and sort them
  feats.simplest.ids <- unique(feats.map[complex.measure, ])
  feats.simplest.ids <- feats.simplest.ids[order(feats.map[4, feats.simplest.ids])]
  counts <- sapply(feats.simplest.ids, function(x) sum(feats.map[complex.measure,] == x))
  feats.simplest <- features[feats.simplest.ids]
  importance <- feats.map[4, feats.simplest.ids, drop = FALSE]

  # Sum  marginal likelihood of identical random effects
  renorms.re <- renorms[re.inds.tot]
  #print(results[[1]]$marg.lik.method)
  #marg.lik.method <- results[[1]]$marg.lik.method
  re.sums <- sum_same_re(re.features, renorms.re, marg.lik.method)#results[[1]]$marg_lik_method)
  # Assign unique random effects
  re.features.new <- re.sums$re.features.new
  # Convert vector to matrix and concatenate with matrix of recalculated marginal likelihoods for the regular features
  re.importance <- matrix(re.sums$re.importance, nrow = 1, ncol = length(re.sums$re.importance))
  importance <- cbind(importance, re.importance)
  # print(length(feats.simplest))
  # print(length(renorms))
  # print(length(re.features.new))
  # print(length(re.features))
  # print(length(importance))

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
    transforms.ts = results[[1]]$transforms.ts
  )
  attr(merged, "class") <- "gmjmcmc_merged"
  return(merged)
}

sum_same_re <- function(re.features, renorms.re, marg_lik_method) {
  string.re.features <- get.string.re.feat(re.features, marg_lik_method)
  summed.inds <- c()
  re.features.new <- list()
  re.importance <- c()
  for (i in 1:length(re.features)) {
    s.re.feat <- string.re.features[i]
    #print(s.re.feat)
    same.re.inds <- which(string.re.features == s.re.feat)
    if(i %in% summed.inds) next
    re.features.new <- c(re.features.new, re.features[i])
    #print(sum(unlist(renorms.re[same.re.inds])))
    same.ind.sum <- sum(unlist(renorms.re[same.re.inds]))
    re.importance <- c(re.importance, same.ind.sum)
    summed.inds <- c(summed.inds, same.re.inds)
  }
  return(list(re.features.new = re.features.new, re.importance = re.importance))
}

get.string.re.feat <- function(re.features, marg_lik_method) {
  if (marg_lik_method == 'nlme') {
    string.re.features <- vector("list")
    for(i in 1:length(re.features)){
      re.feat <- re.features[[i]]
      string.re <- paste(re.feat$cor_arg, re.feat$random_arg)
      string.re.features <- append(string.re.features, string.re)
    }
  }
  else if(marg_lik_method == 'inla') {
    string.re.features <- re.features
  }
  else {
    stop("Marginal likelihood method is not implemented!")
  }
  return(unlist(string.re.features))
}

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
      # Get all random effects for all populations
      re.feats <- results[[i]]$re_populations[[pop]]
      re.features <- append(re.features, re.feats)

      # Get indices for the random effects
      tot.feat.pop <- length(results[[i]]$marg.probs[[pop]])
      tot.feat <- tot.feat + tot.feat.pop
      # Random effects indices always corresponds to the last elements in the marg.probs vector
      re.inds <- (tot.feat - length(re.feats) + 1) : tot.feat
      re.inds.tot <- c(re.inds.tot, re.inds)
    }
  }
  return(list(re.features = re.features, re.inds.tot = re.inds.tot))
}

#' @export merge_results.re
# merge_results.re <- function (results, populations = NULL, complex.measure = NULL, tol = NULL, data = NULL, lw = NULL) {
#   # Default values
#   if (is.null(populations))
#     populations <-"best"
#   if (is.null(complex.measure))
#     complex.measure <- 2
#   if (is.null(tol))
#     tol <- 0.0000001

#   # Check and filter results that did not run successfully
#   results <- filter.results(results)
#   raw.results <- results
#   res.count <- length(results)

#   # Select populations to use
#   res.lengths <- vector("list")
#   for (i in 1:res.count) {
#     res.lengths[[i]] <- length(results[[i]]$populations)
#   }
#   if (populations == "last") pops.use <- res.lengths
#   else if (populations == "all") pops.use <- lapply(res.lengths, function(x) 1:x)
#   else if (populations == "best") pops.use <- lapply(1:res.count, function(x) which.max(unlist(results[[x]]$best.marg)))

#   # Get the population weigths to be able to weight the features
#   pw <- population.weigths(results, pops.use)
#   pop.weights <- pw$weights
  
#   bests <- matrix(data = 0, ncol = length(results), nrow = length(results[[1]]$populations))
#   crit.best <- -Inf
#   pop.best <- 1
#   thread.best <- 1
#   for (i in seq_along(results)) {
#     for (pop in 1:(length(results[[i]]$populations))) {
#       bests[pop, i] <- results[[i]]$best.margs[[pop]]
#       if (results[[i]]$best.margs[[pop]] > crit.best) {
#         crit.best <- results[[i]]$best.margs[[pop]]
#         pop.best <- pop
#         thread.best <- i
#       }
#     }
#   }
  
#   # Collect all features and their renormalized weighted values
#   features <- vector("list")
#   re.features <- vector("list")
#   renorms <- vector("list")
#   weight_idx <- 1
#   for (i in 1:res.count) {
#     results[[i]]$pop.weights <- rep(NA, length(results[[i]]$populations))
#     results[[i]]$model.probs <- list()
#     for (pop in pops.use[[i]]) {
#       features <- append(features, results[[i]]$populations[[pop]])
#       re.features <- append(re.features, results[[i]]$re_populations[[pop]])
#       renorms <- append(renorms, pop.weights[weight_idx] * results[[i]]$marg.probs[[pop]])
#       results[[i]]$pop.weights[pop] <- pop.weights[weight_idx]
#       weight_idx <- weight_idx + 1

#       model.probs <- marginal.probs.renorm(results[[i]]$models[[pop]], "models")
#       results[[i]]$model.probs[[pop]] <- model.probs$probs
#       results[[i]]$models[[pop]] <- results[[i]]$models[[pop]][model.probs$idx]
#     }
#     accept.tot <- results[[i]]$accept.tot
#     best <- results[[i]]$best
#     for (item in names(results[[i]])) {
#       if (!(item %in% (c("accept.tot", "best", "transforms")))) results[[i]][[item]] <- results[[i]][[item]][pops.use[[i]]]
#     }
#     results[[i]]$accept.tot <- accept.tot
#     results[[i]]$best <- best
#   }

#   renorms <- unlist(renorms)
#   na.feats <- which(is.na(renorms))
#   if (length(na.feats) != 0) {
#     warning("Underflow occurred,", length(na.feats), "features removed.\n")
#     renorms <- renorms[-na.feats]
#     features <- features[-na.feats]
#     re.features <- re.features[-na.feats]
#   }

#   renorms <- matrix(renorms, nrow = 1, ncol = length(renorms))

#   merged <- list(
#     features = features,
#     re.features = re.features,
#     marg.probs = renorms,
#     results = results,
#     results.raw = raw.results,
#     pop.best = pop.best,
#     thread.best = thread.best,
#     crit.best = crit.best,
#     reported = pw$best,
#     rep.pop = pw$pop.best,
#     best.log.posteriors = bests,
#     rep.thread = pw$thread.best,
#     transforms = results[[1]]$transforms
#   )
#   attr(merged, "class") <- "gmjmcmc_merged"
#   return(merged)
# }