#' @export 
summary.gmjmcmc.ts <- function (object, pop = "best", tol = 0.0001, labels = FALSE, effects = NULL, data = NULL, ...) {
  transforms.bak <- set.transforms(object$transforms)
  if (pop == "all") {
    results <- list()
    results[[1]] <- object
    merged <- merge_results(results, pop, 2, 0.0000001, data = data)
    
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
  else if (pop == "best") pop <- which.max(unlist(object$best.margs))
  feats.strings <- sapply(object$populations[[pop]], FUN = function(x) print.feature.ts(x = x, labels = labels, round = 2))
  
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

#' @export
string.population.ts <- function(x, round = 2) {
  cbind(sapply(x, print.feature.ts, round = round))
}

#' @export merge_results
merge_results.ts <- function (results, populations = NULL, complex.measure = NULL, tol = NULL, data = NULL, lw) {
  # Default values
  if (is.null(populations))
    populations <-"best"
  if (is.null(complex.measure))
    complex.measure <- 2
  if (is.null(tol))
    tol <- 0.0000001

  # Check and filter results that did not run successfully
  results <- filter.results(results)
  raw.results <- results
  res.count <- length(results)

  # Select populations to use
  res.lengths <- vector("list")
  for (i in 1:res.count) res.lengths[[i]] <- length(results[[i]]$populations)
  if (populations == "last") pops.use <- res.lengths
  else if (populations == "all") pops.use <- lapply(res.lengths, function(x) 1:x)
  else if (populations == "best") pops.use <- lapply(1:res.count, function(x) which.max(unlist(results[[x]]$best.marg)))

  # Get the population weigths to be able to weight the features
  pw <- population.weigths(results, pops.use)
  pop.weights <- pw$weights
  
  bests <- matrix(data = 0, ncol = length(results), nrow = length(results[[1]]$populations))
  crit.best <- -Inf
  pop.best <- 1
  thread.best <- 1
  for (i in seq_along(results)) {
    for (pop in 1:(length(results[[i]]$populations))) {
      bests[pop, i] <- results[[i]]$best.margs[[pop]]
      if (results[[i]]$best.margs[[pop]] > crit.best) {
        crit.best <- results[[i]]$best.margs[[pop]]
        pop.best <- pop
        thread.best <- i
      }
    }
  }
  
  # Collect all features and their renormalized weighted values
  features <- vector("list")
  renorms <- vector("list")
  weight_idx <- 1
  for (i in 1:res.count) {
    results[[i]]$pop.weights <- rep(NA, length(results[[i]]$populations))
    results[[i]]$model.probs <- list()
    for (pop in pops.use[[i]]) {
      features <- append(features, results[[i]]$populations[[pop]])
      renorms <- append(renorms, pop.weights[weight_idx] * results[[i]]$marg.probs[[pop]])
      results[[i]]$pop.weights[pop] <- pop.weights[weight_idx]
      weight_idx <- weight_idx + 1

      model.probs <- marginal.probs.renorm(results[[i]]$models[[pop]], "models")
      results[[i]]$model.probs[[pop]] <- model.probs$probs
      results[[i]]$models[[pop]] <- results[[i]]$models[[pop]][model.probs$idx]
    }
    accept.tot <- results[[i]]$accept.tot
    best <- results[[i]]$best
    for (item in names(results[[i]])) {
      if (!(item %in% (c("accept.tot", "best", "transforms")))) results[[i]][[item]] <- results[[i]][[item]][pops.use[[i]]]
    }
    results[[i]]$accept.tot <- accept.tot
    results[[i]]$best <- best
  }
  renorms <- unlist(renorms)
  na.feats <- which(is.na(renorms))
  if (length(na.feats) != 0) {
    warning("Underflow occurred,", length(na.feats), "features removed.\n")
    renorms <- renorms[-na.feats]
    features <- features[-na.feats]
  }
  feat.count <- length(features)

  # Get complexity for all features
  complex <- complex.features(features)

  ## Detect equivalent features
  # Generate mock data to compare features with
  if (is.null(data)) mock.data <- matrix(runif((feat.count + 2)^2, -100, 100), ncol = feat.count + 2)
  else {
    data2 <- data[(lw + 1) : nrow(data), ]
    mock.data <- check.data(data2, FALSE)
    mock.data.ts <- check.data(data, FALSE)
  }
  
  mock.data.precalc <- precalc.features.ts(mock.data, mock.data.ts, lw, features)[,-(1:2)]

  # Calculate the correlation to find equivalent features
  cors <- cor(mock.data.precalc)

  # A map to link equivalent features together,
  # row 1-3 are the simplest equivalent features based on three different complexity measures
  # row 4 is the total weighted density of those features
  feats.map <- matrix(1:feat.count, 4, feat.count, byrow = TRUE)
  for (i in seq_len(nrow(cors))) {
    equiv.feats <- which(cors[i, ] >= (1 - tol))
    # Compare equivalent features complexity to find most simple
    equiv.complex <- list(width=complex$width[equiv.feats], oc=complex$oc[equiv.feats], depth=complex$depth[equiv.feats])
    equiv.simplest <- lapply(equiv.complex, which.min)
    feats.map[1:3,equiv.feats] <- c(equiv.feats[equiv.simplest$width], equiv.feats[equiv.simplest$oc], equiv.feats[equiv.simplest$depth])
    feats.map[4,equiv.feats] <- sum(renorms[equiv.feats])
  }
  # Select the simplest features based on the specified complexity measure and sort them
  feats.simplest.ids <- unique(feats.map[complex.measure, ])
  feats.simplest.ids <- feats.simplest.ids[order(feats.map[4, feats.simplest.ids])]
  counts <- sapply(feats.simplest.ids, function(x) sum(feats.map[complex.measure,] == x))
  feats.simplest <- features[feats.simplest.ids]
  importance <- feats.map[4, feats.simplest.ids, drop = FALSE]
  merged <- list(
    features = feats.simplest,
    marg.probs = importance,
    counts = counts,
    results = results,
    results.raw = raw.results,
    pop.best = pop.best,
    thread.best = thread.best,
    crit.best = crit.best,
    reported = pw$best,
    rep.pop = pw$pop.best,
    best.log.posteriors = bests,
    rep.thread = pw$thread.best,
    transforms = results[[1]]$transforms
  )
  attr(merged, "class") <- "gmjmcmc_merged"
  return(merged)
}