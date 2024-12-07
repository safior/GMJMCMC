#' @export gmjmcmc.ts
gmjmcmc.ts <- function (
  data_ts,
  loglik.pi = gaussian.loglik,
  loglik.alpha = gaussian.loglik.alpha,
  transforms,
  ###############
  ts_transforms,
  window_list,
  add_lagged_response = TRUE,
  ###############
  P = 10,
  N.init = 100,
  N.final = 100,
  probs = NULL,
  params = NULL,
  sub = FALSE,
  verbose = TRUE
) {

  # Must discard starting values corresponding to the lookback window
  # Must copy the covariates and response and somehow make them available to the algorithm to compute time series features
  if (add_lagged_response) {
    lag_resp <- c(NA, data_ts[1:(nrow(data_ts)-1), 1])
    data_ts <- cbind(data_ts, lagged_response = lag_resp)
    # Remove start to ensure time series transforms are computable.
    # 2 because of NA introduced above
    data_ts <- data_ts[2:nrow(data_ts), ]
  }
  #print(head(data_ts))
  # Find the largest possible lookback over all transforms
  lookback_window <- max(unlist(window_list))
  #print(lookback_window)
  # +1 because we need a full window of covariates
  data <- data_ts[(lookback_window + 1) : nrow(data_ts), ]
  #print(head(data))
  #data[1:start, 1] <- NA
  #print(data)

  # Verify that the data is well-formed
  data <- check.data(data, verbose)
  data_ts <- check.data(data_ts, verbose)

  # Generate default probabilities and parameters if there are none supplied.
  if (is.null(probs)) probs <- gen.probs.gmjmcmc.ts(ts_transforms, transforms)
  if (is.null(params)) params <- gen.params.gmjmcmc(data)
  # Extract labels from column names in dataframe
  labels <- get.labels(data, verbose)
  # Set the transformations option
  set.transforms(transforms)
  set.transforms.ts(ts_transforms)
  # Acceptance probability per population
  accept <- vector("list", P)
  accept <- lapply(accept, function (x) x <- 0)
  # A list of populations that have been visited
  S <- vector("list", P)
  # A list of models that have been visited, refering to the populations
  models <- vector("list", P)
  lo.models <- vector("list", P)
  # A list of all the marginal probabilities for the features, per population
  marg.probs <- vector("list", P)
  # A list of all the marginal probabilities for the models, per population
  model.probs <- vector("list", P)
  # A list of all the indices of the models which the marginal probabilities for the models refer to, per population
  model.probs.idx <- vector("list", P)
  # A list of all the best marginal model likelihoods, per population
  best.margs <- vector("list", P)

  # Create first population
  # How can you add time series features here?
  F.0 <- gen.covariates(ncol(data) - 2)
  if (is.null(params$prel.select))
    S[[1]] <- F.0
  else
    S[[1]] <- F.0[params$prel.select]

  complex <- complex.features(S[[1]])

  ### Main algorithm loop - Iterate over P different populations
  for (p in seq_len(P)) {
    # Set population iteration count
    if (p != P) N <- N.init
    else N <- N.final
    # Precalculate covariates and put them in data.t
    if (length(params$feat$prel.filter) > 0 | p != 1) data.t <- precalc.features.ts(data, data_ts, lookback_window, S[[p]])
    else {
      data.t <- data
      print("ISSUE???")
    }
    
    # Initialize first model of population
    model.cur <- as.logical(rbinom(n = length(S[[p]]), size = 1, prob = 0.5))
    #print(model.cur)
    model.cur.res <- loglik.pre(loglik.pi, model.cur, complex, data.t, params$loglik)
    model.cur <- list(prob = 0, model = model.cur, coefs = model.cur.res$coefs, crit = model.cur.res$crit, alpha = 0)
    best.crit <- model.cur$crit # Reset first best criteria value

    # Run MJMCMC over the population
    if (verbose) print(paste("Population", p, "begin."))
    mjmcmc_res <- mjmcmc.loop(data.t, complex, loglik.pi, model.cur, N, probs, params, sub, verbose)
    if (verbose) cat(paste("\nPopulation", p, "done.\n"))

    # Add the models visited in the current population to the model list
    models[[p]] <- mjmcmc_res$models
    lo.models[[p]] <- mjmcmc_res$lo.models
    # Store marginal likelihoods for current features
    marg.probs[[p]] <- mjmcmc_res$marg.probs
    # Store marginal likelihoods for the visited models
    model.probs[[p]] <- mjmcmc_res$model.probs
    # Store indices for which the marginal likelihoods for the visited models refer to
    model.probs.idx[[p]] <- mjmcmc_res$model.probs.idx
    # Store best marginal model probability for current population
    best.margs[[p]] <- mjmcmc_res$best.crit
    # Print the marginal posterior distribution of the features after MJMCMC
    if (verbose) {
      cat(paste("\rCurrent best crit:", mjmcmc_res$best.crit, "\n"))
      cat("Feature importance:\n")
      print_dist(marg.probs[[p]], sapply(S[[p]], print.feature.ts, labels = labels, round = 2), probs$filter)
    }
    if (params$rescale.large) prev.large <- params$large
    # Generate a new population of features for the next iteration (if this is not the last)
    if (p != P) {
      S[[p + 1]] <- gmjmcmc.transition.ts(S[[p]], F.0, data, data_ts, ts_transforms, window_list, loglik.alpha, marg.probs[[1]], marg.probs[[p]], labels, probs, params$feat, verbose)
      complex <- complex.features(S[[p + 1]])
      if (params$rescale.large) params$large <- lapply(prev.large, function(x) x * length(S[[p + 1]]) / length(S[[p]]))
    }
  }
  # Calculate acceptance rate
  accept.tot <- sum(unlist(accept)) / (N.init * (P - 1) + N.final)
  accept <- lapply(accept, function (x) x / N.init)
  accept[[P]] <- accept[[P]] * N.init / N.final
  # Return formatted results
  results <- list(
    models = models,                   # All models per population
    lo.models = lo.models,             # All local optim models per population
    populations = S,                   # All features per population
    marg.probs = marg.probs,           # Marginal feature probabilities per population
    model.probs = model.probs,         # Marginal feature probabilities per population
    model.probs.idx = model.probs.idx, # Marginal feature probabilities per population
    best.margs = best.margs,           # Best marginal model probability per population
    accept = accept,                   # Acceptance rate per population
    accept.tot = accept.tot,           # Overall acceptance rate
    best = max(unlist(best.margs)),    # Best marginal model probability throughout the run
    transforms = transforms,           # Transformations used by the model
    transforms.ts = ts_transforms,
    lookback_window = lookback_window,
    add_lagged_response = add_lagged_response
  )
  attr(results, "class") <- "gmjmcmc"
  return(results)
}


gmjmcmc.transition.ts <- function(
  S.t,
  F.0,
  data,
  data.ts, 
  ts_transforms,
  window_list,
  loglik.alpha,
  marg.probs.F.0,
  marg.probs,
  labels,
  probs,
  params,
  verbose = TRUE) {
  # Sample which features to keep based on marginal inclusion below probs$filter
  feats.keep <- as.logical(rbinom(n = length(marg.probs), size = 1, prob = pmin(marg.probs / probs$filter, 1)))
  #print(lookback_window)
 
  # Always keep original covariates if that setting is on
  if (params$keep.org) {
    if (params$prel.filter > 0) {
      # Do preliminary filtering if turned on
      feats.keep[(seq_along(F.0))[marg.probs.F.0 > params$prel.filter]] <- T
    } # Keep all if no preliminary filtering
    else feats.keep[seq_along(F.0)] <- T
  }


  # Avoid removing too many features
  if (length(feats.keep) > 0 && mean(feats.keep) < params$keep.min & sum(feats.keep) < params$pop.max/2) {
    feats.add.n <- round((params$keep.min - mean(feats.keep)) * length(feats.keep))
    feats.add <- sample(which(!feats.keep), feats.add.n)
    if((length(feats.add) + sum(feats.keep))>=params$pop.max)
      feats.keep[feats.add] <- T
  }
  
  if(sum(feats.keep)>params$pop.max)
  {
    warning("Number of features to keep greater than pop.max! 
            Continue with pop.max features!
            \n Check your tuning parameters!")
    feats.keep[which(feats.keep==TRUE)[(params$pop.max+1):length(which(feats.keep==TRUE))]] <- FALSE
  }

  # Create a list of which features to replace
  feats.replace <- which(!feats.keep)

  # TODO: Let filtered features become part of new features - tuning parameter
  # Create a list of inclusion probabilities
  marg.probs.use <- c(rep(params$eps, length(F.0)), pmin(pmax(marg.probs, params$eps), (1-params$eps)))

  # Perform the replacements
  if(length(S.t)>params$pop.max)
    feats.replace <- sort(feats.replace,decreasing = T)
  for (i in feats.replace) {
    prev.size <- length(S.t)
    prev.feat.string <- print.feature.ts(S.t[[i]], labels=labels, round = 2)
    if(prev.size>params$pop.max)
    {
      cat("Removed feature", prev.feat.string, "\n")
      #print(length(S.t))
      S.t[[i]] <- NULL
      #print(length(S.t))
      #print(length(marg.probs.use))
    }
    else
    {
      #print(length(S.t))
      #print(length(F.0))
      #out <- gen.feature.ts(c(F.0, S.t), marg.probs.use, data, lookback_window, loglik.alpha, probs, length(F.0), params, verbose)
      S.t[[i]] <- gen.feature.ts(c(F.0, S.t), marg.probs.use, data, data.ts, window_list, loglik.alpha, probs, length(F.0), params, verbose)
      #non_ts[i] <- out$ts.bool
      if (prev.size > length(S.t)) {
        if (verbose) {
          cat("Removed feature", prev.feat.string, "\n")
          cat("Population shrinking, returning.\n")
        }
        return(S.t)
      }
      if (verbose) cat("Replaced feature", prev.feat.string, "with", print.feature.ts(S.t[[i]], labels=labels, round = 2), "\n")
      feats.keep[i] <- T
      marg.probs.use[i] <- mean(marg.probs.use)
    }
  }

  # Add additional features if the population is not at max size
  if (length(S.t) < params$pop.max) {
    for (i in (length(S.t)+1):params$pop.max) {
      prev.size <- length(S.t)
      print(prev.size)
      print(params$pop.max)
      #out <- gen.feature.ts(c(F.0, S.t), marg.probs.use, data, loglik.alpha, probs, length(F.0), params, verbose)
      S.t[[i]] <- gen.feature.ts(c(F.0, S.t), marg.probs.use, data, data.ts, window_list, loglik.alpha, probs, length(F.0), params, verbose)
      #non_ts[i] <- out$ts.bool
      if (prev.size == length(S.t)) {
        if (verbose) cat("Population not growing, returning.\n")
        return(S.t)
      }
      if (verbose) cat("Added feature", print.feature.ts(S.t[[i]], labels=labels, round = 2), "\n")
      marg.probs.use <- c(marg.probs.use, params$eps)
    }
  }
  return(S.t)
}

