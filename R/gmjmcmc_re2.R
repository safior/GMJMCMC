#' @export gmjmcmc.re2
# Main function for using gmjmcmc extension with correlation functions
# The re2 name is because this is the second version.
gmjmcmc.re2 <- function (
  data_ts,
  loglik.pi,
  loglik.alpha = gaussian.loglik.alpha,
  transforms,
  ###############
  ts_transforms, # Time series transforms
  window_list, # List of lookback options for each time series transform, order must match the ts_transforms list
  add_lagged_response = FALSE, # Should probably always be false with correlation features
  ###############
  # Correlation
  ###############
  re_data, # Data used by the correlation features
  rand_effects, # Correlation structures and random effects options
  marg_lik_method, # "inla" or "nlme" currently.
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
  if (add_lagged_response) {
    lag_resp <- c(NA, data_ts[1:(nrow(data_ts)-1), 1])
    data_ts <- cbind(data_ts, lagged_response = lag_resp)
    # Remove start to ensure time series transforms are computable.
    # 2 because of NA introduced above
    data_ts <- data_ts[2:nrow(data_ts), ]
    re_data <- re_data[2:nrow(data_ts), ]
  }
  # Find the largest possible lookback over all transforms
  lookback_window <- max(unlist(window_list))
  # +1 because we need a full window of covariates
  data <- data_ts[(lookback_window + 1) : nrow(data_ts), ]
  re_ncol <- ncol(re_data)
  re_colnames <- colnames(re_data)
  # Adjust re_data as well
  re_data <- re_data[(lookback_window + 1) : nrow(data_ts), ]
  # Below lines ensure that re_data is converted back to a data frame in case it only has one column
  re_data <- as.data.frame(re_data, nrow=nrow(data), ncol=re_ncol)
  colnames(re_data) <- re_colnames

  # Verify that the data is well-formed
  data <- check.data(data, verbose)
  data_ts <- check.data(data_ts, verbose)

  # Generate default probabilities and parameters if there are none supplied.
  if (is.null(probs)) probs <- gen.probs.gmjmcmc.ts(ts_transforms, transforms)
  if (is.null(params)) params <- gen.params.gmjmcmc(data)
  # Extract labels from column names in dataframe
  labels <- get.labels(data, verbose)
  # Set the transformations options
  set.transforms(transforms)
  set.transforms.ts(ts_transforms)
  # Acceptance probability per population
  accept <- vector("list", P)
  accept <- lapply(accept, function (x) x <- 0)
  # A list of populations that have been visited
  S <- vector("list", P)
  S_re <- vector("list", P) 
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

  # Create first population of regular features
  F.0 <- gen.covariates(ncol(data) - 2)
  if (is.null(params$prel.select))
    S[[1]] <- F.0
  else
    S[[1]] <- F.0[params$prel.select]
  
  complex <- complex.features(S[[1]])

  # Create first population of random effects. Empty list given as first argument because no previous population exist.
  S_re[[1]] <- gmjmcmc.transition.re2(list(), params, 0, rand_effects, marg_lik_method, probs)

  ### Main algorithm loop - Iterate over P different populations
  for (p in seq_len(P)) {
    # Set population iteration count
    if (p != P) N <- N.init
    else N <- N.final
    # Precalculate covariates and put them in data.t
    if (length(params$feat$prel.filter) > 0 | p != 1) data.t <- precalc.features.ts(data, data_ts, lookback_window, S[[p]])
    else {
      data.t <- data
    }
    
    # Initialize first model of population for the regular features
    model.cur <- as.logical(rbinom(n = length(S[[p]]), size = 1, prob = 0.5))
    # Initialize first correlation feature of population
    n_re <- length(S_re[[p]])
    re.ind <- sample(0:n_re, size = 1)

    re.pop <- S_re[[p]]

    model.cur.res <- loglik.pre.re2(loglik.pi, model.cur, complex, data.t, params$loglik, NULL, FALSE, re_data, re.pop, re.ind)
    model.cur <- list(prob = 0, model = model.cur, coefs = model.cur.res$coefs, crit = model.cur.res$crit, alpha = 0, 
                      random_effect = re.ind, re.mod = model.cur.res$re.mod)
    best.crit <- model.cur$crit # Reset first best criteria value

    # Run MJMCMC over the population
    if (verbose) print(paste("Population", p, "begin."))
    mjmcmc_res <- mjmcmc.loop.re2(data.t, complex, loglik.pi, model.cur, N, probs, params, sub, verbose, re_data, re.pop, re.ind)
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

    #########################
    # Split marginal probabilities for regular features and correlation features
    n_pop <- length(S[[p]])
    n_total <- length(marg.probs[[p]])
    marg.probs.cov.1 <- marg.probs[[1]][1 : n_pop]
    marg.probs.cov.p <- marg.probs[[p]][1 : n_pop]
    # Last element corresponds to correlation feature
    marg.probs.re.p <- marg.probs[[p]][(n_pop + 1): n_total]

    # Print the marginal posterior distribution of the features after MJMCMC
    if (verbose) {
      cat(paste("\rCurrent best crit:", mjmcmc_res$best.crit, "\n"))
      cat("Feature importance:\n")
      print_dist(marg.probs.cov.p, sapply(S[[p]], print.feature.ts, labels = labels, round = 2), probs$filter)
      ########################
      # Print correlation feature probabilities
      print_dist(marg.probs.re.p, sapply(S_re[[p]], print.feature.re), probs$filter)
    }
    if (params$rescale.large) prev.large <- params$large
    # Generate a new population of features for the next iteration (if this is not the last)
    if (p != P) {
      new_pop <- gmjmcmc.transition.ts(S[[p]], F.0, data, data_ts, ts_transforms, window_list, loglik.alpha, 
                                          marg.probs.cov.1, marg.probs.cov.p, labels, probs, params$feat, verbose)
      S[[p + 1]] <- new_pop$S.t.new
      params$feat <- new_pop$params_feat
      ####################
      # Generate new population of correlation features
      S_re[[p + 1]] <- gmjmcmc.transition.re2(S_re[[p]], params, marg.probs.re.p, rand_effects, marg_lik_method, probs)
      ####################
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
    re_populations = S_re,             # All random effects per population
    marg.lik.method = marg_lik_method, # Method used to calculate marginal likelihood: "inla" or "nlme"
    marg.probs = marg.probs,           # Marginal feature probabilities per population
    model.probs = model.probs,         # Marginal feature probabilities per population
    model.probs.idx = model.probs.idx, # Marginal feature probabilities per population
    best.margs = best.margs,           # Best marginal model probability per population
    accept = accept,                   # Acceptance rate per population
    accept.tot = accept.tot,           # Overall acceptance rate
    best = max(unlist(best.margs)),    # Best marginal model probability throughout the run
    transforms = transforms,           # Transformations used by the model
    transforms.ts = ts_transforms,     # Time series transforms used by the model
    lookback_window = lookback_window, # Largest lookback option from the window_list
    add_lagged_response = add_lagged_response # Whether the lagged response was added as a covariate or not
  )
  attr(results, "class") <- "gmjmcmc"
  return(results)
}


gmjmcmc.transition.re2 <- function(
  S.t.re,
  params,
  marg.probs,
  rand_effects,
  marg_lik_method,
  probs
) {
  # Sample which features are to be kept acoording to their marginal probabilites
  re.keep <- as.logical(rbinom(n = length(marg.probs), size = 1, prob = pmin(marg.probs / probs$filter, 1)))

  # Remove features that are not kept
  null.ind <- which(!(re.keep==1))
  S.t.re[null.ind] <- NULL

  # Fill up population if it is too small
  max_re <- params$max_rand_effects
  # Safety counter
  i <- 0
  while (length(S.t.re) < max_re) {
    re_new <- gen_rand_effects(rand_effects, marg_lik_method, 1)
    # Unlist re_new
    re_new2 <- re_new[[1]]

    if (!in_population(re_new2, S.t.re, marg_lik_method)) {
      #print(re_new2)
      S.t.re <- c(S.t.re, re_new)
    }

    if(i > 100) {
      stop("Check if max_re parameter exceeds the number of unique random effects")
    }
    i <- i + 1
  }
  return(S.t.re)
}

# Checks if the proposed correlation feature is included in the current list of correlation features
in_population <- function(re.new, S.t.re, marg_lik_method) {
  if (marg_lik_method == "nlme") {
    in_pop <- in_population_nlme(re.new, S.t.re)
  }
  else if (marg_lik_method == "inla") {
    in_pop <- in_population_inla(re.new, S.t.re)
  }
  else {
    stop("Invalid method for creating random effects! Only inla and nlme are possible.")
  }
  return(in_pop)
}

# Checks the correlation features in the nlme case
in_population_nlme <- function(re.new, S.t.re) {
  in_pop <- FALSE
  n.pop <- length(S.t.re)
  # If current list of correlation features are empty, then there is nothing to check.
  if (n.pop == 0) {
    return(in_pop)
  }
  # Check if both cor struct and group argument are equal.
  for (i in 1:n.pop) {
    cor.arg.match <- re.new$cor_arg == S.t.re[[i]]$cor_arg
    rand.arg.match <- re.new$random_arg == S.t.re[[i]]$random_arg
    if (cor.arg.match * rand.arg.match) {
      in_pop <- TRUE
      return(in_pop)
    }
  }
  return(in_pop)
}

# Check the correlation features when using inla for marginal likelihood estimation
in_population_inla <- function(re.new, S.t.re) {
  in_pop <- re.new %in% S.t.re
  return(in_pop)
}