# Only one time series transformation is allowed per features, below function checks if the feature 
# has already undergone a time series transformation.
is_ts_feature <- function(feature) {
  is_ts <- FALSE
  feat <- feature[[length(feature)]]
  is_ts <- !is.null(attr(feat, "window"))
  if (!is_ts && is.matrix(feat)) {
      for (i in 2:nrow(feat)) {
        # If we have a nested feature, recurse into it
        if (is.list(feature[[feat[i, 2]]])) {
          is_ts <- is_ts_feature(feature[[feat[i, 2]]])
          if (is_ts) {
            return(is_ts)
          }
        }
    }
  }
  return(is_ts)
}

# Generate a time series feature/lookback modification
gen.time.series.feature <- function (
  features,
  F.0.size,
  window_list,
  marg.probs,
  ts.trans.probs,
  max.width,
  max.size,
  ts.trans.priors) {
  
  # Find all features that have already gone through a time series transformation
  ts_features <- lapply(features, function(x) x <- is_ts_feature(x))

  # Indexes of non time series feature/lookback modification
  non_ts <- !unlist(ts_features)

  # Remove ineligible features
  features <- features[non_ts]

  # Below if check is unecessary if the create.feature function reaturns NULL itself when all features have already 
  # gone through a time series transformation, but this if check is probably faster either way
  # Will never be null if F.0 is concatenated with current population
  if (sum(non_ts) == 0) {
    print("Null features returned")
    return(NULL)
  }

  # Retain only marg.probs for eligible features.
  marg.probs <- marg.probs[non_ts]
  # Sample feature
  feats <- sample.int(n = length(features), size = 1, prob = marg.probs+0.00001)
  # Sample times series transformation
  trans <- sample.int(n = length(ts.trans.probs), size = 1, prob = ts.trans.probs)

  # Get lookback options for the sampled time series transformation
  window_ts <- unlist(window_list[trans])
  # Assign uniform probabilites  over the lookback options
  pr <- rep(1/length(window_ts), length(window_ts))
  # Sample the lookback window
  wind <- sample2(window_ts, size = 1, prob = pr)

  # alphas is not currently implemented for time series transformation
  #alphas <- rep(1, length(feats)+1)

  create.feature.ts(trans, wind, features[feats], ts.trans.priors)#, alphas)
}

# Feature generation including lookback moodifiaction
gen.feature.ts <- function(
  features,
  marg.probs,
  data,
  data.ts,
  window_list,
  loglik.alpha, 
  probs, 
  F.0.size, 
  params, 
  verbose = TRUE) {
  tries <- 0
  feat.ok <- F
  lookback_window <- max(unlist(window_list))
  while (!feat.ok && tries < 50) {
    feat.type <- sample.int(n = 5, size = 1, prob = probs$gen)
    if (feat.type == 1) feat <- gen.multiplication(features, marg.probs)
    if (feat.type == 2) feat <- gen.modification(features, marg.probs, probs$trans, probs$trans_priors)
    if (feat.type == 3) feat <- gen.projection(features, marg.probs, probs$trans, params$L, params$max.proj.size, probs$trans_priors)
    if (feat.type == 4) feat <- gen.new(features, F.0.size)
    # New feature type not included in the base algorithm
    if (feat.type == 5) feat <- gen.time.series.feature(features, F.0.size, window_list, marg.probs, 
                                                        probs$trans_ts, params$L, params$max.proj.size, probs$trans_priors_ts)
    # Check that the feature is not too wide or deep

    if (!(depth.feature(feat) > params$D || width.feature(feat) > params$L)) {
      # Generate alphas using the strategy chosen
      if (params$alpha > 0) {
        feat <- gen.alphas(params$alpha, feat, data, loglik.alpha, verbose)
      }
      if (!is.null(feat)) {
        # Check for linear dependence of new the feature
        if (length(features) == F.0.size) feats <- list()
        else feats <- features[(F.0.size + 1):length(features)]
        if (params$check.col) {
          feat.ok <- !check.collinearity.ts(feat, feats, F.0.size, data, data.ts, lookback_window, params$col.check.mock.data)
        }
        else if (!params$check.col)
          feat.ok <- T
      }
    }
    tries <- tries + 1
    params$eps <- min(params$eps + 0.01, 0.5)
    marg.probs <- pmin(pmax(marg.probs, params$eps), (1 - params$eps))
  }
  if (!feat.ok) return(NULL)
  else return(feat)
}

# Check if there is collinearity present in the current set of features
check.collinearity.ts <- function (proposal, features, F.0.size, data, data.ts, lookback_window, mock) {
  # Add the proposal to the feature list for evaluation
  features[[length(features) + 1]] <- proposal
  # Generate mock data to test with (avoiding too costly computations)
  if (mock)
    mock.data <- matrix(c(runif((F.0.size * 2), -100, 100), rep(1, F.0.size * 2),
                        runif((F.0.size * 2) * (F.0.size), -100, 100)), F.0.size * 2, F.0.size + 2)
  else {
    # Increased multiplicator from 2 to 5 to increase the number of features per population. Seemed like the 
    # mjmcmc part of the algorithm produced alot of identical models. Probably, in general, also makes sense 
    # that time series needs larger multiplicator.
    nr_rows <- min(F.0.size * 5, dim(data)[1])
    # The idea of sampling is to avoid removing features only showing collinearity in the start of the time period
    #sample <- sample.int(nrow(data), nr_rows)
    #mock.data <- check.data(data[sample, ], FALSE)
    #mock.data.ts <- check.data(data.ts[(sample + lookback_window), ], FALSE)
    mock.data <- check.data(data[seq_len(nr_rows), ], FALSE)
    mock.data.ts <- check.data(data.ts[seq_len(nr_rows + lookback_window), ], FALSE)
  }
  # Use the mock data to precalc the features
  mock.data.precalc <- precalc.features.ts(mock.data, mock.data.ts, lookback_window, features)
  # Fit a linear model with the mock data precalculated features
  linearmod <- lm(as.data.frame(mock.data.precalc[, -2]))
  # Check if all coefficients were possible to calculate
  if (sum(is.na(linearmod$coefficients)) == 0) return(FALSE)
  else return(TRUE)
}