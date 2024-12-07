gen.time.series.feature <- function (
  features,
  F.0.size,
  window_list,
  marg.probs,
  ts.trans.probs,
  max.width,
  max.size,
  ts.trans.priors) {
  
  # What if all features are ts features?
  #d <- features[[F.0.size+i-1]]
  # Find all features that have already gone through a time series transformation
  non_ts <- lapply(features, function(x) x <- is.null(attr(x[[length(x)]], "window")))
  non_ts <- unlist(non_ts)
  #print(lapply(features, function(x) x <- print(x[[length(x)]])))
  features <- features[non_ts]
  #print(lapply(features, function(x) x <- print(x[[length(x)]])))
  # Below if check is unecessary if the create.feature function reaturns NULL itself when all features have already 
  # gone through a time series transformation, but this if check is probably faster either way
  # Will never be null because F.0 is concatenated with current population
  if (sum(non_ts) == 0) {
    print(non_ts)
    print("Null features returned")
    return(NULL)
  }

  marg.probs <- marg.probs[non_ts]
  #feat.count <- sample.int(n = (min(max.width, (length(features)))-1), size = 1)
  #feats <- sample.int(n = length(features), size = feat.count, prob = marg.probs+0.00001)
  feats <- sample.int(n = length(features), size = 1, prob = marg.probs+0.00001)
  trans <- sample.int(n = length(ts.trans.probs), size = 1, prob = ts.trans.probs)
  #print(trans)
  #print(ts.trans.probs)

  #print(lookback_window)
  #print(rep(1/lookback_window, lookback_window))
  window_ts <- unlist(window_list[trans])
  pr <- rep(1/length(window_ts), length(window_ts))
  wind <- sample2(window_ts, size = 1, prob = pr)
  #wind <- sample.int(n = lookback_window, size = 1, prob = rep(1/lookback_window, lookback_window))

  #alphas <- rep(1, length(feats)+1)
  #print(lapply(features[F.0.size+1:length(features)], function(x) x <- print(x[[length(x)]])))

  create.feature.ts(trans, wind, features[feats], ts.trans.priors)#, alphas)
}


# Must add another feature type, and ensure that features that have gone through a time series transformation can not do 
# do so again. Should probably be a vector of booleans.
gen.feature.ts <- function (
  features, 
  marg.probs, 
  data, 
  data.ts,
  #non_ts,
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
    #print(marg.probs)
    #print(length(features))
    if (feat.type == 1) feat <- gen.multiplication(features, marg.probs)
    if (feat.type == 2) feat <- gen.modification(features, marg.probs, probs$trans, probs$trans_priors)
    if (feat.type == 3) feat <- gen.projection(features, marg.probs, probs$trans, params$L, params$max.proj.size, probs$trans_priors)
    if (feat.type == 4) feat <- gen.new(features, F.0.size)
    if (feat.type == 5) feat <- gen.time.series.feature(features, F.0.size, window_list, marg.probs, 
                                                        probs$trans_ts, params$L, params$max.proj.size, probs$trans_priors_ts)
    # Check that the feature is not too wide or deep

    #print(print.feature.ts(feat))
    #print(depth.feature(feat))
    #print(width.feature(feat))
    #print(feat.type)
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
    #print(check.collinearity(feat, feats, F.0.size, data, params$col.check.mock.data))
    #print(feat.ok)
    tries <- tries + 1
    params$eps <- min(params$eps + 0.01, 0.5)
    marg.probs <- pmin(pmax(marg.probs, params$eps), (1 - params$eps))
  }

  # out <- list(
  #   feature = NULL,
  #   ts.bool = FALSE
  # )
  # if (feat.ok) {
  #   out$feature <- feat
  #   if (feat.type == 5) out$ts.bool <- TRUE
  #   }
  # return(out)
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
    # The idea of sampling is to avoid removing features only showing collinearity in the start of the time period
    nr_rows <- min(F.0.size * 5, dim(data)[1])
    #sample <- sample.int(nrow(data), nr_rows)
    #print(sample)
    #mock.data <- check.data(data[sample, ], FALSE)
    #mock.data.ts <- check.data(data.ts[(sample + lookback_window), ], FALSE)
    mock.data <- check.data(data[seq_len(nr_rows), ], FALSE)
    mock.data.ts <- check.data(data.ts[seq_len(nr_rows + lookback_window), ], FALSE)
    #print(mock.data)
    #print(mock.data.ts)
  }
  # Use the mock data to precalc the features
  mock.data.precalc <- precalc.features.ts(mock.data, mock.data.ts, lookback_window, features)
  # Fit a linear model with the mock data precalculated features
  linearmod <- lm(as.data.frame(mock.data.precalc[, -2]))
  # Check if all coefficients were possible to calculate
  if (sum(is.na(linearmod$coefficients)) == 0) return(FALSE)
  else return(TRUE)
}