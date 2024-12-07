#' @export gen.probs.gmjmcmc
gen.probs.gmjmcmc.ts <- function (ts.transforms, transforms) {
  if (!is.character(ts.transforms))
    stop("The argument transforms must be a character vector specifying the transformations.")

  # Get probs for mjmcmc
  probs <- gen.probs.gmjmcmc(transforms)

  ## Feature generation probabilities
  transcount <- length(ts.transforms)
  trans_ts <- rep(1 / transcount, transcount)  # probability for each different time-series transformation
  trans_priors_ts <- rep(1, transcount)        # Default values assigned to each transformation to be used as "operation count".

  probs$trans_ts <- trans_ts
  probs$trans_priors_ts <- trans_priors_ts

  # Update transformation probabilities
  gen <- c(0.20, 0.20, 0.20, 0.20, 0.20)
  probs$gen <- gen

  return(probs)
}

#' @export gen.params.gmjmcmc.ts
gen.params.gmjmcmc.ts <- function(data) {
  # Get mjmcmc params
  params <- gen.params.mjmcmc(data)

  # Subtract 1 as opposed to 2, since lagged response is added as feature
  ncov <- ncol(data) - 1

  feat_params <- list(D = 5, L = 15,                                # Hard limits on feature complexity
                      alpha = 0,                                    # alpha strategy (0 = None, 1,2,3 = strategies as per Hubin et al.) TODO: Fully Bayesian
                      pop.max = min(100,as.integer(ncov * 1.5)),    # Max features population size
                      keep.org = FALSE,                             # Always keep original covariates in every population
                      prel.filter = 0,                           # Filtration threshold for first population (i.e. filter covariates even if keep.org=TRUE)
                      keep.min = 0.8,                               # Minimum proportion of features to always keep [0,1]
                      eps = 0.05,                                   # Inclusion probability limit for feature generation
                      check.col = TRUE,                             # Whether the colinearity should be checked
                      col.check.mock.data = FALSE,                  # Use mock data when checking for colinearity during feature generation
                      max.proj.size = 15)                           # Maximum projection size
  #print(feat_params$pop.max)
  params$feat <- feat_params
  params$rescale.large <- FALSE
  params$prel.select <- NULL                                        # Specify which covariates to keep in the first population. See Issue #15.

  return(params)
}