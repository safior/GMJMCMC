# Function for precalculating features for a new feature population
precalc.features.ts <- function (data_, data_ts, lookback_window, features) {
  precalc <- matrix(NA, nrow(data_), length(features) + 2)
  precalc[, 1:2] <- data_[, 1:2]
  for (f in seq_along(features)) {
    feature <- features[[f]]
    feature_string <- print.feature.ts(feature, dataset = TRUE)
    if (!is.null(attr(feature[[length(feature)]], "window"))) {
        data <- data_ts
        #print(feature_string)
        calc <- eval(parse(text = feature_string))
        
        start <- lookback_window
        end <- nrow(data_ts)
        calc <- calc[(start+1):end]
    }
    else {
        data <- data_
        calc <- eval(parse(text = feature_string))
    }
    #print(dim(precalc))
    #print(length(calc))
    precalc[, (f + 2)] <- calc
  }
  # Replace any -Inf and Inf values caused by under- or overflow
  precalc <- replace.infinite.data.frame(precalc)
  return(precalc)
}

#' @export set.transforms.ts
set.transforms.ts <- function (transforms) {
  old_transforms <- getOption("gmjmcmc-transformations-ts")
  options("gmjmcmc-transformations-ts" = transforms)
  return(old_transforms)
}