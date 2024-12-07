create.feature.ts <- function (transform, window, features, trans.priors, alphas=NULL) {
  # Given no alphas, assume no intercept and unit coefficients
  if (is.null(alphas)) alphas <- c(0, rep(1, length(features)))
  if (length(alphas) != (length(features) + 1)) stop("Invalid alpha/feature count")
  # Calculate the depth, operation count and width of the new feature
  if (transform == 0) {
    depth <- 1 + depth.feature(features[[1]]) + depth.feature(features[[2]])
    oc <- 1 + oc.feature(features[[1]]) + oc.feature(features[[2]])
    width <- 2 + width.feature(features[[1]]) + width.feature(features[[2]])
  }
  else {
    depth <- 0 # Assume 0 depth to find the deepest included feature
    oc <- length(features) - 1 # Every + is an operation, and the outer transform is also one
    oc <- oc + trans.priors[transform]
    width <- length(features) # Width is the number of features
    for (i in 1:(length(features))) {
      locdepth <- depth.feature(features[[i]])
      lococ <- oc.feature(features[[i]])
      locwidth <- width.feature(features[[i]])
      width <- width + locwidth
      oc <- oc + lococ
      if (locdepth > depth) depth <- locdepth
    }
    depth <- depth + 1 # Add 1 to depth to account for the outer transformation
  }

  # Generate the new feature matrix
  newFeature <- list(matrix(c(transform, depth, rep(NA,length(alphas)-2),
                      width, 1:(length(features)), alphas), length(alphas)))
  attr(newFeature[[1]], "oc") <- oc
  attr(newFeature[[1]], "window") <- window
  feature <- append(newFeature, features, 0)
  class(feature) <- "feature"
  #print(feature[[length(feature)]])
  return(feature)
}

#' @export print.feature.ts
print.feature.ts <- function (x, dataset = FALSE, alphas = FALSE, labels = FALSE, round = FALSE, ...) {
  fString <- ""
  feat <- x[[length(x)]]
  # This is a more complex feature
  if (is.matrix(feat)) {
    transforms <- getOption("gmjmcmc-transformations")
    transforms.ts <- getOption("gmjmcmc-transformations-ts")
    if (is.null(transforms)) stop("Please set the gmjmcmc-transformations option to your non-linear functions (see ?set.transforms).")
    # Assume that we are not doing multiplication
    op <- "+"
    # Add the outer transform is there is one
    #print(feat)
    #print(transforms)
    if (feat[1, 1] > 0) {
      if (is.null(attr(feat, "window"))) fString <- paste0(fString, transforms[[feat[1, 1]]], "(")
      else {
        w <- attr(feat, "window")
        fString <- paste0(fString, transforms.ts[[feat[1, 1]]], sprintf("(%s, ", w))
      }
    }
    # If g = 0, we are doing multiplication
    else {
      op <- "*"
      fString <- paste0(fString, "(")
    }
    # If we are printing rounded features for neat output, round all alphas
    if (round) {
      feat[, 3] <- round(feat[, 3], round)
    }
    for (j in seq_len(nrow(feat))) {
      # No plus or multiplication sign on the last one
      if (j == nrow(feat)) op <- ""
      # If this is an intercept just add it in
      if (j == 1 && feat[j, 3] != 0) {
        if (!alphas) fString <- paste0(fString, feat[j, 3], op)
        else fString <- paste0(fString, "?", op)
      }
      # Otherwise this is a feature or covariate, do a recursive conversion
      if (j != 1) {
        # Process alphas, which are only present if there is more than one term in the feature
        # this implies that the feature is not a multiplication (i.e. only one _term_).
        if ((nrow(feat) > 2 || feat[1, 3] != 0) && feat[1, 1] > 0) {
          if (alphas) fString <- paste0(fString, "?*")
          else fString <- paste0(fString, feat[j,3], "*")
        }
        fString <- paste0(fString, print.feature(x[[feat[j, 2]]], dataset, alphas, labels, round), op)
      }
    }
    fString <- paste0(fString, ")")
  }
  # This is a plain covariate
  else if (is.numeric(feat)) {
    if (dataset) fString <- paste0("data[,", feat + 2, "]")
    else if (labels[1] != F) fString <- labels[feat]
    else fString <- paste0("x", feat)
  } else stop("Invalid feature structure")
  return(fString)
}