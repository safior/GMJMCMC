loglik.pre.re2 <- function (loglik.pi, model, complex, data, params = NULL, visited.models = NULL, sub = FALSE,#visited.models = visited.models, sub = sub, 
                            re_data, re.pop, re.ind) {  

    mod.check <- c(model, re.ind)
    if (!is.null(visited.models) && has_key(visited.models, mod.check)) {
      if (!sub) {
        return(visited.models[[mod.check]])
      } else {
        params$coefs <- visited.models[[mod.check]]$coefs
        params$crit <- visited.models[[mod.check]]$crit
        params$re.mod <- visited.models[[mod.check]]$re.mod
      }
  }

  # Find current random effect
  if(re.ind > 0) {
    re.feat <- re.pop[[re.ind]]
  }
  else if (re.ind == 0) {
    re.feat <- NULL
  }
  # Get the complexity measures for just this model
  complex <- list(width = complex$width[model], oc = complex$oc[model], depth = complex$depth[model])
  # Call the model estimator with the data and the model, note that we add the intercept to every model
  model.res <- loglik.pi(data[, 1], data[, -1], c(T, model), complex, params, re_data, re.feat)
  # Check that the critical value is acceptable
  if (!is.numeric(model.res$crit) || is.nan(model.res$crit)) model.res$crit <- -.Machine$double.xmax
  # Alpha cannot be calculated if the current and proposed models have crit which are -Inf or Inf
  if (is.infinite(model.res$crit)) {
    if (model.res$crit > 0)  model.res$crit <- .Machine$double.xmax
    else model.res$crit <- -.Machine$double.xmax
  }
  return(model.res)
}