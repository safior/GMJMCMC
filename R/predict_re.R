#' @export predict.gmjmcmc.inla
# Prediction with inla model fitting
predict.gmjmcmc.inla <- function(object, data, lookback_window=0, link = function(x) x, 
                                  quantiles = c(0.025, 0.5, 0.975),  pop = NULL, tol =  0.0000001, ...) {
  lw <- object$lookback_window

  x <- format_imputed(object, data)
  if(class(object) == "gmjmcmc") {
    object <- merge_results.re2(list(object), data = x, populations = pop, tol = tol, lw = lw, marg.lik.method = object$marg.lik.method)
  }
  
  n.row <- nrow(x) - lw
  n.pred <- sum(is.na(data[, 1]))

  preds <- list()
  for (i in seq_along(object$results)) {
    preds[[i]] <- list()
    for (j in seq_along(object$results[[i]]$populations)) {
      # Select the models and features to predict from at this iteration
      models <- object$results[[i]]$models[[j]]
      model.probs <- object$results[[i]]$model.probs[[j]]

      yhat <- matrix(0, nrow=n.row, ncol=length(models))

      for (k in seq_along(models)) {
        # Models which have 0 weight are skipped since they may also be invalid, and would not influence the predictions.
        if (models[[k]]$crit == -.Machine$double.xmax) next
        inla.mod <- models[[k]]$re.mod
        #inla.res <- inla.mod$summary.random

        #yhat[, k] <- inla.mod$summary.fitted.values[, 'mean']#inla.mod$summary.fitted.values[(n.row-n.pred+1):n.row, 'mean']
        yhat[, k] <- inla.mod[, 'mean']
      }
      mean.pred <- rowSums(yhat %*% diag(as.numeric(model.probs)))
      pred.quant <- apply(yhat, 1, weighted.quantiles, weights=model.probs, prob=quantiles)

      preds[[i]][[j]] <- list(mean=mean.pred, quantiles=pred.quant, weight=object$results[[i]]$pop.weights[j])
    }
  }

  aggr <- list()
  aggr$mean <- 0 * preds[[1]][[1]]$mean
  aggr$quantiles <- 0 * preds[[1]][[1]]$quantiles
  for (i in seq_along(preds)) {
    for (j in seq_along(preds[[i]])) {
      aggr$mean <- aggr$mean + preds[[i]][[j]]$mean * object$results[[i]]$pop.weights[j]
      aggr$quantiles <- aggr$quantiles + preds[[i]][[j]]$quantiles * object$results[[i]]$pop.weights[j]
    }
  }
  ret.list <- list(aggr = aggr, pred = preds)
  return(ret.list)
}

#' @export predict.gmjmcmc.nlme
predict.gmjmcmc.nlme <- function (object, x, re.data, link = function(x) x, quantiles = c(0.025, 0.5, 0.975),  pop = NULL, tol =  0.0000001, ...) {
  transforms.bak <- set.transforms(object$transforms)
  transforms.ts.bak <- set.transforms.ts(object$transforms.ts)

  x <- format_imputed(object, x)

  merged <- merge_results.re2(list(object), data = x, populations = pop, tol = tol, lw = object$lookback_window, marg.lik.method = object$marg.lik.method)
  set.transforms(transforms.bak)
  set.transforms.ts(transforms.ts.bak)
  return(predict.gmjmcmc_merged.nlme(merged, x, re.data, link, quantiles))
}

#' @export predict.gmjmcmc_merged.nlme
predict.gmjmcmc_merged.nlme <- function (object, x, re.data, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), pop = NULL, tol =  0.0000001, ...) {
  
  lookback_window <- object$results.raw[[1]]$lookback_window
  if (object$results.raw[[1]]$add_lagged_response) {
    lag_resp <- c(NA, x[1:(nrow(x)-1), 1])
    x <- cbind(x, lagged_response = lag_resp)
    x <- x[2:nrow(x), ]
    re.data <- re.data[2:nrow(x), ]
  }
  # Since lagged response can be a covariate we need to remove y here, instead of as input to the 
  # function, as for the predict.gmjmcmc_merged function
  x <- x[, -1]

  # Get data in which all possible features are computable for
  x2 <- x[(lookback_window + 1) : nrow(x), ]
  x <- format_imputed(object, x)
  x2 <- format_imputed(object, x2)

  # Match re.data to covariate data
  re.data <- re.data[(lookback_window + 1) : nrow(x), ]

  transforms.bak <- set.transforms(object$transforms)
  transforms.ts.bak <- set.transforms.ts(object$results.raw[[1]]$transforms.ts)
  if(!is.null(pop))
    object <- merge_results.ts(object$results.raw, pop, 2, tol, data = x, lw = lookback_window)
  
  preds <- list()
  for (i in seq_along(object$results)) {
    preds[[i]] <- list()
    for (j in seq_along(object$results[[i]]$populations)) {
      # Select the models and features to predict from at this iteration
      models <- object$results[[i]]$models[[j]]
      features <- object$results[[i]]$populations[[j]]
      model.probs <- object$results[[i]]$model.probs[[j]]

      # Precalculate the features for the new data (c(0,1...) is because precalc features thinks there is an intercept and y col).
      x.precalc <- precalc.features.ts(cbind(0, 1, x2), cbind(0, 1, x), lookback_window, features)[, -1]

      yhat <- matrix(0, nrow=nrow(x2), ncol=length(models))

      for (k in seq_along(models)) {
        # Models which have 0 weight are skipped since they may be invalid, and would not influence the predictions either way.
        if (models[[k]]$crit == -.Machine$double.xmax) next
        nlme.mod <- models[[k]]$re.mod
        # Skip NULL models
        if (is.null(models[[k]]$re.mod)) next
        # Outcommented code predicts using nlme model object
        #mod.length <- length(models[[k]]$model)
        #mod <- models[[k]]$model[-mod.length]
        #mod.data <- x.precalc[, as.logical(c(TRUE, mod)), drop=FALSE]
        #colnames(mod.data) <- names(nlme.mod$coefficients)
        #new.data <- as.data.frame(cbind(mod.data, re.data))
        yhat[, k] <- nlme.mod #predict(nlme.mod, new.data)
      }

      mean.pred <- rowSums(yhat %*% diag(as.numeric(model.probs)))
      pred.quant <- apply(yhat, 1, weighted.quantiles, weights=model.probs, prob=quantiles)

      preds[[i]][[j]] <- list(mean=mean.pred, quantiles=pred.quant, weight=object$results[[i]]$pop.weights[j])
    }
  }

  aggr <- list()
  aggr$mean <- 0 * preds[[1]][[1]]$mean
  aggr$quantiles <- 0 * preds[[1]][[1]]$quantiles
  for (i in seq_along(preds)) {
    for (j in seq_along(preds[[i]])) {
      aggr$mean <- aggr$mean + preds[[i]][[j]]$mean * object$results[[i]]$pop.weights[j]
      aggr$quantiles <- aggr$quantiles + preds[[i]][[j]]$quantiles * object$results[[i]]$pop.weights[j]
    }
  }
  set.transforms(transforms.bak)
  set.transforms.ts(transforms.ts.bak)
  return(list(aggr = aggr, preds = preds))
}