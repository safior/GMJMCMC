predict.gmjmcmc.ts <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975),  pop = NULL,tol =  0.0000001, ...) {
  transforms.bak <- set.transforms(object$transforms)
  transforms.ts.bak <- set.transforms.ts(object$transforms.ts)

  if (object$add_lagged_response) {
    lag_resp <- c(NA, x[1:(nrow(x)-1), 1])
    x <- cbind(x, lagged_response = lag_resp)
    x <- x[2:nrow(x), ]
  }

  x <- x[, -1]
  
  if(!is.null(attr(object,which = "imputed")))
  {
    df <- data.frame(x)
    na.matr <- data.frame(1*(is.na(df)))
    cm <- colMeans(na.matr)
    na.matr <- na.matr[,attr(object,which = "imputed")]
    names(na.matr) <- paste0("mis_",names(na.matr))
    for (i in which(cm!=0)){
      df[[i]][is.na(df[[i]])] <- median(df[[i]], na.rm = TRUE)
    }
    x <- as.matrix(data.frame(df,na.matr))
    rm(df)
    rm(na.matr)
  } else x <- as.matrix(x)
  
  merged <- merge_results.ts(list(object), data = cbind(1,x), populations = pop,tol = tol, lw = object$lookback_window)
  set.transforms(transforms.bak)
  set.transforms.ts(transforms.ts.bak)
  return(predict.gmjmcmc_merged.ts(merged, x, link, quantiles))
}

#' @export
predict.gmjmcmc_merged.ts <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), pop = NULL,tol =  0.0000001, ...) {
  
  lookback_window <- object$results.raw[[1]]$lookback_window
  if(!is.null(attr(object,which = "imputed")))
  {
    df <- data.frame(x)
    na.matr <- data.frame(1*(is.na(df)))
    cm <- colMeans(na.matr)
    na.matr <- na.matr[,attr(object,which = "imputed")]
    names(na.matr) <- paste0("mis_",names(na.matr))
    for (i in which(cm!=0)){
      df[[i]][is.na(df[[i]])] <- median(df[[i]], na.rm = TRUE)
    }
    x <- as.matrix(data.frame(df,na.matr))
    rm(df)
    rm(na.matr)
  } else {
    #print(lookback_window + 1)
    x2 <- x[(lookback_window + 1) : nrow(x), ]
    x2 <- as.matrix(x2)
    x <- as.matrix(x)
  }

  transforms.ts.bak <- set.transforms.ts(object$results.raw[[1]]$transforms.ts)
  transforms.bak <- set.transforms(object$transforms)
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
        # Models which have 0 weight are skipped since they may also be invalid, and would not influence the predictions.
        if (models[[k]]$crit == -.Machine$double.xmax) next
        yhat[, k] <- link(x.precalc[, c(TRUE, models[[k]]$model), drop=FALSE] %*% models[[k]]$coefs)
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
  set.transforms.ts(transforms.ts.bak)
  set.transforms(transforms.bak)
  return(list(aggr = aggr, preds = preds))
}