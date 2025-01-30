gmjmcmc.parallel.ts <- function(runs = 2, cores = getOption("mc.cores", 2L), 
                        merge.options = list(populations = "best", complex.measure = 2, tol = 0.0000001), data_ts, loglik.pi = gaussian.loglik, loglik.alpha = gaussian.loglik.alpha, transforms, ts_transforms, window_list, add_lagged_response, ...) {
    options("gmjmcmc-transformations" = transforms)
    options("gmjmcmc-transformations-ts" = ts_transforms)

    results <- rmclapply(seq_len(runs), 
                args = list(data_ts = data_ts, loglik.pi = loglik.pi, loglik.alpha = loglik.alpha, transforms = transforms, ts_transforms = ts_transforms, window_list = window_list, add_lagged_response = add_lagged_response, ...), 
                mc.cores = cores,
                fun = gmjmcmc.ts)

    class(results) <- "gmjmcmc_parallel"
    lookback_window <- results[[1]]$lookback_window
    #if (add_lagged_response) {
    ##    lag_resp <- c(NA, data_ts[1:(nrow(data_ts)-1), 1])
     #   x <- cbind(x, lagged_response = lag_resp)
     #   x <- x[2:nrow(x), ]
    #}
    #x <- x[, -1]
    merged <- merge_results.ts(results,
                merge.options$populations,
                merge.options$complex.measure,
                merge.options$tol,
                data = cbind(1, data_ts),
                lw = lookback_window)
  return(merged)
}
