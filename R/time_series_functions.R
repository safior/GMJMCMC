#' @export lagged
# Lagged function (lag/backshift operator)
lagged <- function(lag, x) {
    new_end <- length(x) - lag
    return(c(rep(NA, lag), x[1:new_end]))
}

# moving average function
mov_avg <- function(window, x) {
    n <- length(x)
    m_av <- rep(NA, n)
    for (i in (window+1):n) {
        j <- i - window
        m_av[i] <- mean(x[j : i])
    }
    return(m_av)
}