#' @export lagged
lagged <- function(lag, x) {
    new_end <- length(x) - lag
    return(c(rep(NA, lag), x[1:new_end]))
}

mov_avg <- function(window, x) {
    n <- length(x)
    m_av <- rep(NA, n)
    for (i in (window+1):n) {
        j <- i - window
        m_av[i] <- mean(x[j : i])
    }
    return(m_av)
}

#seq <- 1:100
#print(seq)
#l <- lagged(5, seq)
#print(length(l))
#print(l)
#ma <- mov_avg(5, seq)
#print(length(ma))
#print(ma)