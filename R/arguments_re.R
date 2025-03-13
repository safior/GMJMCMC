#' @export gen.params.gmjmcmc.ts
gen.params.gmjmcmc.re <- function(data) {
  # Get mjmcmc params
  params <- gen.params.gmjmcmc.ts(data)

  params$max_rand_effects <- min(5, ceiling(ncol(data)/2))

  return(params)
}