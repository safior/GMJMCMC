#' @export gen.params.gmjmcmc.ts
gen.params.gmjmcmc.re <- function(data) {
  # Get mjmcmc params
  params <- gen.params.gmjmcmc.ts(data)

  # Parameter below should probably be set manually since correlation features do not undergo transformations
  params$max_rand_effects <- min(3, ceiling(ncol(data)/2))

  return(params)
}