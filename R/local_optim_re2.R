# Title     : Local optimization
# Objective : Local optimizers for the mode jumps
# Created by: jonlachmann
# Created on: 2021-02-11

simulated.annealing.re2 <- function (model, data, loglik.pi, indices, complex, params, loglikparams, kernel=NULL, visited.models=NULL, sub = FALSE, 
                                    re_data, re.pop, re.ind, re.ind.fix) {
  # Initialize a list to keep models that we visit in
  models <- vector("list", 0)

  # Select which kernel to use for the random steps
  if (is.null(kernel)) kernel <- sample.int(n = 6, size = 1, prob = params$kern$probs)

  temp <- params$t.init # Initial temperature

  # Calculate current likelihood
  model.res <- loglik.pre.re2(loglik.pi, model, complex, data, loglikparams, visited.models, sub, 
                              re_data = re_data, re.pop = re.pop, re.ind = re.ind)#, re.ind.fix = re.ind.fix)
  model.lik <- model.res$crit
  models[[length(models) + 1]] <- list(prob=NA, model=c(model, re.ind), coefs=model.res$coefs, crit=model.lik, 
                                        alpha=NA, re.mod = model.res$re.mod)
  # print(paste("SA Start:", model.lik))
  n_re <- length(re.pop)
  while (temp > params$t.min) {
    # Make M tries at current temperature
    for (m in 1:params$M) {
      # Get a modified model as proposal and calculate its likelihood
      prop <- gen.proposal.re2(model, params$kern, kernel, indices, n_re = n_re, re.ind = re.ind, re.ind.fix = re.ind.fix)
      prop.re.ind <- prop$re.ind
      proposal <- xor(model, prop$swap)

      model.proposal <- loglik.pre.re2(loglik.pi, proposal, complex, data, loglikparams, visited.models = visited.models, sub = sub, 
                                      re_data, re.pop, prop.re.ind)

      proposal.lik <- model.proposal$crit
      # Append random effects index to proposal
      # Store the model that we have calculated
      models[[length(models) + 1]] <- list(prob=NA, model=c(proposal, re.ind), coefs=model.proposal$coefs, crit=proposal.lik, 
                                          alpha=NA, re.mod = model.proposal$re.mod)
      # Calculate move probability for negative steps (Bolzmann distribution, see Blum and Roli p. 274)
      if (proposal.lik > model.lik) alpha <- 1
      else alpha <- min(1, exp((proposal.lik - model.lik) / temp))
      # Accept move with probability alpha
      if (runif(1) < alpha) {
        model <- proposal
        re.ind <- prop.re.ind
        model.lik <- proposal.lik
      }
    }
    # Update temperature
    temp <- temp * exp(-params$dt)
  }
  # print(paste("SA Finish:", model.lik))
  return(list(model=model, kern=kernel, models=models, re.ind = re.ind))
}

greedy.optim.re2 <- function (model, data, loglik.pi, indices, complex, params, loglikparams, kernel=NULL, visited.models = NULL, sub = FALSE,
                              re_data = re_data, re.pop = re.pop, re.ind = re.ind, re.ind.fix = re.ind.fix) {
  # Initialize a list to keep models that we visit in
  models <- vector("list", 0)

  # Select which kernel to use for the random steps
  if (is.null(kernel)) kernel <- sample.int(n = 6, size = 1, prob = params$kern$probs)

  # Calculate current likelihood
  model.res <- loglik.pre.re2(loglik.pi, model, complex, data, loglikparams, visited.models, sub, 
                              re_data, re.pop, re.ind)
  model.lik <- model.res$crit
  models[[length(models)+1]] <- list(prob=NA, model=c(model, re.ind), coefs=model.res$coefs, crit=model.lik, 
                                      alpha=NA, re.mod = model.res$re.mod)

  n_re <- length(re.pop)
  # Run the algorithm for the number of steps specified

  for (i in 1:params$steps) {
    # For each step, do the specified number of tries
    proposal.best <- NULL
    proposal.lik.best <- -Inf
    for (j in 1:params$tries) {
      # Get a modified model as proposal and calculate its likelihood
      prop <- gen.proposal.re2(model, params$kern, kernel, indices, n_re = n_re, re.ind = re.ind, re.ind.fix=FALSE)
      prop.re.ind <- prop$re.ind
      proposal <- xor(model, prop$swap)
      #print(model)
      #print(prob$swap)

      model.proposal <- loglik.pre.re2(loglik.pi, proposal, complex, data, loglikparams, visited.models, sub, 
                                        re_data, re.pop, prop.re.ind)
      proposal.lik <- model.proposal$crit

      # Append random effects index to proposal
      # Store the model that we have calculated
      # Store random effects below
      models[[length(models)+1]] <- list(prob=NA, model=c(proposal, prop.re.ind), coefs=model.proposal$coefs, crit=proposal.lik, 
                                          alpha=NA, re.mod = model.proposal$re.mod)
      if (proposal.lik > proposal.lik.best) {
        proposal.best <- proposal
        re.ind.best <- prop.re.ind
        proposal.lik.best <- proposal.lik
      }
    }
    # Accept every improvement
    if (proposal.lik.best > model.lik) {
      model <- proposal.best
      re.ind <- re.ind.best
      model.lik <- proposal.lik.best
    }
  }
  return(list(model=model, kern=kernel, models=models, re.ind = re.ind))
}

local.optim.re2 <- function (model, data, loglik.pi, indices, complex, type, params, kernel=NULL, visited.models = NULL, sub = FALSE, 
                            re_data, re.pop, re.ind, re.ind.fix) {
  if (type == 1) {
    return(simulated.annealing.re2(model, data, loglik.pi, indices, complex, params$sa, params$loglik, kernel, visited.models = NULL, sub = FALSE, 
            re_data, re.pop, re.ind, re.ind.fix))
  }
  if (type == 2) {
    return(greedy.optim.re2(model, data, loglik.pi, indices, complex, params$greedy, params$loglik, kernel, visited.models = NULL, sub = FALSE,
            re_data, re.pop, re.ind, re.ind.fix))
  }
  if (type == 3) {
    return("not implemented")
  }
  stop("Invalid local optimizer chosen")
}
