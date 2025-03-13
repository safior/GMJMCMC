mjmcmc.loop.re2 <- function (data, complex, loglik.pi, model.cur, N, probs, params, sub = FALSE, verbose = TRUE, re_data, re.pop, re.ind) {
  # Acceptance count
  accept <- 0
  # Number of covariates or features (+ 1 because only 1 random effect is allowed?) + number of random effects
  covar_count <- ncol(data) - 2 + 1 #+1 does not work and seems complicated to make work
  # A list of models that have been visited
  models <- vector("list", N)
  # Initialize a vector to contain local opt visited models
  lo.models <- vector("list", 0)
  # Get model without random effect and convert random effect index to logical vector for use in pip_estimate
  mod.without.re <- model.cur$model
  n_re <- length(re.pop)
  re.log <- ind.to.log(re.ind, n_re)
  # Merge fixed effects vector with random effects vector
  model.cur$model <- c(model.cur$model, re.ind)
  # Initialize list for keeping track of unique visited models
  visited.models <- hashmap()
  visited.models[[model.cur$model]] <- list(crit = model.cur$crit, coefs = model.cur$coefs, re.mod = model.cur$re.mod)
  #visited.models <- list(models = matrix(model.cur$model, 1, covar_count), crit = model.cur$crit, count = 1)
  best.crit <- model.cur$crit # Set first best criteria value
  best.coefs <- model.cur$coefs

  progress <- 0
  mcmc_total <- as.numeric(c(mod.without.re, re.log))
  #print(mcmc_total)
  for (i in seq_len(N)) {
    if (verbose && N > 40 && i %% floor(N / 40) == 0) progress <- print_progressbar(progress, 40)

    if (i > params$burn_in) {
      pip_estimate <- mcmc_total / i
      # Select current random effect from pip_estimate vector
      mod.length <- length(model.cur$model)
      #print(mod.length)
      #print("fefef")
      re.ind <- model.cur$model[mod.length]
      cur.inds <- c(1 : (mod.length-1), ((mod.length-1) + re.ind))
      #print(1 : (mod.length-1))
      #print(re.ind)
      #print(cur.inds)
      #print(pip_estimate)
      pip_estimate <- pip_estimate[cur.inds]
      #print(length(pip_estimate))
      #print(pip_estimate)
      #print("ukukukuk")
    }
    else pip_estimate <- rep(1 / covar_count, covar_count)

    proposal <- mjmcmc.prop.re2(data, loglik.pi, model.cur, complex, pip_estimate, probs, params, visited.models, sub = sub, re_data, re.pop)#, re.ind)
    if (proposal$crit > best.crit) {
      best.crit <- proposal$crit
      #print(best.crit)
      if (verbose) cat(paste("\rNew best population crit:", best.crit, "\n"))
    }

    # If we did a large jump and visited models to save
    if (!is.null(proposal$models)) {
      lo.models <- c(lo.models, proposal$models)
      for (mod in seq_along(proposal$models)) {
        visited.models[[proposal$models[[mod]]$model]] <- list(crit = proposal$models[[mod]]$crit, coefs = proposal$models[[mod]]$coefs, re.mod = proposal$models[[mod]]$re.mod)
      }
      proposal$models <- NULL
    }
    visited.models[[proposal$model]] <- list(crit = proposal$crit, coefs = proposal$coefs, re.mod = proposal$re.mod)

    # If we did a large jump and visited models to save
    # if (!is.null(proposal$models)) {
    #   lo.models <- c(lo.models, proposal$models)
    #   # If we are doing subsampling and want to update best mliks
    #   if (sub) {
    #     for (mod in seq_along(proposal$models)) {
    #       # Check if we have seen this model before
    #       mod.idx <- vec_in_mat(visited.models$models[seq_len(visited.models$count), , drop = FALSE], proposal$models[[mod]]$model)
    #       if (mod.idx == 0) {
    #         # If we have not seen the model before, add it
    #         visited.models$count <- visited.models$count + 1
    #         visited.models$crit <- c(visited.models$crit, proposal$models[[mod]]$crit)
    #         visited.models$models <- rbind(visited.models$models, proposal$models[[mod]]$model)
    #       } # This is a model seen before, set the best of the values available
    #       else visited.models$crit[mod.idx] <- max(proposal$models[[mod]]$crit, visited.models$crit[mod.idx])
    #     }
    #   }
    #   proposal$models <- NULL
    # }
    # if (sub) {
    #   # Check if we have seen this model before
    #   mod.idx <- vec_in_mat(visited.models$models[seq_len(visited.models$count), , drop = FALSE], proposal$model)
    #   if (mod.idx == 0) {
    #     # If we have not seen the model before, add it
    #     visited.models$count <- visited.models$count + 1
    #     visited.models$crit <- c(visited.models$crit, proposal$crit)
    #     visited.models$models <- rbind(visited.models$models, proposal$model)
    #   } # This is a model seen before, set the best of the values available
    #   else visited.models$crit[mod.idx] <- max(proposal$crit, visited.models$crit[mod.idx])
    # }

    if (log(runif(1)) <= proposal$alpha) {
      model.cur <- proposal
      accept <- accept + 1
    }
    #print(mcmc_total)
    #print(model.cur$model)
    # Convert random effect index to logical vector for use in pip_estimate
    mod.length <- length(model.cur$model)
    mod.without.re <- model.cur$model[1:(mod.length - 1)]
    mod.log <- c(mod.without.re, ind.to.log(model.cur$model[mod.length], n_re))
    #print("thththth")
    #print(mod.log)
    mcmc_total <- mcmc_total + mod.log
    # Add the current model to the list of visited models
    models[[i]] <- model.cur
  }

  # Calculate and store the marginal inclusion probabilities and the model probabilities
  formatted.models <- reformat.re.models(c(models, lo.models), length(re.pop))
  marg.probs <- marginal.probs.renorm(formatted.models, type = "both")

  return(list(
    models = models,
    accept = accept,
    lo.models = lo.models,
    best.crit = best.crit,
    marg.probs = marg.probs$probs.f,
    model.probs = marg.probs$probs.m,
    model.probs.idx = marg.probs$idx
  ))
}

reformat.re.models <- function(models, n.re) {
  mm.length <- length(models[[1]]$model)#['model'])
  #print(mm.length)
  models <- lapply(models, reformat.re.model, n.re, mm.length)
  return(models)
}

reformat.re.model <- function(model, n.re, mm.length) {
  model.model <- model$model
  #mm.length <- length(model.model)
  re.ind <- model.model[mm.length]
  re.ind.log <- ind.to.log(re.ind, n.re)
  model.model <- c(model.model[1:(mm.length - 1)], re.ind.log)
  model$model <- model.model
  return(model)
}

#' Subalgorithm for generating a proposal and acceptance probability in (G)MJMCMC
#'
#' @param data The data to use in the algorithm
#' @param loglik.pi The the (log) density to explore
#' @param model.cur The current model to make the proposal respective to
#' @param complex The complexity measures used when evaluating the marginal likelihood
#' @param pip_estimate The current posterior inclusion probability estimate, used for proposals
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#' @param visited.models A list of the previously visited models to use when subsampling and avoiding recalculation
#'
#'
#' @noRd
#'
mjmcmc.prop.re2 <- function (data, loglik.pi, model.cur, complex, pip_estimate, probs, params, visited.models=NULL, sub = FALSE, re_data, re.pop) {#}, re.ind) {
  model_length <- length(model.cur$model)
  model.cur.mod <- model.cur$model[1:(model_length-1)]
  re.ind.start <- model.cur$model[model_length]
  n_re <- length(re.pop)
  
  l <- runif(1)
  if (l < probs$large) {
    ### Large jump

    ### Select kernels to use for the large jump
    q.l <- sample.int(n = 4, size = 1, prob = probs$large.kern) # Select large jump kernel
    q.o <- sample.int(n = 2, size = 1, prob = probs$localopt.kern) # Select optimizer function
    q.r <- sample.int(n = 2, size = 1, prob = probs$random.kern) # Select randomization kernel

    # Generate and do large jump
    #print(pip_estimate)
    large.jump <- gen.proposal.re2(model.cur.mod, params$large, q.l, NULL, pip_estimate, n_re = n_re, re.ind = re.ind.start, re.ind.fix = FALSE) # Get the large jump
    chi.0.star <- xor(model.cur.mod, large.jump$swap) # Swap large jump indices
    #print(model.cur.mod)
    #print(large.jump$swap)
    re.ind.lj <- large.jump$re.ind
    # If random effect is changed, do not allow change in optimization and randomization
    re.ind.fix <- !(re.ind.start == re.ind.lj)

    # Optimize to find a mode
    #print(length(model.cur.mod))
    #print(length(large.jump$swap))
    localopt <- local.optim.re2(chi.0.star, data, loglik.pi, !large.jump$swap, complex, q.o, params, re_data = re_data, re.pop = re.pop, re.ind = re.ind.lj, re.ind.fix = re.ind.fix) # Do local optimization
    chi.k.star <- localopt$model
    re.ind.lo <- localopt$re.ind
    #print(re.ind.lo)

    # Randomize around the mode
    #proposal <- gen.proposal.re2(chi.k.star, params$random, q.r, !large.jump$swap, pip_estimate, prob=TRUE, n_re = n_re, re.ind = re.ind.lo, re.ind.fix = re.ind.fix)
    #print((pip_estimate * 0 + 1 - params$random$prob))
    proposal <- gen.proposal.re2(chi.k.star, list(neigh.size = length(pip_estimate), 
                                neigh.min = 1, neigh.max = length(pip_estimate)), q.r, NULL, 
                                (pip_estimate * 0 + 1 - params$random$prob), prob=TRUE,
                                n_re = n_re, re.ind = re.ind.lo, re.ind.fix = re.ind.fix)
    proposal$model <- xor(chi.k.star, proposal$swap)
    re.ind.rand <- proposal$re.ind

    # Do a backwards large jump and add in the kernel used in local optim to use the same for backwards local optim.
    chi.0 <- xor(proposal$model, large.jump$swap)
    # Backwards large jump for random effect is the start re if re is a part of the large jump, if not then the 
    # backwards re is the re after forward randomization
    if (re.ind.fix) {
      re.ind.back <- re.ind.start
    }
    else {
      re.ind.back <- re.ind.rand
    }

    # Do a backwards local optimization
    localopt2 <- local.optim.re2(chi.0, data, loglik.pi, !large.jump$swap, complex, q.o, params, kernel = localopt$kern, re_data = re_data, re.pop = re.pop, re.ind = re.ind.back, re.ind.fix = re.ind.fix)
    # re.ind is last element of chi.k
    chi.k <- localopt2$model
    chi.k.re.ind <- localopt2$re.ind

    ### Calculate acceptance probability
    # Set up the parameters that were used to generate the proposal
    prop.params <- list(neigh.min = params$random$min, neigh.max = params$random$max, neigh.size = proposal$S)

    # Calculate current model probability given proposal
    model.cur$prob <- prob.proposal.re2(c(proposal$model, proposal$re.ind), c(chi.k, chi.k.re.ind), q.r, prop.params, pip_estimate, n_re) # Get probability of gamma given chi.k

    # Store models visited during local optimization
    # re.ind is the lasst element of localopt models
    proposal$models <- c(localopt$models, localopt2$models)
  } else {
    ### Regular MH step
    # Select MH kernel
    q.g <- sample.int(n = 6, size = 1, prob = probs$mh)
    # Generate the proposal
    #print("efefef")
    #print(pip_estimate)
    proposal <- gen.proposal.re2(model.cur.mod, params$mh, q.g, NULL, pip_estimate, prob = TRUE, n_re = n_re, re.ind = re.ind.start, re.ind.fix = FALSE)
    proposal$model <- xor(proposal$swap, model.cur.mod)
    #re.ind.new <- proposal$re.ind

    # Calculate current model probability given proposal
    model.cur$prob <- prob.proposal.re2(c(proposal$model, proposal$re.ind), c(model.cur.mod, re.ind.start), q.g, params$mh, pip_estimate, n_re)
  }
  # Calculate log likelihoods for the proposed model
  proposal.res <- loglik.pre.re2(loglik.pi, proposal$model, complex, data, params$loglik, visited.models=visited.models, sub = sub, 
                                  re_data = re_data, re.pop = re.pop, re.ind = proposal$re.ind)
  proposal$crit <- proposal.res$crit

  # If we are running with subsampling, check the list for a better mlik
  # if (!is.null(visited.models)) {
  #   #print("ererere")
  #   #print(proposal$model)
  #   #print(proposal$re.ind)
  #   #print(model.cur.mod)
  #   #print(re.ind.start)
  #   #print(re.pop)
  #   #print("fgfgfg")
  #   #print(visited.models$models[1:visited.models$count,,drop=FALSE])
  #   mod.idx <- vec_in_mat(visited.models$models[1:visited.models$count,,drop=FALSE], c(proposal$model, proposal$re.ind))
  #   if (mod.idx != 0) proposal$crit <- max(proposal$crit, visited.models$crit[mod.idx])
  # }

  # Calculate acceptance probability for proposed model
  proposal$alpha <- min(0, (proposal$crit + model.cur$prob) - (model.cur$crit + proposal$prob))

  ### Format results and return them
  proposal$swap <- NULL; proposal$S <- NULL
  proposal$coefs <- proposal.res$coefs
  # Append random effect ind to model
  proposal$model <- c(proposal$model, proposal$re.ind)
  proposal$re.mod <- proposal.res$re.mod
  #print(proposal$coefs)
  return(proposal)
}
