# Function to generate a proposed model given a current one
gen.proposal.re2 <- function (model, params, type, indices=NULL, probs=NULL, prob=FALSE, n_re, re.ind, re.ind.fix) {

  # Split features and random effects
  # Determine if current model includes random effect (ittle bit hackish to use existing gen.proposal function) # Set re.ind to 0 to be compatible with binary input expected by gen.proposal
  re.incl <- re.ind > 0
  model <- c(model, re.incl)
  model_length <- length(model)
  # Use index functionality to disallow change of random effect
  if (!is.null(indices)) {
    indices[model_length] <- !re.ind.fix
  }

  # If re.ind is not changed, then we set the new re.ind equal to the old re.ind
  re.ind.new <- re.ind

  proposal <- gen.proposal(model, params, type, indices, probs, prob)

  if (type < 5) {
    if(proposal$swap[model_length] != re.incl) {
      re.ind.ops <- 0:n_re
      # +1 below because we allow removal of the random effect
      re.ind.new <- sample(re.ind.ops[-(re.ind+1)], 1)
      # Unsure if below line is correct
      proposal$prob <- proposal$prob + log(1/n_re)
    }
  }
  else if (type == 5) {
    if (proposal$swap[model_length]) {
      re.ind.ops <- 1:n_re
      re.ind.new <- sample(re.ind.ops, 1)
      # Unsure if below line is correct
      proposal$prob <- proposal$prob + log(1/n_re)
    }
  }
  else if (type == 6) {
    if (!proposal$swap[model_length]) {
      re.ind.new <- 0
    }
  }
  else {
    stop("Something went wrong")
  }

  # Remove random effect index in swap to ensure comptatibility with mjmcmc_re2
  proposal$swap <- proposal$swap[1:(model_length-1)]
  proposal$re.ind <- re.ind.new
  return(proposal)
}

# Calculate the probaility of getting a specified proposal given the current model (i.e. a pdf function)
prob.proposal.re2 <- function (proposal, current, type, params, probs=NULL, n_re) {
  model_length <- length(proposal)

  proposal.re <- proposal[model_length]
  proposal[model_length] <- proposal.re > 0

  current.re <- current[model_length]
  current[model_length] <- current.re > 0

  # Get the difference between the two models
  swaps <- xor(proposal, current)
  if (type < 5) {
    # Prepare parameters for probability calculation
    if (is.null(probs)) probs <- rep(1, length(proposal))
    if (type == 2 || type == 4) {
      params$neigh.min <- params$neigh.size
      params$neigh.max <- params$neigh.size
    }
    prob <- model.proposal.1_4.prob(swaps, probs, params$neigh.size, params$neigh.min, params$neigh.max)
    # If the proposal random effect differs from the current random effect we must adjust the probability
    if (proposal.re != current.re) {
      prob <- prob + log(1/n_re)
    }
  } else if (type == 5) {
    # Generate a proposal of type 5 (addition of a covariate)
    prob <- model.proposal.5_6.prob(current, addition=TRUE)
    # If proposal adds random effect we must adjust the probability
    # Note that it should not be necessary to check that the current model does not include a random effect, but this is done 
    # as a sanity check
    if (proposal[model_length]==1 && current[model_length]==0) {
      prob <- prob + log(1/n_re)
    }
  } else if (type == 6) {
    # Generate a proposal of type 6 (subtraction of a covariate)
    prob <- model.proposal.5_6.prob(current, addition=FALSE)
  }
  return(prob)
}
