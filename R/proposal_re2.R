# Function to generate a proposed model given a current one
# Update of gen.proposal function to accomodate correlation feature.
# Naming of re2 is due to this being the second version.
# Could be applied, with some modifications, to integer valued elements in general.
gen.proposal.re2 <- function (model, params, type, indices=NULL, probs=NULL, prob=FALSE, n_re, re.ind, re.ind.fix) {

  # Split regular features and correlation features
  # Determine if current model includes correlation feature.
  # Adjust model representation to be compatible with binary input expected by gen.proposal
  # Set correlation feature element equal to 1 or 0 depending on whether a correlation feature is
  # included in current model or not
  re.incl <- re.ind > 0
  model <- c(model, re.incl)
  model_length <- length(model)
  # Use indices vector to allow or disallow change of correlation feature
  if (!is.null(indices)) {
    indices[model_length] <- !re.ind.fix
  }

  # If re.ind is not changed, then we set the new re.ind equal to the old re.ind
  re.ind.new <- re.ind

  # Generate proposal of fully binary vector using the base gen.proposal function.
  proposal <- gen.proposal(model, params, type, indices, probs, prob)

  if (type < 5) {
    if(proposal$swap[model_length] != re.incl) {
      re.ind.ops <- 0:n_re
      # Sample from correlation features not included in current model
      # +1 below because we allow the removal of the correlation feature
      re.ind.new <- sample(re.ind.ops[-(re.ind+1)], 1)
      # Adjust probability to account for sampling over remaining correlation features, and no correlation feature.
      proposal$prob <- proposal$prob + log(1/n_re)
    }
  }
  else if (type == 5) {
    # If correlation feature was not included, but is to be included by the proposal.
    if (proposal$swap[model_length]) {
      re.ind.ops <- 1:n_re
      # Sample over correlation feature options
      re.ind.new <- sample(re.ind.ops, 1)
      # Adjust probability to account for sampling over remaining correlation features, and no correlation feature.
      proposal$prob <- proposal$prob + log(1/n_re)
    }
  }
  else if (type == 6) {
    # If correlation feature is to be removed from the proposal
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

# Unused
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
