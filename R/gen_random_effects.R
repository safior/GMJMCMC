# Generate n random effects
gen_rand_effects <- function(re_ops, marg_lik_method, n) {
  if(marg_lik_method == "inla") {
    rand_effects <- gen_rand_effects_inla(re_ops, n)
  }
  else if(marg_lik_method == "nlme") {
    rand_effects <- gen_rand_effects_nlme(re_ops, n)
  }
  else {
    stop("Invalid method for creating random effects! Only inla and nlme are possible.")
  }
  return(rand_effects)
}

# Generate random effects for inla
gen_rand_effects_inla <- function(rand_effects, n) {
  rand_effects <- sample(rand_effects, size = n, replace=FALSE)
  return(as.list(rand_effects))
}

# Generate a single random effect
gen_rand_effects_nlme <- function(re_ops, n) {
  rand_effects <- list()
  # Should make sure all random effects in each population are unique
  for (i in 1:n) {
    rand_prop <- gen_rand_effect_nlme(re_ops)
    rand_effects[[i]] <- rand_prop
  }
  return(rand_effects)
}

# Generate random effect for nlme.
gen_rand_effect_nlme <- function(rand_effects_ops) {
  cor_struct <- create_cor_struct(rand_effects_ops)
  cor_arg <- get_cor_arg(cor_struct)

  group_ops <- rand_effects_ops$group_ops
  random_arg <- get_re_groups(group_ops)
  return(list(cor_arg = cor_arg, random_arg = random_arg))
}

# Generate a cor struct from the possibilites given by the cor_ops paramter. Used with nlme
create_cor_struct <- function(cor_ops) {
    cor_feat <- list()

    cor_structs <- cor_ops$cor_structs
    cs_ind <- sample.int(length(cor_structs), size = 1)
    cor_feat$cor_struct <- cor_structs[cs_ind]

    params <- cor_ops$params[[cs_ind]]
    n_params <- length(params)
    names <- names(params)
    ps <- rep(NULL, n_params)
    for (i in 1:n_params) {
        v <- sample(params[[i]], size = 1)
        ps[i] <- paste(names[i], "=", v, sep = "")
    }
    cor_feat$params <- ps

    cor_feat$cor_group <- cor_ops$cor_group
    cor_feat$cor_data <- cor_ops$cor_data

    return(cor_feat)
}

# Get correlation argument for nlme
get_cor_arg <- function(cor_feat) {
    cor_struct <- cor_feat$cor_struct
    
    pars <- cor_feat$params
    params <- ""
    for (i in 1:length(pars)) {
        params <- paste(params, pars[i], sep = ", ")
    }

    cor_data <- cor_feat$cor_data
    cor_group <- cor_feat$cor_group

    if(cor_group == "") {
      arg_string <- paste(cor_struct, "(", "form = ~", cor_data, params, ")", sep="")
    }
    else {
      arg_string <- paste(cor_struct, "(", "form = ~", cor_data, " | ", cor_group, params, ")", sep="")
    }

    return(arg_string)
}

# Get a random formula for nlme. Currently only works with one formula
get_re_groups <- function(group_ops, covariates) {
    n_group_ops <- length(names(group_ops))
    groups <- c()
    for (i in 1:n_group_ops) {
        g <- sample(group_ops[[i]], size = 1)
        groups <- c(groups, g)
    }
    random_arg <- paste(groups, collapse = "+")
    #random_arg <- paste0("~", groups)
    #print(random_arg)
    return(random_arg)
}