# Likelihood via nlme
mixed.model.loglik.nlme <- function (y, x, model, complex, params, re_data, re) {
    # In case the model does not have any correlation feature we assign null arguments
    if (is.null(re)) {
        re <- list(cor_arg = NULL, random_arg = "")
    }
    cor_arg <- re$cor_arg
    cor_arg <- eval(parse(text = cor_arg))
    if (sum(model) > 1) {
        x.model <- x[, as.logical(model)]
        data <- data.frame(y, x = x.model[, -1])
        covariates <- paste0(names(data)[2:dim(data)[2]], collapse = "+")
        formula <- paste0("y ~ 1 + ", covariates)
    } else { 
        #model without fixed effects
        data <- data.frame(y)
        covariates <- ""
        formula <- paste0("y ~ 1")
    }
    formula <- as.formula(formula)

    data <- cbind(data, re_data)

    if (length(params$r) == 0) params$r <- 1/dim(x)[1] # default value or r
    # Log prior
    lp <- log(params$r) * (sum(complex$depth * complex$oc))

    # Returned if model fails.
    noise <- rnorm(1, mean = 0, sd = 10)
    ret_list <- list(crit = -10000 + noise + lp, coefs = rep(0, dim(x)[2]), re.mod = NULL)

    random_arg <- re$random_arg
    # Try-catch block because nlme estimation of specified model might fail, e.g. convergence error
    tryCatch({
        # If no grouping structure
        if (random_arg=="") {
            mm <- gls(formula,
                correlation = cor_arg,
                data = data, 
                method = 'ML')
        }
        # If grouping structure
        else {
            re_formula <- transform_re_formula(random_arg, covariates)
            mm <- lme(formula,
                correlation = cor_arg,
                random = re_formula,
                data = data, 
                method = 'ML')
        }

        # Laplace approximation
        mloglik <- as.numeric(logLik(mm)) - log(length(y)) * (dim(data)[2] - (1+dim(re_data)[2])) / 2
        # In-sample predictions
        preds <- predict(mm, data)
        # To predict out-of-sample the most efficient approach would be to predict here. 
        # The tested alternative is to return the full model, this works but is very memory intensive, 
        # and is commented out

        ret_list <- list(crit = mloglik + lp, coefs = coef(mm), re.mod = preds) #re.mod = mm)
        }, error = function(e) {
            # Print formula an error in case model fitting fails
            print(formula)
            cat("An error occurred:", conditionMessage(e), "\n")
     })
    return(ret_list)
}

# Format grouping structure formula for nlme
transform_re_formula <- function(group_formula, covariates) {
    if (grepl("model_", group_formula, fixed = TRUE)) {
        if (covariates=="") {
            group_formula <- sub("model_", "1", group_formula)
        } else{
            group_formula <- sub("model_", covariates, group_formula)
        }
    }
    group_formula <- paste("~", group_formula)
    return(as.formula(group_formula))
}


# Likelihood via inla
loglik.inla <- function (y, x, model, complex, params, re_data, re) {
    
    if(sum(model)>1) {
        x.model <- x[, as.logical(model)]
        data <- data.frame(y, x = x.model[, -1])
        covariates <- paste0(names(data)[2:dim(data)[2]], collapse = "+")
        formula <- paste0("y ~ 1 + ", covariates)
        # Inla does not accept added columns (I think?), therefore we pre-add them. 
        # This was intended for random slope model, however, the current framework probably results 
        # in overparametrization, and should therefore be improved before it is used. 
        # E.g. a model of the form y ~ x_1 + x_2 + f(group_id, (x_1+x_2), model = iid) 
        model_ <- rowSums(matrix(x.model[, 2:ncol(x.model)], nrow = nrow(x)))
  } else {
        data <- data.frame(y)
        covariates <- ""
        formula <- paste0("y ~ 1")
        model_ <- rep(1, length(y))
    }
    if (!is.null(re)) {
        formula <- paste(formula, "+", re)
    }
    formula <- as.formula(formula)
    data <- cbind(data, re_data, model_)

    if (length(params$r) == 0)  params$r <- 1/dim(x)[1]
    # Log prios
    lp <- log(params$r) * (sum(complex$depth * complex$oc))

    # In case model fitting fails, a default list is returned
    noise <- rnorm(1, mean = 0, sd = 1)
    ret_list <- list(crit = -10000 + noise + lp, coefs = rep(0, dim(x)[2]), re.mod = NULL)
   # Try-catch block in case model fitting fails
   tryCatch({
        # Family is now hardcoded to Gaussian, but could of course be an input parameter
        mod <- inla(data = data, formula = formula, family = "gaussian", 
                    silent = 1L, safe = F)#, verbose = FALSE, debug = FALSE) 
        # Get marginal likelihood from fitted model object
        mloglik <- mod$mlik[1]
    
        # Predictions below are in-sample and out-of-sample combined, because prediction in inla must be done while 
        # fitting the model. Returning full model is VERY memory intensive and is not a good solution, 
        # although functional as long as enough memory is available.
        ret_list <- list(crit = mloglik + lp, coefs = mod$summary.fixed$mode, re.mod = mod$summary.fitted.values) #re.mod = mod)
        }, error = function(e) {
        # Print formula an error in case model fitting fails.
        print(formula)
        cat("An error occurred:", conditionMessage(e), "\n")
   })

    return(ret_list)
}

# Unused
replace_mod_inla <- function(cor, covariates) {
    if (grepl("model_", cor, fixed = TRUE)) {
        cor <- sub("model_", covariates, cor)
    }
    return(cor)
}