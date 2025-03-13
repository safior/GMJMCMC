#estimator function with lme4
mixed.model.loglik.lme4 <- function (y, x, model, complex, params, cor_data, random_effects)
{
    if (sum(model) > 1) {
        x.model = x[,model]
        data <- data.frame(y, x = x.model[,-1], dr = cor_data)
        covariates <- paste0(names(data)[2:(dim(data)[2]-1)], collapse = "+")
        random_effs <- create_random_effect_lme4(cor_data)
        mm <- lmer(as.formula(paste0("y ~", covariates,
                    random_effs)), data = data, REML = FALSE)
    } else { #model without fixed effects
        data <- data.frame(y, dr = params$dr)
        mm <- lmer(as.formula(paste0("y ~", random_effs)), data = data, REML = FALSE)
    }
    # Laplace approximation
    mloglik <- as.numeric(logLik(mm)) - log(length(y)) * (dim(data)[2] - 2)
    # logarithm of model prior
    if (length(params$r) == 0) params$r <- 1/dim(x)[1] # default value or r
    lp <- log.prior(params, complex)

    return(list(crit = mloglik + lp, coefs = fixef(mm)))
}

mixed.model.loglik.nlme <- function (y, x, model, complex, params, re_data, re) {
    if (is.null(re)) {
        re <- list(cor_arg = NULL, random_arg = "")
    }
    #print(re$cor_arg)
    #print(re)
    cor_arg <- re$cor_arg
    #print(cor_arg)
    #print(re)
    cor_arg <- eval(parse(text = cor_arg))
    if (sum(model) > 1) {
        #print(model)
        x.model <- x[, as.logical(model)]
        data <- data.frame(y, x = x.model[, -1])
        covariates <- paste0(names(data)[2:dim(data)[2]], collapse = "+")
        formula <- paste0("y ~ 1 + ", covariates)
    } else { #model without fixed effects
        data <- data.frame(y)
        covariates <- ""
        formula <- paste0("y ~ 1")
    }
    formula <- as.formula(formula)

    data <- cbind(data, re_data)

    if (length(params$r) == 0) params$r <- 1/dim(x)[1] # default value or r
    lp <- log(params$r) * (sum(complex$depth * complex$oc))#log.prior(params, complex)

    # Returned if model fails.
    # Legg til støy under
    noise <- rnorm(1, mean = 0, sd = 1)
    ret_list <- list(crit = -10000 + noise + lp, coefs = rep(0, dim(x)[2]), re.mod = NULL)
    random_arg <- re$random_arg
    tryCatch({
        if (random_arg=="") {
            mm <- gls(formula,
                correlation = cor_arg,
                data = data, 
                method = 'ML')
        }
        else {
            re_formula <- transform_re_formula(random_arg, covariates)
            #print(re_formula)
            mm <- lme(formula,
                correlation = cor_arg,
                random = re_formula,
                data = data, 
                method = 'ML')
        }

        # Laplace approximation
        mloglik <- as.numeric(logLik(mm)) - log(length(y)) * (dim(data)[2] - 1)
        # logarithm of model prior

        ret_list <- list(crit = mloglik + lp, coefs = coef(mm), re.mod = mm)
        }, error = function(e) {
            #print(formula)
            #print(data)
            #print(model)
            #print(x)
            #print(x.model)
            #print(covariates)
            #print(re$cor_arg)
            #print(re$random_arg)
            # You can also print a message or log the error if needed
            cat("An error occurred:", conditionMessage(e), "\n")
     })
    return(ret_list)
}

transform_re_formula <- function(group_formula, covariates) {
    if (grepl("model_", group_formula, fixed = TRUE)) {
        if (covariates=="") {
            group_formula <- sub("model_", "1", group_formula)
        } else{
            group_formula <- sub("model_", covariates, group_formula)
        }
    }
    group_formula <- paste("~", group_formula)
    #print(group_formula)
    return(as.formula(group_formula))
}



mixed.model.loglik.inla <- function (y, x, model, complex, params, re_data, re) {
    
    if(sum(model)>1) {
        x.model <- x[, as.logical(model)]
        data <- data.frame(y, x = x.model[, -1])
        #data1 = data.frame(y, as.matrix(x[,model]), params$dr)
        covariates <- paste0(names(data)[2:dim(data)[2]], collapse = "+")
        formula <- paste0("y ~ 1 + ", covariates)
        # Inla does not accept added columns (I think?), therefore we pre-add them.
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
    #formula_inla <- replace_mod_inla(formula, covariates)
    #formula <- as.formula(formula_inla)
    formula <- as.formula(formula)

    data <- cbind(data, re_data, model_)
    #print(data)

    if (length(params$r) == 0)  params$r <- 1/dim(x)[1]
    lp <- log(params$r) * (sum(complex$depth * complex$oc))

    noise <- rnorm(1, mean = 0, sd = 1)
    ret_list <- list(crit = -10000 + noise + lp, coefs = rep(0, dim(x)[2]), re.mod = NULL)
   #Error handling for unstable libraries that one does not trust 100%
   tryCatch({
     mod <- inla(family = "gaussian", silent = 1L, safe = F, data = data, formula = formula)
     #print(mod$summary.random)
     #print(mod$model.random)
     mloglik <- mod$mlik[1]
     ret_list <- list(crit = mloglik + lp, coefs = mod$summary.fixed$mode, re.mod = mod)
     #print(formula)
     }, error = function(e) {
    #print(formula_inla)
     # You can also print a message or log the error if needed
     cat("An error occurred:", conditionMessage(e), "\n")
   })

    return(ret_list)
}

replace_mod_inla <- function(cor, covariates) {
    #print(cor)
    if (grepl("model_", cor, fixed = TRUE)) {
        #covariates <- paste("\"", covariates, "\"")
        #covariates <- paste("(", covariates, ")")
        cor <- sub("model_", covariates, cor)
    }
    return(cor)
}