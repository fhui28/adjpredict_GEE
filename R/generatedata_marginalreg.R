#' @title Simulate data from a independent cluster marginal regression model. 
#' 
#' @description
#' Uses [simstudy::simstudy()] under the hood to simulate the responses. 
#"
#' @param formula_X One-sided formula for the mean of the marginal regression model.
#' @param data A data frame from which containing relevant covariate information, and, importantly, **must** also contain a column called "id" which indexes the cluster of each row in the data fame.
#' @param fixed_effects A vector of fixed effect coefficients, where the length must be equal to the number of columns of the model matrix implied by \code{formula_X}.
#' @param margcor The marginal correlation matrix used in the marginal regression model simulation. 
#' @param dispersion The dispersion parameter used in the marginal regression model simulation. Currently this is only employed for the Gamma distribution
#' @param family The distribution of the responses for the model. Currently,
#' 
#' @return A data frame of simulated independent cluster responses, the true marginal mean (\code{mean}), and the so-called "true working response" (\code{xi}; this can be safely ignored).
#' 
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>

library(tweedie)
library(bindata)
library(simstudy)

generatedata_marginalreg <- function(formula_X, 
                                     data, 
                                     fixed_effects, 
                                     margcor, 
                                     dispersion = 1, 
                                     family) {
    
    if(!(family$family %in% c("poisson", "binomial", "gaussian", "Gamma"))) 
        stop("Specified family not permitted. Sorry!")
    
    if(is.null(data$id))
        stop("data must contain an id column identifying the cluster id")
    data$id <- as.numeric(factor(data$id))
    
    true_linpred <- as.vector(model.matrix(formula_X, data = data) %*% fixed_effects)
    simdat <- NULL
    for(k0 in 1:length(unique(dat$id))) {
        sel_indices <- which(dat$id == unique(dat$id)[k0])
        
        
        if(response_type$family == "gaussian")
            simdat <- c(simdat, rmvnorm(1, mean = true_linpred[sel_indices], sigma = margcor))
        if(family$family == "binomial") {
            simdat <- c(simdat , bindata::rmvbin(1, margprob = family$linkinv(true_linpred[sel_indices]), sigma = margcor))
            }
        if(family$family == "poisson") {
            simdat <- c(simdat , simstudy::genCorGen(n = 1, 
                                                     nvars = length(sel_indices), 
                                                     params1 = family$linkinv(true_linpred[sel_indices]), 
                                                     dist = "poisson",
                                                     corMatrix = margcor)$X)
            }
        if(family$family == "Gamma") {
            simdat <- c(simdat , simstudy::genCorGen(n = 1, 
                                                     nvars = length(sel_indices), 
                                                     params1 = family$linkinv(true_linpred[sel_indices]), 
                                                     params2 = rep(dispersion, length(sel_indices)),
                                                     dist = "gamma",
                                                     corMatrix = margcor)$X)
            }
    }
    
    xi_true <- (family$mu.eta(true_linpred) * model.matrix(formula_X, data = data)) %*% fixed_effects + (simdat - family$linkinv(true_linpred)) 
    return(data.frame(dat, 
                      resp = simdat, 
                      mean = family$linkinv(true_linpred), 
                      xi = as.vector(xi_true))
           )
    }


