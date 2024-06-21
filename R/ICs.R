#' @title Computer a few information criterion to select the working correlation structure in an independent cluster GEE.
#' 
#' @description
#' This function is a modification of the [geepack::QIC.geeglm()] function to include more criterion for selecting the working correlation structure. Please see the corresponding help file for more information.
#' 
#' @param object A fitted GEE object obtained from applying [geepack::geeglm()].
#' @param ... Optionally more fitted GEE object obtained from applying [geepack::geeglm()].
#' @param tol The tolerance used for matrix inversion.
#' @param env Environment.
#' 
#' @return For each fitted GEE object, the following three information criterion are computed:
#' \item{QIC: }{The quasi information criterion of [https://doi.org/10.1111/j.0006-341X.2001.00120.x].}
#' \item{CIC: }{The correlation information criterion of [https://doi.org/10.1002/sim.3489].}
#' \item{PAC: }{The criterion of [https://doi.org/10.1007/s00362-017-0881-0].}
#' 
#' @note Acknowledgement goes to authors of the \code{geepack} package for the original code!
#'  
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>

library(geepack)

getICs <- function(object, 
                   ..., 
                   tol = .Machine$double.eps, 
                   env = parent.frame()) {
    
    if(!("geeglm" %in% class(object)) ) {
        stop("QIC requires a geeglm object as input.")
        }
    
    # Setup functions
    invert <- if ("MASS" %in% loadedNamespaces()) {
        MASS::ginv
        } 
    else {
        solve
        }
    
    # Create function to make the computations
    computeqic <- function(object) {
        ## Fitted and observed values for quasi likelihood
        mu <- object$fitted.values
        y  <- object$y
        
        ## Quasi Likelihood 
        ## quasi.R <- sum((y*log(mu.R)) - mu.R) # poisson()$dev.resids - scale and weights = 1
        type <- family(object)$family
        quasi <- switch(type,
                        poisson  = sum((y * log(mu)) - mu),
                        gaussian = sum(((y - mu)^2)/-2),
                        binomial = sum(y * log(mu / (1 - mu)) + log(1 - mu)),
                        Gamma    = sum(-y/(mu - log(mu))),
                        stop("Error: distribution not recognized."))
        
        ## Fit model with independence correlation structure
        object$call$corstr <- "independence"
        object$call$zcor <- NULL
        
        model.indep <- eval(object$call, envir=env) 
        
        # Trace term (penalty for model complexity)
        AIinverse <- invert(model.indep$geese$vbeta.naiv, tol=tol)
        Vr <- object$geese$vbeta
        trace <- sum(diag(AIinverse %*% Vr))
        params <- length(coef(object)) # Mean parameters in the model
        
        kpm <- params+length(object$geese$alpha)
        
        
        # QIC
        QIC  <- -2*(quasi - trace)
        
        
        # PAC -- Note this only works if all clusters have the same size
        SSt <- simplify2array(by(as.vector(object$y - object$fitted), 
                                 INDICES = object$data$id, 
                                 FUN = function(x) tcrossprod(x), simplify = TRUE))  
        numerator <- det(apply(SSt, c(1,2), mean))
        rm(SSt)
        Vi_fn <- function(k1) {
            sel_indices <- which(object$data$id == k1)
            rootA_test <- object$family$variance(object$family$linkinv(object$linear.predictors[sel_indices]))
            rootA_test <- diag(x = sqrt(rootA_test), nrow = length(sel_indices))
            
            Vi <- .calc_marginal_correlation(structure = object$corstr, num_train = length(sel_indices), num_test = 0, alpha = object$geese$alpha)
            return(rootA_test %*% Vi %*% rootA_test)
        }
        denominator <- det(Reduce("+", lapply(unique(object$data$id), Vi_fn))/length(unique(object$data$id)))
        PAC <- abs(numerator/denominator - 1) 
        
        
        output <- c(QIC, trace, PAC)
        names(output) <- c("QIC", "CIC", "PAC")
        output
        }
    
    
    if (length(list(...))) {
        # Make the computations
        results <- lapply(list(object, ...), computeqic)
        
        # Check same data size
        check <- sapply(list(object, ...), function(x) {
            length(x$y)
        })
        
        if(any(check != check[1]))
            warning("models are not all fitted to the same number of observations")
        
        # Merge the results together in a data.matrix
        res <- do.call("rbind", results)
        
        # Set the row names corresponding to the models
        Call <- match.call()
        Call$k <- NULL
        row.names(res) <- as.character(Call[-1L])
        res
        } 
    else {
        computeqic(object)
        }
    
    }


#' @title Compute (even) more criterion for selecting the working correlation structure in an independent cluster GEE.
#' 
#' @description
#' Uses the [glmtoolbox::GHYC()], [glmtoolbox::RJC()], [glmtoolbox::AGPC()] and [glmtoolbox::SGPC()] functions to compute these criterion. Please see the corresponding help file for more information.
#' 
#' @param object A fitted GEE object obtained from applying [glmtoolbox::glmgee()].
#' @param ... Optionally more fitted GEE object obtained from applying [glmtoolbox::glmgee()].
#' 
#' @return For each fitted GEE object, the following three information criterion are computed:
#' \item{GHYC: }{The criterion of [https://doi.org/10.1080/03610926.2010.501938].}
#' \item{RJC: }{The Rotnitzky-Jewell criterion, which is explained in more detail in [https://doi.org/10.1198/000313007X245122].}
#' \item{AGPC: }{The Akaike-type penalized Gaussian pseudo-likelihood criterion of [https://doi.org/10.1002/sim.4300].}
#' \item{SGPC: }{The Schwarz-type penalized Gaussian pseudo-likelihood criterion of [https://doi.org/10.1002/sim.4300]. Note in the manuscript, we refer to this as the GBIC or Gaussian pseudo-likelihood BIC.}
#' 
#' @note Acknowledgement goes to authors of the \code{glmtoolbox} package for the original code! 
#'  
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>

getmoreICs <- function(object, ...) {
    
    if(!("glmgee" %in% class(object)) ) {
        stop("Requires glmgee object as input.")
        }
    
    
    out <- data.frame(
        GHYC = glmtoolbox::GHYC(object, ..., verbose = FALSE)$GHYC,
        RJC = glmtoolbox::RJC(object, ..., verbose = FALSE)$RJC,
        AGPC = glmtoolbox::AGPC(object, ..., verbose = FALSE)$AGPC,
        SGPC = glmtoolbox::SGPC(object, ..., verbose = FALSE)$SGPC 
        )
    
    return(out)
    }
    
    

