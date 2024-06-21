#' @title Construct standard and adjusted predictions for independent cluster GEEs.
#' 
#' @description
#' For independent cluster GEEs, and prediction at future time points for the given clusters, the standard approach adopted by many practitioners and statistical software is to simply base the prediction on the marginal mean model alone. 
#' On the other hand, adjusted predictions are built by viewing the GEE as solving an iterative working linear model, and then borrowing ideas from universal kriging to construct a predictor that exploits working cross-correlations between the current and new observations within the same cluster. 
#' 
#' @param object A fitted GEE object obtained from applying [geepack::geeglm()].
#' @param newdata New data frame for prediction. It must *contain* a column called \code{id} identifying the clusters.
#' @param intervals Should uncertainty intervals also be produced for the point predictions? Defaults to \code{FALSE}.
#' @param intervals_quantiles Quantiles to use for constructing uncertainty intervals. Default to 95\% intervals formed by taking the lower and upper 2.5\% quantiles. 
#' 
#' @details
#' Uncertainty intervals for standard GEE predictions are constructed in two ways: 1) using the asymptotic normality of the GEE coefficient to construct Wald-type interval on the linear predictor scale, and then inverting the limits; 2) using the asymptotic normality of the GEE coefficient to construct Wald-type interval directly on the response scale. Both are inspired by, but do not explicit use, [glmtoolbox::predict.glmgee()] function. 
#' 
#' Uncertainty intervals for adjusted GEE predictions are constructed by attempting to construct an approximate measure of mean squared error of prediction on the response scale, and then assuming a normal distribution as the basis for constructing the uncertainty interval. We do **not** make any guarantees about the performance of these intervals! Please see the associated manuscript for more details.
#' 
#' For all uncertainty intervals, truncation is applied to the intervals limits as appropriate to ensure they stay within the range of the response.
#'  
#' @note
#' Computationally, the [Matrix::Matrix()] package is not used as this function generally assumes the GEEs are not fitted to dataset with extremely large cluster sizes (meaning it not worth employing the package to avoid the overhead).
#' 
#' @return A data frame with the same number of nrows as \code{newdata}, with the following columns:
#' \item{id: }{Same as \code{newdata$id}.}
#' \item{standard: }{Standard GEE predictions}
#' \item{adjusted: }{Adjusted GEE predictions.}
#' \item{standard_lower/standard_upper: }{If appropriate, uncertainty intervals for the standard GEE predictions. These are constructed using the asymptotic normality of the GEE coefficient to construct Wald-type interval on the linear predictor scale, and then inverting the limits.}
#' \item{standard2_lower/standard2_upper: }{If appropriate, uncertainty intervals for the standard GEE predictions. These are constructed using the asymptotic normality of the GEE coefficient to construct Wald-type interval directly on the response scale.}
#' \item{adjusted_lower/adjusted_upper: }{If appropriate, uncertainty intervals for the adjusted GEE predictions.}
#' 
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>


library(mvtnorm)

gee_predictions <- function(object, 
                            newdata, 
                            intervals = FALSE, 
                            intervals_quantiles = c(0.025,0.975)) {
    
    if(!("geeglm" %in% class(object)) ) {
        stop("Currently only a geeglm object is permitted as input.")
        }
    if(is.null(newdata$id))
        stop("newdata must contain a column called id which identifies the clusters.")
    
    
    ##--------------------------------
    #' ## Standard standard GEE predictions for test data 
    ##--------------------------------
    get_alpha <- qnorm(intervals_quantiles)
    y_test_standard <- y_test_adjusted <- UI_adjusted <- UI_xi <- NULL
    raw_res <- residuals(object, type = "pearson")
    linpred_test <- predict(object, newdata = newdata, type = "link")
    
    if(intervals) {
        #' ## Construct model model matrix for test data
        tt <- delete.response(terms(object))
        m <- model.frame(tt, data = newdata, xlev = object$xlevels)
        newX <- model.matrix(tt, m, contrasts.arg = object$contrasts)
        rm(tt, m)
        
        #' ## Standard GEE prediction uncertainty intervals -- Wald interval on linear predictor and then apply inverse link to the limits
        varhat <- vcov(object)
        se <- as.vector(sqrt(apply(tcrossprod(newX, varhat) * newX, 1, sum)))
        UI_standard <- cbind(newX %*% object$coefficients + get_alpha[1] * se, newX %*% object$coefficients + get_alpha[2] * se)  
        UI_standard <- object$family$linkinv(UI_standard)
        
        
        #' ## Standard GEE prediction uncertainty intervals -- Wald interval on response directly
        se <- se * abs(object$family$mu.eta(linpred_test))
        UI_standard2 <- cbind(object$family$linkinv(linpred_test) + get_alpha[1] * se, object$family$linkinv(linpred_test) + get_alpha[2] * se)  
        if(object$family$family %in% c("Gamma", "poisson", "binomial", "inverse.gaussian")) 
            UI_standard2[UI_standard2[,1] < 0, 1] <- 0
        if(object$family$family == "binomial") 
            UI_standard2[UI_standard2[,2] > 1,2] <- 1
        }

    
    ##--------------------------------
    #' ## Adjusted GEE predictions -- If any predictions are not valid, then revert to the original prediction
    ##--------------------------------
    for(k1 in unique(newdata$id)) {
        sel_train_indices <- which(object$data$id == k1)
        num_train <- length(sel_train_indices)
        sel_test_indices <- which(newdata$id == k1)
        num_test <- length(sel_test_indices)
        if(num_train == 0) 
            stop("Adjusted GEE predictions can not be constructed for IDs that are not in the current data i.e., in object$data.")
        
        y_test_standard <- c(y_test_standard, object$family$linkinv(linpred_test[sel_test_indices]))
        

        #' ## Construct marginal working covariance matrix for some test data
        adjpred_fn <- function(linpred, raw_res, alpha) { # Function can accept vectors for linpred and raw_res, but only a scalar for alpha
            fullmargcorr <- .calc_marginal_correlation(structure = object$corstr, 
                                                      num_train = num_train, 
                                                      num_test = num_test, 
                                                      alpha = alpha)
            trainmargcov_inv <- solve(fullmargcorr[1:num_train, 1:num_train])
                
            rootA_test <- diag(x = sqrt(object$family$variance(object$family$linkinv(linpred))), nrow = num_test)
            adjusted_pred <- as.vector((crossprod(fullmargcorr[1:num_train, num_train + 1:num_test], trainmargcov_inv) %*% raw_res))
                
            make_adjusted_pred <- object$family$linkinv(linpred) + rootA_test %*% adjusted_pred
            check_issues <-  sapply(make_adjusted_pred, object$family$validmu)
            make_adjusted_pred[!check_issues] <- object$family$linkinv(linpred)[!check_issues] 
            
            
            if(intervals) {
                T1 <- diag(rootA_test %*% fullmargcorr[num_train + 1:num_test, num_train + 1:num_test] %*% rootA_test)
                T2 <- diag(rootA_test %*% crossprod(fullmargcorr[1:num_train, num_train + 1:num_test], trainmargcov_inv) %*% fullmargcorr[1:num_train, num_train + 1:num_test] %*% rootA_test)
                
                rootA_train <- diag(x = sqrt(object$family$variance(object$fitted[sel_train_indices])), nrow = num_train)
                Vinv_train <- solve(rootA_train) %*% trainmargcov_inv %*% solve(rootA_train)
                nu_traintest <- rootA_train %*% fullmargcorr[1:num_train, num_train + 1:num_test] %*% rootA_test
                
                D_test <- object$family$mu.eta(linpred) * newX[sel_test_indices,,drop=FALSE]
                D_train <- object$family$mu.eta(object$linear.predictors[sel_train_indices]) * object$geese$X[sel_train_indices,,drop=FALSE]
                tau_traintest <- t(D_test) - crossprod(D_train, Vinv_train) %*% nu_traintest
                rm(D_train, nu_traintest, Vinv_train)
                T3 <- diag(crossprod(tau_traintest, vcov(object)) %*% tau_traintest)
                
                root_MSEP <- sqrt(T1 - T2 + T3)
                xi_test <- D_test %*% object$coefficients + rootA_test %*% adjusted_pred
                UI <- cbind(as.vector(make_adjusted_pred + get_alpha[1] * root_MSEP), 
                          as.vector(make_adjusted_pred + get_alpha[2] * root_MSEP))
            
                if(object$family$family %in% c("Gamma", "poisson", "binomial", "inverse.gaussian")) 
                    UI[UI[,1] < 0, 1] <- 0
                if(object$family$family == "binomial") 
                    UI[UI[,2] > 1, 2] <- 1
                }
            if(!intervals) {
                UI <- UI_xi <- NULL
                }
                
            
            return(list(point_prediction = as.vector(make_adjusted_pred), interval = UI))
            }    
        
        
        out <- try(adjpred_fn(linpred = linpred_test[sel_test_indices], 
                          raw_res = raw_res[sel_train_indices], 
                          alpha = object$geese$alpha), silent = TRUE)
        if(inherits(out, "try-error")) {
            y_test_adjusted <- c(y_test_adjusted, rep(NA, num_test))
            if(intervals) {
                UI_adjusted <- rbind(UI_adjusted, matrix(NA, nrow = num_test, ncol = 2))
                }
            next;
            }
            
        
        y_test_adjusted <- c(y_test_adjusted, out$point_prediction)
        if(intervals) {
            UI_adjusted <- rbind(UI_adjusted, out$interval)
            }
        
        rm(out)
        }
    
    
    if(!intervals)
        out <- data.frame(id = newdata$id, standard = y_test_standard, adjusted = y_test_adjusted)
    if(intervals) {
        out <- data.frame(id = newdata$id, 
                          standard = y_test_standard, 
                          UI_standard = UI_standard, 
                          UI_standard2 = UI_standard2,
                          adjusted = y_test_adjusted, 
                          UI_adjusted = UI_adjusted)
        colnames(out)[3:4] <- c("standard_lower", "standard_upper")
        colnames(out)[5:6] <- c("standard2_lower", "standard2_upper")
        colnames(out)[8:9] <- c("adjusted_lower", "adjusted_upper")
        }

    
    return(out)
    }
