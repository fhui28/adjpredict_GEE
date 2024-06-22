#' @title Main simulation function for comparing different types of predictions for independent cluster data.
#' 
#' @description
#' The simulation compares predictions on test data, considering the following methods: 
#'  1. Independent cluster GEEs fitted to training data, and predictions made using standard and adjusted GEE prediction approaches. Different information criteria are also used to select the working correlation matrix;
#'  2. An oracle independent cluster GEE where all (training and test) data are used and the true marginal correlation matrix structure is assumed to be known. In this case, the predictions are simply the corresponding fitted values;
#'  3. An independent cluster generalized linear mixed models (GLMM) assuming a random intercept and slope for time.
#'
#' @param k0 The seed used to simulate independent clustered data.
#' @param formula_X One-sided formula for the mean of the marginal regression model for generating independent clustered data. Also used as the basis for forming the mean models in the fitted GEE and
#' @param data A data frame from which containing relevant covariate information. Importantly, it **must** contain a column called \code{id} which indexes the cluster of each row in the data fame, and a column called \code{time} which indexes the time point of each row in the dataset.
#' @param fixed_effects A vector of fixed effect coefficients, where the length must be equal to the number of columns of the model matrix implied by \code{formula_X}.
#' @param margcor The marginal correlation matrix used in the marginal regression model simulation. 
#' @param family The distribution of the responses for the model. 
#' @param training_timepoints A vector identifying the observations in \code{data} corresponding to the training set. This is done by looking at \code{data$time} and subsetting to those corresponding to \code{training_timepoints}.
#' @param test_timepoints A vector identifying the observations in \code{data} corresponding to the test set. This is done by looking at \code{data$time} and subsetting to those corresponding to \code{test_timepoints}.
#' @param true_correlation_structure A string identifying the structure of the true marginal correlation matrix as given by \code{margcor}.
#' 
#' @details
#' Please see the associated manuscript for details of the simulation study. Briefly, independent clustered data are generated from a marginal regression model, and each cluster is then split in training (current) and test (new or future) data. Models are then fitted to the current data and predictions constructed to the future set. Predictive performance is then assessed in terms of root-mean-squared-error (RMSE) and mean absolute error (MAE) averaged across time points and clusters, along with pearson and Spearman rank correlation, and Tjur R2 between the true and predicted responses. Note Tjur R2 only makes sense for binary responses
#' 
#' Currently the code only works assuming equal cluster sizes (time points). But the proposed prediction method does work for non-equal cluster sizes.
#' 
#' @return A list containing the following elements:
#' \item{out: }{A summary for the metrics asessing point predictive performance.}
#' \item{out_UI: }{A summary for the metrics asessing uncertainty interval performance.}
#' \item{selected_correlation_structure: }{For each information criterion used for independent cluster GEEs, the working correlation structure that was chosen.}


library(tidyverse)
library(geepack)
library(glmtoolbox)
library(glmmTMB)

sim_fn <- function(k0, 
                   formula_X,
                   data,
                   fixed_effects,
                   margcor,
                   family,
                   training_timepoints, 
                   test_timepoints, 
                   true_correlation_structure) {
    
    message("Dataset ", k0)
    set.seed(k0)
    
    
    ##--------------------------
    #' ## Simulate cluster from independent cluster marginal regression model
    ##--------------------------
    full_dat <- generatedata_marginalreg(formula_X = formula_X,
                                         data = dat,
                                         fixed_effects = fixed_effects,
                                         margcor = margcor,
                                         family = family)
    
    train_dat <- full_dat %>% 
        dplyr::filter(time %in% training_timepoints)
    test_dat <- full_dat %>% 
        dplyr::filter(time %in% test_timepoints)
    
    
    ##--------------------------
    #' ## Method 1: Fit GEE to training data, also selecting the best working correlation from three choices available in geeglm
    ##--------------------------
    correlation_structure_choices <- c("independence", "ar1", "exchangeable")
    fit_gees <- lapply(correlation_structure_choices, function(x) {
        geeglm(update.formula(formula_X, resp ~ .), 
               id = id, 
               data = train_dat,
               family = family, 
               corstr = x, 
               scale.fix = TRUE) 
        })
    names(fit_gees) <- correlation_structure_choices
    pick_best <- getICs(fit_gees[["independence"]], fit_gees[["ar1"]], fit_gees[["exchangeable"]])
    
    
    correlation_structure_choices2 <- c("Independence", "AR-M-dependent(1)", "Exchangeable") 
    fit_gees2 <- lapply(correlation_structure_choices2, function(x) {
        glmtoolbox::glmgee(update.formula(formula_X, resp ~ .), 
               id = id, 
               data = train_dat,
               family = family, 
               corstr = x, 
               scale.fix = TRUE) 
        })
    pick_best <- bind_cols(pick_best, getmoreICs(fit_gees2[[1]], fit_gees2[[2]], fit_gees2[[3]]))
    rm(fit_gees2, correlation_structure_choices2)
    

    ##--------------------------
    #' ## Method 2: Oracle GEE -- Know all data and know true correlation structure
    ##--------------------------
    full_fit <- geeglm(update.formula(formula_X, resp ~ .), 
                       id = id, 
                       data = full_dat, 
                       family = family, 
                       corstr = true_correlation_structure, 
                       scale.fix = TRUE)
    
    
    ##--------------------------
    #' ## Make GEE predictions
    ##--------------------------
    best_predictions_QIC <- gee_predictions(object = fit_gees[[which.min(pick_best[,1])]], newdata = test_dat %>% select(-resp), intervals = FALSE)
    best_predictions_CIC <- gee_predictions(object = fit_gees[[which.min(pick_best[,2])]], newdata = test_dat %>% select(-resp), intervals = FALSE)
    best_predictions_PAC <- gee_predictions(object = fit_gees[[which.min(pick_best[,3])]], newdata = test_dat %>% select(-resp), intervals = FALSE)
    best_predictions_GHYC <- gee_predictions(object = fit_gees[[which.min(pick_best[,4])]], newdata = test_dat %>% select(-resp), intervals = FALSE)
    best_predictions_RJC <- gee_predictions(object = fit_gees[[which.min(pick_best[,5])]], newdata = test_dat %>% select(-resp), intervals = FALSE)
    best_predictions_AGPC <- gee_predictions(object = fit_gees[[which.min(pick_best[,6])]], newdata = test_dat %>% select(-resp), intervals = FALSE)
    best_predictions_SGPC <- gee_predictions(object = fit_gees[[which.min(pick_best[,7])]], newdata = test_dat %>% select(-resp), intervals = TRUE)

    
    ##--------------------------
    #' ## Method 3: GLMM fit
    ##--------------------------
    fit_glmm <- glmmTMB(update.formula(formula_X, resp ~ . + (1 + time|id)), 
                        data = train_dat, 
                        family = family)   
    y_test_glmm <- predict(fit_glmm, newdata = test_dat, type = "response") 
    
    
    
    ##--------------------------
    #' ## Wrangle output
    ##--------------------------
    out <- data.frame(#id = test_dat$id, 
                      timepoints = test_timepoints,
                      standardpred_QIC = best_predictions_QIC$standard, 
                      adjustedpred_QIC = best_predictions_QIC$adjusted, 
                      standardpred_CIC = best_predictions_CIC$standard, 
                      adjustedpred_CIC = best_predictions_CIC$adjusted, 
                      standardpred_PAC = best_predictions_PAC$standard, 
                      adjustedpred_PAC = best_predictions_PAC$adjusted, 
                      standardpred_GHYC = best_predictions_GHYC$standard, 
                      adjustedpred_GHYC = best_predictions_GHYC$adjusted,
                      standardpred_RJC = best_predictions_RJC$standard, 
                      adjustedpred_RJC = best_predictions_RJC$adjusted, 
                      standardpred_AGPC = best_predictions_AGPC$standard, 
                      adjustedpred_AGPC = best_predictions_AGPC$adjusted, 
                      standardpred_SGPC = best_predictions_SGPC$standard, 
                      adjustedpred_SGPC = best_predictions_SGPC$adjusted, 
                      full_gee = full_fit$fitted.values[full_dat$time %in% test_timepoints],
                      glmm = y_test_glmm, 
                      true = test_dat$resp) %>% 
        pivot_longer(-c(true, timepoints), values_to = "predictions", names_to = "method") %>% 
        mutate(method = fct_inorder(method)) %>% 
        group_by(method) %>%
        summarise(RMSE = sqrt(mean((predictions - true)^2)),
                  MAE = mean(abs(predictions - true)),
                  pearson_cor = cor(true, predictions),
                  spearman_cor = cor(true, predictions, method = "spearman"),
                  tjurR2 = mean(predictions[true == 1]) - mean(predictions[true == 0])
                  )

    
    # Uncertainty interval performance
    out_UI <- data.frame(
        timepoints = test_timepoints,
        standardpred_coverage = (best_predictions_SGPC$standard_lower < test_dat$resp & best_predictions_SGPC$standard_upper > test_dat$resp),
        standardpred_width = (best_predictions_SGPC$standard_upper - best_predictions_SGPC$standard_lower),
        standard2pred_coverage = (best_predictions_SGPC$standard2_lower < test_dat$resp & best_predictions_SGPC$standard2_upper > test_dat$resp),
        standard2pred_width = (best_predictions_SGPC$standard2_upper - best_predictions_SGPC$standard2_lower),
        adjustedpred_coverage = (best_predictions_SGPC$adjusted_lower < test_dat$resp & best_predictions_SGPC$adjusted_upper > test_dat$resp),
        adjustedpred_width = (best_predictions_SGPC$adjusted_upper - best_predictions_SGPC$adjusted_lower),
        resp = test_dat$resp) %>%
        select(-c(timepoints, resp)) %>% 
        colMeans(., na.rm = TRUE)

    set.seed(NULL)
    return(list(out = out %>% arrange(method) %>% as.data.frame, 
                out_UI = out_UI,
                selected_correlation_structure = apply(pick_best, 2, which.min)))
    }



