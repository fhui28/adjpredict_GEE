.calc_marginal_correlation <- function(structure, num_train, num_test, alpha = NULL) {
    structure <- match.arg(structure, choices = c("independence", "exchangeable","ar1", "toeplitz"))
    
    if(structure == "independence") {
        fullmargcov <- diag(nrow = num_train + num_test)
        }
    if(structure == "exchangeable") {
        fullmargcov <- diag(nrow = num_train + num_test)
        fullmargcov[lower.tri(fullmargcov)] <- alpha
        fullmargcov[upper.tri(fullmargcov)] <- alpha
        }
    if(structure == "ar1") {
        num_full <- num_train + num_test
        P <- abs(outer(1:num_full, 1:num_full, FUN = "-"))
        fullmargcov <- alpha^P
        }
    if(structure == "toeplitz") {
        if(length(alpha) != (num_train + num_test))
            stop("The length of alpha should be equal to num_train + num_test")
        fullmargcov <- toeplitz(alpha)
        }
    
    return(fullmargcov)
    }



#' @title This function is solely used to process the simulation results for the associated manuscript. Otherwise it can be safely ignored.
make_summary <- function(results_long) {
    out <- results_long %>%     
        group_by(selection_method, type) %>% 
        reframe(RMSE_mean = mean(RMSE, na.rm = TRUE),
                RMSE_sd = sd(RMSE, na.rm = TRUE),
                MAE_mean = mean(MAE, na.rm = TRUE),
                MAE_sd = sd(MAE, na.rm = TRUE),
                pearson_cor_mean = mean(pearson_cor, na.rm = TRUE),
                pearson_cor_sd = sd(pearson_cor, na.rm = TRUE),
                spearman_cor_mean = mean(spearman_cor, na.rm = TRUE),
                spearman_cor_sd = sd(spearman_cor, na.rm = TRUE),
                tjurR2_mean = mean(tjurR2, na.rm = TRUE),
                tjurR2_sd = sd(tjurR2, na.rm = TRUE)
        ) 
    
    return(out)
    }



