#' ---
#' title: Simulation study to assess adjusted predictions in GEEs (Gamma distributed responses).
#' abstract: Currently set up with N = 25 clusters and m = 10 time points, using the first five time points as the current dataset and predicting the second set of five time points (new dataset). The code can be straightforwardly modified to use N = 50 or 100 clusters.
#' author: Francis KC Hui
#' date: Code started Nov 2023
#' ---

rm(list = ls())
library(tidyverse)
library(colorspace)
library(patchwork)
library(mvtnorm)
library(doParallel)
library(abind)
library(foreach)

here::i_am("simulations/gamma_n25m10.R")
library(here)

source(here("R", "generatedata_marginalreg.R"))
source(here("R", "ICs.R"))
source(here("R", "utils.R"))
source(here("R", "predictions.R"))

source(here("simulations", "simulationfunction.R"))

registerDoParallel(cores = detectCores()-3)


##-------------------------
#' # Setting 1: True marginal correlation is AR1
##-------------------------
set.seed(112023)
num_clus <- 25
clus_size <- 10
response_type <- Gamma(link = "log")
dat <- data.frame(id = rep(1:num_clus, each = clus_size), time = 1:clus_size)
true_correlation_structure <- "ar1"


#' ## Generate the X covariates, which stay the same throughout the simulation
partX <- rmvnorm(nrow(dat), sigma = matrix(c(0.5,0.2,0.2,0.5), nrow = 2))
X <- data.frame(treatment = rep(rbinom(num_clus, size = 1, prob = 0.5), each = clus_size),
                x1 = partX[,1],
                x2 = partX[,2])
true_fixed_effects <- c(1, 0.5, 0.1, 0.25, -0.25, -0.1)
dat <- cbind(dat, X)
true_R0 <- .calc_marginal_correlation(structure = true_correlation_structure, 
                                      num_train = 5, 
                                      num_test = 5, 
                                      alpha = 0.5) 
rm(X, partX)



#' ## Do the simulation!
#' Dataset 273, 722 hangs due to problem fitting the GEE (at least on the session information given below)
formula_X <- ~ 1 + treatment*time + x1 + x2 
all_results <- foreach(k0 = (1:1000)[-c(273,722)]) %dopar% sim_fn(k0 = k0,
                                                  formula_X = formula_X,
                                                  fixed_effects = true_fixed_effects,
                                                  data = dat,
                                                  margcor = true_R0,
                                                  family = response_type,
                                                  training_timepoints = 1:5,
                                                  test_timepoints = 6:10,
                                                  true_correlation_structure = true_correlation_structure)
save(all_results, file = here("simulations", "gammaresp_AR1margcor_n25m10.RData"))



#' ## Correlation structure selection results
#correlation_structure_choices <- c("independence", "ar1", "exchangeable")
sapply(all_results, function(x) x$selected_correlation_structure) %>% 
    t %>% 
    apply(., 2, table)
true_correlation_structure


#' ## Prediction results
all_results_long <- do.call(rbind, lapply(all_results, function(x) x$out)) %>% 
    mutate(dataset = rep(1:length(all_results), each = nrow(.)/length(all_results))) %>% 
    mutate(selection_method = rep(rep(c("QIC", "CIC", "PAC", "GYHC", "RJC", "AGPC", "SGPC", "Oracle", "GLMM"), c(2,2,2,2,2,2,2,1,1)), length(all_results))) %>%
    mutate(selection_method = fct_inorder(selection_method)) %>% 
    mutate(type = rep(c(rep(c("Standard","Adjusted"), 7), "Oracle", "GLMM"), length(all_results))) %>% 
    mutate(type = fct_inorder(type)) %>% 
    mutate(dataset = factor(dataset)) %>% 
    dplyr::select(-method) %>% 
    relocate(dataset, selection_method, type) 


make_summary(results_long = all_results_long)


all_results_long <- all_results_long %>% 
    filter(selection_method %in% c("SGPC", "Oracle")) %>% 
    select(-c(MAE, tjurR2)) %>% 
    filter(RMSE < 25) # Remove crazy outliers in the standard GEE prediction


p1 <- ggplot(all_results_long %>% 
                 pivot_longer(RMSE:spearman_cor, names_to = "criteria") %>% 
                 mutate(criteria = fct_recode(criteria, "Pearson correlation" = "pearson_cor", "Spearman correlation" = "spearman_cor")) %>% 
                 mutate(criteria = fct_inorder(criteria)),
             aes(x = type, y = value)) +
    geom_boxplot(notch = FALSE) +
    facet_grid(criteria ~ ., scales = "free") +
    labs(x = "Method", y = "Value", color = "Prediction type", title = "AR(1) marginal correlation") + 
    theme_bw() +
    theme(legend.position = "bottom")

p1

sapply(all_results, function(x) x$out_UI) %>% 
    t %>% 
    summary



##-------------------------
#' # Setting 1: True marginal correlation is exchangeable
##-------------------------
set.seed(112023)
num_clus <- 25
clus_size <- 10
response_type <- Gamma(link = "log")
dat <- data.frame(id = rep(1:num_clus, each = clus_size), time = 1:clus_size)
true_correlation_structure <- "exchangeable"


#' ## Generate the X covariates, which stay the same throughout the simulation
partX <- rmvnorm(nrow(dat), sigma = matrix(c(0.5,0.2,0.2,0.5), nrow = 2))
X <- data.frame(treatment = rep(rbinom(num_clus, size = 1, prob = 0.5), each = clus_size),
                x1 = partX[,1],
                x2 = partX[,2])
true_fixed_effects <- c(1, 0.5, 0.1, 0.25, -0.25, -0.1)
dat <- cbind(dat, X)
true_R0 <- .calc_marginal_correlation(structure = true_correlation_structure, 
                                      num_train = 5, 
                                      num_test = 5, 
                                      alpha = 0.5) 
rm(X, partX)



#' ## Do the simulation!
#' Dataset 273, 722 hangs due to fitting GEE (at least on the session information given below)
formula_X <- ~ 1 + treatment*time + x1 + x2 
all_results <- foreach(k0 = (1:1000)[-c(273,722)]) %dopar% sim_fn(k0 = k0,
                                                     formula_X = formula_X,
                                                     fixed_effects = true_fixed_effects,
                                                     data = dat,
                                                     margcor = true_R0,
                                                     family = response_type,
                                                     training_timepoints = 1:5,
                                                     test_timepoints = 6:10,
                                                     true_correlation_structure = true_correlation_structure)
save(all_results, file = here("simulations", "gammaresp_exchmargcor_n25m10.RData"))




#' ## Correlation structure selection results
#correlation_structure_choices <- c("independence", "ar1", "exchangeable")
sapply(all_results, function(x) x$selected_correlation_structure) %>% 
    t %>% 
    apply(., 2, table)
true_correlation_structure


#' ## Prediction results
all_results_long <- do.call(rbind, lapply(all_results, function(x) x$out)) %>% 
    mutate(dataset = rep(1:length(all_results), each = nrow(.)/length(all_results))) %>% 
    mutate(selection_method = rep(rep(c("QIC", "CIC", "PAC", "GYHC", "RJC", "AGPC", "SGPC", "Oracle", "GLMM"), c(2,2,2,2,2,2,2,1,1)), length(all_results))) %>%
    mutate(selection_method = fct_inorder(selection_method)) %>% 
    mutate(type = rep(c(rep(c("Standard","Adjusted"), 7), "Oracle", "GLMM"), length(all_results))) %>% 
    mutate(type = fct_inorder(type)) %>% 
    mutate(dataset = factor(dataset)) %>% 
    dplyr::select(-method) %>% 
    relocate(dataset, selection_method, type) 


make_summary(results_long = all_results_long)


all_results_long <- all_results_long %>% 
    filter(selection_method %in% c("SGPC", "Oracle")) %>% 
    select(-c(MAE, deviance, tjurR2)) %>% 
    filter(RMSE < 25) # Remove crazy standard GEE predictions


p2 <- ggplot(all_results_long %>% 
                 pivot_longer(RMSE:spearman_cor, names_to = "criteria") %>% 
                 mutate(criteria = fct_recode(criteria, "Pearson correlation" = "pearson_cor", "Spearman correlation" = "spearman_cor")) %>% 
                 mutate(criteria = fct_inorder(criteria)),
             aes(x = type, y = value)) +
    geom_boxplot(notch = FALSE) +
    facet_grid(criteria ~ ., scales = "free") +
    labs(x = "Method", y = "Value", color = "Prediction type", title = "Exchangeable marginal correlation") + 
    theme_bw() +
    theme(legend.position = "bottom")

p2

sapply(all_results, function(x) x$out_UI) %>% 
    t %>% 
    summary



##-------------------------
#' # Setting 2: True correlation is Toeplitz (not contained in the candidate set of working correlation structures)
##-------------------------
set.seed(112023)
num_clus <- 25
clus_size <- 10
response_type <- Gamma(link = "log")
dat <- data.frame(id = rep(1:num_clus, each = clus_size), time = 1:clus_size)
true_correlation_structure <- "toeplitz"


#' ## Generate the X covariates, which stay the same throughout the simulation
partX <- rmvnorm(nrow(dat), sigma = matrix(c(0.5,0.2,0.2,0.5), nrow = 2))
X <- data.frame(treatment = rep(rbinom(num_clus, size = 1, prob = 0.5), each = clus_size),
                x1 = partX[,1],
                x2 = partX[,2])
true_fixed_effects <- c(1, 0.5, 0.1, 0.25, -0.25, -0.1)
dat <- cbind(dat, X)
true_R0 <- .calc_marginal_correlation(structure = true_correlation_structure, 
                                      num_train = 5, 
                                      num_test = 5, 
                                      alpha = c(1,seq(0.7,0.2,length=5),c(0.2,0.2,0.2,0.2))) 
rm(X, partX)


#' ## Do the simulation!
#' Dataset 273 hangs due to fitting GEE (at least on the session information given below)
formula_X <- ~ 1 + treatment*time + x1 + x2 
all_results <- foreach(k0 = (1:1000)[-c(273)]) %dopar% sim_fn(k0 = k0,
                                                     formula_X = formula_X,
                                                     fixed_effects = true_fixed_effects,
                                                     data = dat,
                                                     margcor = true_R0,
                                                     family = response_type,
                                                     training_timepoints = 1:5,
                                                     test_timepoints = 6:10,
                                                     true_correlation_structure = true_correlation_structure)
save(all_results, file = here("simulations", "gammaresp_toepmargcor_n25m10.RData"))



#' ## Correlation structure selection results
#correlation_structure_choices <- c("independence", "ar1", "exchangeable")
sapply(all_results, function(x) x$selected_correlation_structure) %>% 
    t %>% 
    apply(., 2, table)
true_correlation_structure


#' ## Prediction results
all_results_long <- do.call(rbind, lapply(all_results, function(x) x$out)) %>% 
    mutate(dataset = rep(1:length(all_results), each = nrow(.)/length(all_results))) %>% 
    mutate(selection_method = rep(rep(c("QIC", "CIC", "PAC", "GYHC", "RJC", "AGPC", "SGPC", "Oracle", "GLMM"), c(2,2,2,2,2,2,2,1,1)), length(all_results))) %>%
    mutate(selection_method = fct_inorder(selection_method)) %>% 
    mutate(type = rep(c(rep(c("Standard","Adjusted"), 7), "Oracle", "GLMM"), length(all_results))) %>% 
    mutate(type = fct_inorder(type)) %>% 
    mutate(dataset = factor(dataset)) %>% 
    dplyr::select(-method) %>% 
    filter(type != "Oracle") %>% 
    relocate(dataset, selection_method, type) 


make_summary(results_long = all_results_long)


all_results_long <- all_results_long %>% 
    filter(selection_method %in% c("SGPC")) %>% 
    select(-c(MAE, tjurR2)) %>% 
    filter(RMSE < 25) # Remove crazy standard GEE predictions  


p3 <- ggplot(all_results_long %>% 
                 pivot_longer(RMSE:spearman_cor, names_to = "criteria") %>% 
                 mutate(criteria = fct_recode(criteria, "Pearson correlation" = "pearson_cor", "Spearman correlation" = "spearman_cor")) %>% 
                 mutate(criteria = fct_inorder(criteria)),
             aes(x = type, y = value)) +
    geom_boxplot(notch = FALSE) +
    facet_grid(criteria ~ ., scales = "free") +
    labs(x = "Method", y = "Value", color = "Prediction type", title = "Toeplitz marginal correlation") + 
    theme_bw() +
    theme(legend.position = "bottom")

p3

sapply(all_results, function(x) x$out_UI) %>% 
    t %>% 
    summary



##-------------------------
#' # Example 2: True correlation is unstructured (not contained in the candidate set of working correlation structures)
##-------------------------
set.seed(112023)
num_clus <- 25
clus_size <- 10
response_type <- Gamma(link = "log")
dat <- data.frame(id = rep(1:num_clus, each = clus_size), time = 1:clus_size)
true_correlation_structure <- "unstructured"


#' ## Generate the X covariates, which stay the same throughout the simulation
partX <- rmvnorm(nrow(dat), sigma = matrix(c(0.5,0.2,0.2,0.5), nrow = 2))
X <- data.frame(treatment = rep(rbinom(num_clus, size = 1, prob = 0.5), each = clus_size),
                x1 = partX[,1],
                x2 = partX[,2])
true_fixed_effects <- c(1, 0.5, -0.5, 0.25, -0.25, -0.1)
dat <- cbind(dat, X)
SigmaR0 <- .calc_marginal_correlation(structure = "ar1",
                                      num_train = 5,
                                      num_test = 5,
                                      alpha = 0.5)
true_R0 <- rWishart(1, df = 10, Sigma = SigmaR0)[,,1] %>% 
    cov2cor
rm(X, partX)




#' ## Do the simulation!
formula_X <- ~ 1 + treatment*time + x1 + x2 
all_results <- foreach(k0 = (1:1000)) %dopar% sim_fn(k0 = k0,
                                                     formula_X = formula_X,
                                                     fixed_effects = true_fixed_effects,
                                                     data = dat,
                                                     margcor = true_R0,
                                                     family = response_type,
                                                     training_timepoints = 1:5,
                                                     test_timepoints = 6:10,
                                                     true_correlation_structure = "ar1")
save(all_results, file = here("simulations", "gaussianresp_unstrucmargcor_n25m10.RData"))



#' ## Correlation structure selection results
#correlation_structure_choices <- c("independence", "ar1", "exchangeable")
sapply(all_results, function(x) x$selected_correlation_structure) %>% 
    t %>% 
    apply(., 2, table)


#' ## Prediction results
all_results_long <- do.call(rbind, lapply(all_results, function(x) x$out)) %>% 
    mutate(dataset = rep(1:length(all_results), each = nrow(.)/length(all_results))) %>% 
    mutate(selection_method = rep(rep(c("QIC", "CIC", "PAC", "GYHC", "RJC", "AGPC", "SGPC", "Oracle", "GLMM"), c(2,2,2,2,2,2,2,1,1)), length(all_results))) %>%
    mutate(selection_method = fct_inorder(selection_method)) %>% 
    mutate(type = rep(c(rep(c("Standard","Adjusted"), 7), "Oracle", "GLMM"), length(all_results))) %>% 
    mutate(type = fct_inorder(type)) %>% 
    mutate(dataset = factor(dataset)) %>% 
    dplyr::select(-method) %>% 
    filter(type != "Oracle") %>%  # Sine it is assumed we can not fit this structure
    relocate(dataset, selection_method, type) 


make_summary(results_long = all_results_long)


all_results_long <- all_results_long %>% 
    filter(selection_method %in% c("SGPC")) %>% 
    select(-c(MAE, tjurR2)) 


p4 <- ggplot(all_results_long %>% 
                 pivot_longer(RMSE:spearman_cor, names_to = "criteria") %>% 
                 mutate(criteria = fct_recode(criteria, "Pearson correlation" = "pearson_cor", "Spearman correlation" = "spearman_cor")) %>% 
                 mutate(criteria = fct_inorder(criteria)),
             aes(x = type, y = value)) +
    geom_boxplot(notch = FALSE) +
    facet_grid(criteria ~ ., scales = "free") +
    labs(x = "Method", y = "Value", color = "Prediction type", title = "Unstructured marginal correlation") + 
    theme_bw() +
    theme(legend.position = "bottom")

p4


sapply(all_results, function(x) x$out_UI) %>% 
    t %>% 
    summary







set.seed(112023)
num_clus <- 25
clus_size <- 10
response_type <- Gamma(link = "log")
dat <- data.frame(id = rep(1:num_clus, each = clus_size), time = 1:clus_size)
true_correlation_structure <- "unstructured"


#' ## Generate the X covariates, which stay the same throughout the simulation
partX <- rmvnorm(nrow(dat), sigma = matrix(c(0.5,0.2,0.2,0.5), nrow = 2))
X <- data.frame(treatment = rep(rbinom(num_clus, size = 1, prob = 0.5), each = clus_size),
                x1 = partX[,1],
                x2 = partX[,2])
true_fixed_effects <- c(1, 0.5, 0.1, 0.25, -0.25, -0.1)
dat <- cbind(dat, X)
SigmaR0 <- .calc_marginal_correlation(structure = "ar1",
                                      num_train = 5,
                                      num_test = 5,
                                      alpha = 0.5)
true_R0 <- rWishart(1, df = 10, Sigma = SigmaR0)[,,1] %>% 
    cov2cor
rm(X, partX)


#' ## Do the simulation!
#' Dataset 273 hangs due to fitting GEE (at least on the session information given below) 
formula_X <- ~ 1 + treatment*time + x1 + x2 
all_results <- foreach(k0 = (1:1000)[-c(273)]) %dopar% sim_fn(k0 = k0,
                                                     formula_X = formula_X,
                                                     fixed_effects = true_fixed_effects,
                                                     data = dat,
                                                     margcor = true_R0,
                                                     family = response_type,
                                                     training_timepoints = 1:5,
                                                     test_timepoints = 6:10,
                                                     true_correlation_structure = "ar1")
save(all_results, file = here("simulations", "gammaresp_unstrucmargcor_n25m10.RData"))



#' ## Correlation structure selection results
#correlation_structure_choices <- c("independence", "ar1", "exchangeable")
sapply(all_results, function(x) x$selected_correlation_structure) %>% 
    t %>% 
    apply(., 2, table)
true_correlation_structure


#' ## Prediction results
all_results_long <- do.call(rbind, lapply(all_results, function(x) x$out)) %>% 
    mutate(dataset = rep(1:length(all_results), each = nrow(.)/length(all_results))) %>% 
    mutate(selection_method = rep(rep(c("QIC", "CIC", "PAC", "GYHC", "RJC", "AGPC", "SGPC", "Oracle", "GLMM"), c(2,2,2,2,2,2,2,1,1)), length(all_results))) %>%
    mutate(selection_method = fct_inorder(selection_method)) %>% 
    mutate(type = rep(c(rep(c("Standard","Adjusted"), 7), "Oracle", "GLMM"), length(all_results))) %>% 
    mutate(type = fct_inorder(type)) %>% 
    mutate(dataset = factor(dataset)) %>% 
    dplyr::select(-method) %>% 
    filter(type != "Oracle") %>% 
    relocate(dataset, selection_method, type) 


make_summary(results_long = all_results_long)


all_results_long <- all_results_long %>% 
    filter(selection_method %in% c("SGPC")) %>% 
    select(-c(MAE, deviance, tjurR2)) %>% 
    filter(RMSE < 30) # Remove crazy standard GEE predictions


p4 <- ggplot(all_results_long %>% 
                 pivot_longer(RMSE:spearman_cor, names_to = "criteria") %>% 
                 mutate(criteria = fct_recode(criteria, "Pearson correlation" = "pearson_cor", "Spearman correlation" = "spearman_cor")) %>% 
                 mutate(criteria = fct_inorder(criteria)),
             aes(x = type, y = value)) +
    geom_boxplot(notch = FALSE) +
    facet_grid(criteria ~ ., scales = "free") +
    labs(x = "Method", y = "Value", color = "Prediction type", title = "Unstructured marginal correlation") + 
    theme_bw() +
    theme(legend.position = "bottom")

p4

sapply(all_results, function(x) x$out_UI) %>% 
    t %>% 
    summary



##----------------------
sessioninfo::session_info()
##----------------------
# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 (2021-11-01)
# os       Linux Mint 21.1
# system   x86_64, linux-gnu
# ui       RStudio
# language en_AU:en
# collate  en_AU.UTF-8
# ctype    en_AU.UTF-8
# tz       Australia/Sydney
# date     2024-06-23
# rstudio  2024.04.0+735 Chocolate Cosmos (desktop)
# pandoc   NA
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package       * version    date (UTC) lib source
# abind         * 1.4-5      2016-07-21 [1] CRAN (R 4.1.2)
# backports       1.5.0      2024-05-23 [1] CRAN (R 4.1.2)
# bigmemory       4.6.1      2022-05-02 [1] CRAN (R 4.1.2)
# bigmemory.sri   0.1.6      2022-11-09 [1] CRAN (R 4.1.2)
# bindata       * 0.9-20     2021-01-29 [1] CRAN (R 4.1.2)
# boot            1.3-28     2021-05-03 [4] CRAN (R 4.1.1)
# broom           1.0.5      2023-06-09 [1] CRAN (R 4.1.2)
# class           7.3-20     2022-01-13 [4] CRAN (R 4.1.2)
# cli             3.6.2      2023-12-11 [1] CRAN (R 4.1.2)
# coda            0.19-4.1   2024-01-31 [1] CRAN (R 4.1.2)
# codetools       0.2-18     2020-11-04 [4] CRAN (R 4.0.3)
# colorspace    * 2.1-0      2023-01-23 [1] CRAN (R 4.1.2)
# data.table      1.14.8     2023-02-17 [1] CRAN (R 4.1.2)
# doParallel    * 1.0.17     2022-02-07 [1] CRAN (R 4.1.2)
# dplyr         * 1.1.4      2023-11-17 [1] CRAN (R 4.1.2)
# e1071           1.7-13     2023-02-01 [1] CRAN (R 4.1.2)
# emmeans         1.9.0      2023-12-18 [1] CRAN (R 4.1.2)
# estimability    1.4.1      2022-08-05 [1] CRAN (R 4.1.2)
# fansi           1.0.6      2023-12-08 [1] CRAN (R 4.1.2)
# farver          2.1.2      2024-05-13 [1] CRAN (R 4.1.2)
# fastglm         0.0.3      2022-05-23 [1] CRAN (R 4.1.2)
# forcats       * 1.0.0      2023-01-29 [1] CRAN (R 4.1.2)
# foreach       * 1.5.2      2022-02-02 [1] CRAN (R 4.1.2)
# Formula         1.2-5      2023-02-24 [1] CRAN (R 4.1.2)
# geepack       * 1.3.9      2022-08-16 [1] CRAN (R 4.1.2)
# generics        0.1.3      2022-07-05 [1] CRAN (R 4.1.2)
# ggplot2       * 3.5.1      2024-04-23 [1] CRAN (R 4.1.2)
# glmmTMB       * 1.1.9      2024-03-20 [1] CRAN (R 4.1.2)
# glmtoolbox    * 0.1.9      2023-10-10 [1] CRAN (R 4.1.2)
# glue            1.7.0      2024-01-09 [1] CRAN (R 4.1.2)
# gtable          0.3.5      2024-04-22 [1] CRAN (R 4.1.2)
# here          * 1.0.1      2020-12-13 [1] CRAN (R 4.1.2)
# hms             1.1.3      2023-03-21 [1] CRAN (R 4.1.2)
# iterators     * 1.0.14     2022-02-05 [1] CRAN (R 4.1.2)
# labeling        0.4.3      2023-08-29 [1] CRAN (R 4.1.2)
# lattice         0.20-45    2021-09-22 [4] CRAN (R 4.1.1)
# lifecycle       1.0.4      2023-11-07 [1] CRAN (R 4.1.2)
# lme4            1.1-35.3   2024-04-16 [1] CRAN (R 4.1.2)
# lubridate     * 1.9.2      2023-02-10 [1] CRAN (R 4.1.2)
# magrittr        2.0.3      2022-03-30 [1] CRAN (R 4.1.2)
# MASS            7.3-55     2022-01-13 [4] CRAN (R 4.1.2)
# Matrix          1.6-4      2023-11-30 [1] CRAN (R 4.1.2)
# mgcv            1.9-0      2023-07-11 [1] CRAN (R 4.1.2)
# minqa           1.2.5      2022-10-19 [1] CRAN (R 4.1.2)
# multcomp        1.4-23     2023-03-09 [1] CRAN (R 4.1.2)
# munsell         0.5.1      2024-04-01 [1] CRAN (R 4.1.2)
# mvtnorm       * 1.2-5      2024-05-21 [1] CRAN (R 4.1.2)
# nlme            3.1-155    2022-01-13 [4] CRAN (R 4.1.2)
# nloptr          2.0.3      2022-05-26 [1] CRAN (R 4.1.2)
# numDeriv        2016.8-1.1 2019-06-06 [1] CRAN (R 4.1.2)
# patchwork     * 1.2.0      2024-01-08 [1] CRAN (R 4.1.2)
# pillar          1.9.0      2023-03-22 [1] CRAN (R 4.1.2)
# pkgconfig       2.0.3      2019-09-22 [1] CRAN (R 4.1.2)
# proxy           0.4-27     2022-06-09 [1] CRAN (R 4.1.2)
# purrr         * 1.0.2      2023-08-10 [1] CRAN (R 4.1.2)
# R6              2.5.1      2021-08-19 [1] CRAN (R 4.1.2)
# Rcpp            1.0.12     2024-01-09 [1] CRAN (R 4.1.2)
# RcppZiggurat    0.1.6      2020-10-20 [1] CRAN (R 4.1.2)
# readr         * 2.1.4      2023-02-10 [1] CRAN (R 4.1.2)
# Rfast           2.0.8      2023-07-03 [1] CRAN (R 4.1.2)
# rlang           1.1.3      2024-01-10 [1] CRAN (R 4.1.2)
# rprojroot       2.0.4      2023-11-05 [1] CRAN (R 4.1.2)
# rstudioapi      0.16.0     2024-03-24 [1] CRAN (R 4.1.2)
# sandwich        3.0-2      2022-06-15 [1] CRAN (R 4.1.2)
# scales          1.3.0      2023-11-28 [1] CRAN (R 4.1.2)
# sessioninfo     1.2.2      2021-12-06 [1] CRAN (R 4.1.2)
# simstudy      * 0.7.1      2023-11-23 [1] CRAN (R 4.1.2)
# statmod         1.5.0      2023-01-06 [1] CRAN (R 4.1.2)
# stringi         1.8.4      2024-05-06 [1] CRAN (R 4.1.2)
# stringr       * 1.5.1      2023-11-14 [1] CRAN (R 4.1.2)
# survival        3.2-13     2021-08-24 [4] CRAN (R 4.1.1)
# TH.data         1.1-2      2023-04-17 [1] CRAN (R 4.1.2)
# tibble        * 3.2.1      2023-03-20 [1] CRAN (R 4.1.2)
# tidyr         * 1.3.1      2024-01-24 [1] CRAN (R 4.1.2)
# tidyselect      1.2.1      2024-03-11 [1] CRAN (R 4.1.2)
# tidyverse     * 2.0.0      2023-02-22 [1] CRAN (R 4.1.2)
# timechange      0.2.0      2023-01-11 [1] CRAN (R 4.1.2)
# TMB             1.9.11     2024-04-03 [1] CRAN (R 4.1.2)
# tweedie       * 2.3.5      2022-08-17 [1] CRAN (R 4.1.2)
# tzdb            0.3.0      2022-03-28 [1] CRAN (R 4.1.2)
# utf8            1.2.4      2023-10-22 [1] CRAN (R 4.1.2)
# uuid            1.1-0      2022-04-19 [1] CRAN (R 4.1.2)
# vctrs           0.6.5      2023-12-01 [1] CRAN (R 4.1.2)
# withr           3.0.0      2024-01-16 [1] CRAN (R 4.1.2)
# xtable          1.8-4      2019-04-21 [1] CRAN (R 4.1.2)
# zoo             1.8-11     2022-09-17 [1] CRAN (R 4.1.2)
# 
# [1] /home/fkch/R/x86_64-pc-linux-gnu-library/4.1
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
