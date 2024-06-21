#' ---
#' title: Application of adjusted predictions from a GEE to the spruces dataset
#' author: Francis KC Hui
#' date: Code started Feb 2024
#' ---

rm(list = ls())
library(tidyverse)
library(colorspace)
library(patchwork)
library(glmtoolbox)
library(geepack)
library(ggforce)
library(glmmTMB)
library(splines)

here::i_am("application_spruces/application.R")
library(here)

source(here("R", "ICs.R"))
source(here("R", "utils.R"))
source(here("R", "predictions.R"))



##-------------------------
#' # Load in and wrangle/explore data
#' Note from Example 1 of [https://journal.r-project.org/articles/RJ-2023-056/] it was shown that a Gamma distribution with a quadratic mean-variance function and an AR marginal correlation structure was the best selected GEE 
##-------------------------
data(spruces)
?spruces

ggplot(spruces, aes(x = days, y = size, group = tree, color = treat)) +
    geom_line(alpha = 0.7) +
    theme_bw()


table(spruces$tree)
table(spruces$days)
head(spruces)

spruces <- spruces %>% 
    rename(id = "tree") 
#' Note data are balanced: 
#' In the first group, 54 trees were grown under an ozone-enriched atmosphere, ozone exposure at 70 parts per billion. 
#' In the second group, 25 trees were grown under normal conditions. 
#' The size of each tree was observed 13 times across time, that is, 152, 174, 201, 227, 258, 469, 496, 528, 556, 579, 613, 639 and 674 day



##-------------------------
#' # Fit GEE to entire dataset and construct predictions for future time points
##-------------------------
fit_gees <- geepack::geeglm(formula = size ~ ns(days, df = 4) + treat,
                            id = id,
                            data = spruces,
                            family = Gamma(link = "log"),
                            corstr = "ar1",
                            scale.fix = FALSE)
summary(fit_gees)


#' ## Add ten extra days beyond the final day in the data
expanded_spruces <- spruces %>% 
    group_by(id) %>% 
    complete(days = seq(700, by = 25, length = 10)) %>% 
    mutate(treat = treat %>% na.omit %>% unique) %>% 
    ungroup %>% 
    arrange(id, days)
expanded_spruces <- expanded_spruces %>% 
    arrange(id) %>% 
    bind_cols(., 
              gee_predictions(fit_gees, newdata = expanded_spruces, intervals = TRUE) %>% select(-id)) %>% 
    mutate(current = rep(c(TRUE,FALSE), c(13,10)) %>% rep(., expanded_spruces$id %>% unique %>% length))

head(expanded_spruces)


ggplot(expanded_spruces %>% 
           pivot_longer(c(standard, adjusted), names_to = "type") %>% 
           mutate(type = fct_inorder(type)) %>% 
           mutate(type = fct_recode(type, 
                                    "Standard GEE predictions" = "standard", 
                                    "Adjusted GEE predictions" = "adjusted")), 
       aes(x = days, y = size, group = id, color = treat)) +
    geom_line(alpha = 0.25) +
    geom_point(aes(x = days, y = value, color = treat, group = id), alpha = 0.8) +
    geom_line(aes(x = days, y = value, color = treat, group = id), alpha = 0.8) +
    facet_wrap(. ~ type, nrow = 1) + 
    scale_color_discrete_diverging(palette = "Blue-Yellow 3") +
    #scale_color_viridis_d() +
    scale_y_log10() +
    labs(x = "Days", y = "Size", color = "Treatment") +
    theme_bw() +
    theme(legend.position = "bottom")


ggsave(file = here("application_spruces", "spruce_predictions.pdf"), width = 8, height = 5)



ggplot(expanded_spruces %>% 
           select(id:treat, standard, standard2_lower, standard2_upper),
       aes(x = days, y = size, group = id, color = treat, ymin = standard2_lower, ymax = standard2_upper)) +
    geom_line(alpha = 0.25) +
    geom_ribbon(aes(fill = treat),linetype = 3, show.legend = FALSE) +
    geom_point(aes(x = days, y = standard, color = treat, group = id), alpha = 0.8) +
    geom_line(aes(x = days, y = standard, color = treat, group = id), alpha = 0.8) +
    scale_color_discrete_diverging(palette = "Blue-Yellow 3") +
    scale_fill_manual(values = c("cadetblue1", "lightyellow")) +
    labs(x = "Days", y = "Size", color = "Treatment", fill = "Treatment", title = "Standard GEE predictions") +
    theme_bw() +
    theme(legend.position = "bottom")


ggsave(file = here("application_spruces","spruce_standardpredictions.pdf"), height = 8, width = 8)




p <- ggplot(expanded_spruces %>% 
                select(id:treat, adjusted, adjusted_lower, adjusted_upper),
            aes(x = days, y = size, color = treat, ymin = adjusted_lower, ymax = adjusted_upper)) +
    geom_line(alpha = 0.25) + 
    geom_ribbon(aes(fill = treat), linetype = 3, show.legend = FALSE) +
    geom_point(aes(x = days, y = adjusted, color = treat), alpha = 0.8) +
    geom_line(aes(x = days, y = adjusted, color = treat), alpha = 0.8) +
    facet_wrap_paginate(. ~ id, ncol = 2, nrow = 2, page = 20) +
    scale_color_discrete_diverging(palette = "Blue-Yellow 3") +
    scale_fill_manual(values = c("cadetblue1", "lightyellow")) +
    labs(x = "Days", y = "Size", color = "Treatment", fill = "Treatment", title = "Adjusted GEE predictions") +
    theme_bw() +
    theme(legend.position = "bottom") 


pdf(file = here("application_spruces","spruces_adjustedpredictions.pdf"), width = 8, height = 8)
for(k0 in 1:n_pages(p)) {
    p <- ggplot(expanded_spruces %>% 
                    select(id:treat, adjusted, adjusted_lower, adjusted_upper),
                aes(x = days, y = size, color = treat, ymin = adjusted_lower, ymax = adjusted_upper)) +
        geom_line(alpha = 0.25) + 
        geom_ribbon(aes(fill = treat), linetype = 3, show.legend = FALSE) +
        geom_point(aes(x = days, y = adjusted, color = treat), alpha = 0.8) +
        geom_line(aes(x = days, y = adjusted, color = treat), alpha = 0.8) +
        facet_wrap_paginate(. ~ id, ncol = 2, nrow = 2, page = k0) +
        scale_color_discrete_diverging(palette = "Blue-Yellow 3") +
        scale_fill_manual(values = c("cadetblue1", "lightyellow")) +
        labs(x = "Days", y = "Size", color = "Treatment", fill = "Treatment", title = "Adjusted GEE predictions") +
        theme_bw() +
        theme(legend.position = "bottom")
    print(p)
}
dev.off()


##-------------------------
#' # Perform a rolling origin CV, starting from five time points for each tree
##-------------------------
rollingCV <- function(starting_time_point = 5,
                      testing_period = 3,
                      response_type = Gamma(link = "log"),
                      GEE_formula = size ~ ns(days, df = 4) + treat) {

    
    z_alpha <- qnorm(0.975)
    training_time_points <- starting_time_point:(min(table(spruces$id)) - testing_period)
        
    inner_fn <- function(k0) {
        message("Current fold/data comprises the first ", training_time_points[k0], " time points...")
        
        #' ## Split data into training and test
        train_dat <- spruces %>%
            group_by(id) %>% 
            filter(days %in% days[1:training_time_points[k0]]) %>%
            ungroup
        test_dat <- spruces %>% 
            group_by(id) %>% 
            filter(days %in% days[training_time_points[k0] + (1:testing_period)]) %>% 
            ungroup

    
        #' ## Fit GEEs to training , selecting the marginal correlation structure via GBIC (SGPC), and predict to test
        correlation_structure_choices <- c("Independence", "AR-M-dependent(1)", "Exchangeable") 
        fit_gees <- lapply(correlation_structure_choices, function(x) {
            glmtoolbox::glmgee(formula = GEE_formula,
                               id = id, 
                               data = train_dat,
                               family = response_type, 
                               corstr = x, 
                               scale.fix = FALSE) 
            })
        pick_best <- bind_cols(QIC = glmtoolbox::QIC(fit_gees[[1]], fit_gees[[2]], fit_gees[[3]], verbose = FALSE)$QIC,
                               CIC = glmtoolbox::CIC(fit_gees[[1]], fit_gees[[2]], fit_gees[[3]], verbose = FALSE)$CIC,
                               PAC = glmtoolbox::PAC(fit_gees[[1]], fit_gees[[2]], fit_gees[[3]], verbose = FALSE)$PAC,
                               getmoreICs(fit_gees[[1]], fit_gees[[2]], fit_gees[[3]])
                               )
        
        correlation_structure_choices <- c("independence", "ar1", "exchangeable")
        fit_gees <- geepack::geeglm(formula = GEE_formula,
                                    id = id,
                                    data = train_dat,
                                    family = response_type,
                                    corstr = correlation_structure_choices[which.min(pick_best[,7])],
                                    scale.fix = FALSE)
        
        best_predictions_SGPC <- gee_predictions(object = fit_gees, newdata = test_dat %>% select(-size), intervals = TRUE)
        
        
        #' ## Fit GLMM to training data, assuming either random intercept or random intercept plus slope for time. Then construct conditional predictions to test data
        fit_glmm <- glmmTMB(update.formula(GEE_formula, ~ . + (1|id)), 
                            data = train_dat, 
                            family = response_type)   
        y_test_glmm <- predict(fit_glmm, newdata = test_dat, type = "response", se.fit = TRUE) 
        y_test_glmm$lower_limit <- y_test_glmm$fit - z_alpha * y_test_glmm$se.fit
        y_test_glmm$upper_limit <- y_test_glmm$fit + z_alpha * y_test_glmm$se.fit
        
        fit_glmm2 <- glmmTMB(update.formula(GEE_formula, ~ . + (1 + days|id)), 
                            data = train_dat, 
                            family = response_type)   
        y_test_glmm2 <- predict(fit_glmm2, newdata = test_dat, type = "response", se.fit = TRUE) 
        y_test_glmm2$lower_limit <- y_test_glmm2$fit - z_alpha * y_test_glmm2$se.fit
        y_test_glmm2$upper_limit <- y_test_glmm2$fit + z_alpha * y_test_glmm2$se.fit
        
        
        #' ## Wrangle output
        #' Aggregated across all clusters
        out <- data.frame(
            timepoints = test_dat$days,
            standardpred_SGPC = best_predictions_SGPC$standard, 
            adjustedpred_SGPC = best_predictions_SGPC$adjusted, 
            glmm_randomintercept = y_test_glmm$fit, 
            glmm_randomintercepttime = y_test_glmm2$fit, 
            true = test_dat$size) %>% 
            pivot_longer(-c(true, timepoints), values_to = "predictions", names_to = "method") %>% 
            mutate(method = fct_inorder(method)) %>% 
            group_by(method) %>%
            summarise(RMSE = sqrt(mean((predictions - true)^2, na.rm = TRUE)),
                      #MAE = mean(abs(predictions - true), na.rm = TRUE),
                      pearson_cor = cor(true, predictions),
                      spearman_cor = cor(true, predictions, method = "spearman")
            ) %>% 
            arrange(method) %>% 
            as.data.frame %>% 
            mutate(num_currenttimepoints = training_time_points[k0])
        
        
        #' Uncertainty interval performance
        interval_score_fn <- function(lower_limit, upper_limit, true) {
            out <- upper_limit - lower_limit + 2/0.05*(lower_limit - true)*as.numeric(true < lower_limit) + 2/0.05*(true - upper_limit)*as.numeric(true > upper_limit)
            return(out)
            }
        
        out_UI <- data.frame(
            timepoints = test_dat$days,
            standardpred_coverage = (best_predictions_SGPC$standard_lower < test_dat$size & best_predictions_SGPC$standard_upper > test_dat$size),
            standardpred_width = (best_predictions_SGPC$standard_upper - best_predictions_SGPC$standard_lower),
            standardpred_intervalscore = interval_score_fn(lower_limit = best_predictions_SGPC$standard_lower, 
                                                           upper_limit = best_predictions_SGPC$standard_upper, 
                                                           true = test_dat$size),
            standard2pred_coverage = (best_predictions_SGPC$standard2_lower < test_dat$size & best_predictions_SGPC$standard2_upper > test_dat$size),
            standard2pred_width = (best_predictions_SGPC$standard2_upper - best_predictions_SGPC$standard2_lower),
            standard2pred_intervalscore = interval_score_fn(lower_limit = best_predictions_SGPC$standard2_lower, 
                                                           upper_limit = best_predictions_SGPC$standard2_upper, 
                                                           true = test_dat$size),
            adjustedpred_coverage = (best_predictions_SGPC$adjusted_lower < test_dat$size & best_predictions_SGPC$adjusted_upper > test_dat$size),
            adjustedpred_width = (best_predictions_SGPC$adjusted_upper - best_predictions_SGPC$adjusted_lower),
            adjustedpred_intervalscore = interval_score_fn(lower_limit = best_predictions_SGPC$adjusted_lower, 
                                                           upper_limit = best_predictions_SGPC$adjusted_upper, 
                                                           true = test_dat$size),
            glmm_randomintercept_coverage = (y_test_glmm$lower_lim < test_dat$size & y_test_glmm$upper_lim > test_dat$size),
            glmm_randomintercept_width = (y_test_glmm$upper_lim - y_test_glmm$lower_lim),
            glmm_randomintercept_intervalscore = interval_score_fn(lower_limit = y_test_glmm$lower_limit, 
                                                                   upper_limit = y_test_glmm$upper_limit, 
                                                                   true = test_dat$size),
            glmm_randomintercepttime_coverage = (y_test_glmm2$lower_lim < test_dat$size & y_test_glmm2$upper_lim > test_dat$size),
            glmm_randomintercepttime_width = (y_test_glmm2$upper_lim - y_test_glmm2$lower_lim),
            glmm_randomintercepttime_intervalscore = interval_score_fn(lower_limit = y_test_glmm2$lower_limit, 
                                                                       upper_limit = y_test_glmm2$upper_limit, 
                                                                       true = test_dat$size),
            true = test_dat$size) %>% 
            mutate(num_currenttimepoints = training_time_points[k0]) %>% 
            colMeans(., na.rm = TRUE)


        return(list(out = out, 
                    out_UI = out_UI, 
                    selected_correlation_structure = apply(pick_best, 2, which.min)))
        }
    

    out <- lapply(1:length(training_time_points), inner_fn)
    return(out)
    }



allresults <- rollingCV()


sapply(allresults, function(x) x$selected_correlation_structure) 
#' Note GBIC (SGPC) always choose the AR1 structure; this is not surprising given the results of [https://journal.r-project.org/articles/RJ-2023-056/]


lapply(allresults, function(x) x$out) %>% 
    do.call(rbind, .) %>% 
    as.data.frame %>% 
    select(-c(deviance, MAE)) %>% 
    group_by(method) %>% 
    summarise(across(RMSE:spearman_cor, ~ mean(.x)))


s <- lapply(allresults, function(x) x$out) %>% 
    do.call(rbind, .) %>% 
    as.data.frame %>% 
    pivot_longer(RMSE:spearman_cor, names_to = "metric") %>% 
    mutate(method = fct_inorder(method)) %>% 
    mutate(method = fct_recode(method,
                               "GEE (standard)" = "standardpred_SGPC", 
                               "GEE (adjusted)" = "adjustedpred_SGPC",
                               "GLMM (random int.)" = "glmm_randomintercept",
                               "GLMM (random int. and slope)" = "glmm_randomintercepttime")) %>% 
    mutate(metric = fct_inorder(metric)) %>% 
    mutate(metric = fct_recode(metric, "Pearson correlation" = "pearson_cor", "Spearman correlation" = "spearman_cor"))

s2 <- lapply(allresults, function(x) x$out_UI) %>% 
    do.call(rbind, .) %>% 
    as.data.frame() %>% 
    select(c(num_currenttimepoints, ends_with("score"))) %>% 
    pivot_longer(-num_currenttimepoints, names_to = "method") %>% 
    mutate(method = fct_inorder(method)) %>% 
    mutate(method = fct_recode(method,
                               "GEE (standard) version 2" = "standardpred_intervalscore", # These two are swapped around as the standard2pred_intervalscore does better and so we present that belo
                               "GEE (standard)" = "standard2pred_intervalscore", 
                               "GEE (adjusted)" = "adjustedpred_intervalscore",
                               "GLMM (random int.)" = "glmm_randomintercept_intervalscore",
                               "GLMM (random int. and slope)" = "glmm_randomintercepttime_intervalscore")) %>% 
    mutate(metric = "Interval Score") %>% 
    relocate(method, num_currenttimepoints, metric, value)

s <- bind_rows(s, s2) %>% 
    arrange(num_currenttimepoints, method) %>% 
    mutate(metric = fct_inorder(metric)) %>% 
    filter(method != "GEE (standard) version 2")



#' Plotting how metrics computed across all clusters varies with number of time points in current dataset
p1 <- ggplot(s, aes(x = num_currenttimepoints, y = value, color = method, shape = method)) +
    geom_point(size = 3) +
    geom_line() +
    facet_wrap(. ~ metric, nrow = 2, scales = "free") +
    scale_color_discrete_qualitative() +
    scale_y_log10() +
    labs(x = "No. of time points in current dataset", y = "Metric", color = "Method", shape = "Method") +
    theme_bw() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 1,byrow = TRUE))
p1


ggsave(p1, file = here("application_spruces", "rollingCVresults.pdf"), width = 10, height = 8)




##----------------------
sessioninfo::session_info()
##----------------------
# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 (2021-11-01)
# os       Linux Mint 21.1
# system   x86_64, linux-gnu
# ui       RStudio
# language en_AU:en
# collate  en_AU.UTF-8
# ctype    en_AU.UTF-8
# tz       Australia/Sydney
# date     2024-06-21
# rstudio  2024.04.0+735 Chocolate Cosmos (desktop)
# pandoc   NA
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package      * version    date (UTC) lib source
# backports      1.5.0      2024-05-23 [1] CRAN (R 4.1.2)
# boot           1.3-28     2021-05-03 [4] CRAN (R 4.1.1)
# broom          1.0.5      2023-06-09 [1] CRAN (R 4.1.2)
# cli            3.6.2      2023-12-11 [1] CRAN (R 4.1.2)
# coda           0.19-4.1   2024-01-31 [1] CRAN (R 4.1.2)
# codetools      0.2-18     2020-11-04 [4] CRAN (R 4.0.3)
# colorspace   * 2.1-0      2023-01-23 [1] CRAN (R 4.1.2)
# dplyr        * 1.1.4      2023-11-17 [1] CRAN (R 4.1.2)
# emmeans        1.9.0      2023-12-18 [1] CRAN (R 4.1.2)
# estimability   1.4.1      2022-08-05 [1] CRAN (R 4.1.2)
# fansi          1.0.6      2023-12-08 [1] CRAN (R 4.1.2)
# farver         2.1.2      2024-05-13 [1] CRAN (R 4.1.2)
# forcats      * 1.0.0      2023-01-29 [1] CRAN (R 4.1.2)
# Formula        1.2-5      2023-02-24 [1] CRAN (R 4.1.2)
# geepack      * 1.3.9      2022-08-16 [1] CRAN (R 4.1.2)
# generics       0.1.3      2022-07-05 [1] CRAN (R 4.1.2)
# ggforce      * 0.4.1      2022-10-04 [1] CRAN (R 4.1.2)
# ggplot2      * 3.5.1      2024-04-23 [1] CRAN (R 4.1.2)
# glmmTMB      * 1.1.9      2024-03-20 [1] CRAN (R 4.1.2)
# glmtoolbox   * 0.1.9      2023-10-10 [1] CRAN (R 4.1.2)
# glue           1.7.0      2024-01-09 [1] CRAN (R 4.1.2)
# gtable         0.3.5      2024-04-22 [1] CRAN (R 4.1.2)
# here         * 1.0.1      2020-12-13 [1] CRAN (R 4.1.2)
# hms            1.1.3      2023-03-21 [1] CRAN (R 4.1.2)
# labeling       0.4.3      2023-08-29 [1] CRAN (R 4.1.2)
# lattice        0.20-45    2021-09-22 [4] CRAN (R 4.1.1)
# lifecycle      1.0.4      2023-11-07 [1] CRAN (R 4.1.2)
# lme4           1.1-35.3   2024-04-16 [1] CRAN (R 4.1.2)
# lubridate    * 1.9.2      2023-02-10 [1] CRAN (R 4.1.2)
# magrittr       2.0.3      2022-03-30 [1] CRAN (R 4.1.2)
# MASS           7.3-55     2022-01-13 [4] CRAN (R 4.1.2)
# Matrix       * 1.6-4      2023-11-30 [1] CRAN (R 4.1.2)
# mgcv           1.9-0      2023-07-11 [1] CRAN (R 4.1.2)
# minqa          1.2.5      2022-10-19 [1] CRAN (R 4.1.2)
# multcomp       1.4-23     2023-03-09 [1] CRAN (R 4.1.2)
# munsell        0.5.1      2024-04-01 [1] CRAN (R 4.1.2)
# mvtnorm      * 1.2-5      2024-05-21 [1] CRAN (R 4.1.2)
# nlme           3.1-155    2022-01-13 [4] CRAN (R 4.1.2)
# nloptr         2.0.3      2022-05-26 [1] CRAN (R 4.1.2)
# numDeriv       2016.8-1.1 2019-06-06 [1] CRAN (R 4.1.2)
# patchwork    * 1.2.0      2024-01-08 [1] CRAN (R 4.1.2)
# pillar         1.9.0      2023-03-22 [1] CRAN (R 4.1.2)
# pkgconfig      2.0.3      2019-09-22 [1] CRAN (R 4.1.2)
# pkgload        1.3.2      2022-11-16 [1] CRAN (R 4.1.2)
# polyclip       1.10-4     2022-10-20 [1] CRAN (R 4.1.2)
# purrr        * 1.0.2      2023-08-10 [1] CRAN (R 4.1.2)
# R6             2.5.1      2021-08-19 [1] CRAN (R 4.1.2)
# ragg           1.2.5      2023-01-12 [1] CRAN (R 4.1.2)
# Rcpp           1.0.12     2024-01-09 [1] CRAN (R 4.1.2)
# RcppZiggurat   0.1.6      2020-10-20 [1] CRAN (R 4.1.2)
# readr        * 2.1.4      2023-02-10 [1] CRAN (R 4.1.2)
# Rfast          2.0.8      2023-07-03 [1] CRAN (R 4.1.2)
# rlang          1.1.3      2024-01-10 [1] CRAN (R 4.1.2)
# rprojroot      2.0.4      2023-11-05 [1] CRAN (R 4.1.2)
# rstudioapi     0.16.0     2024-03-24 [1] CRAN (R 4.1.2)
# sandwich       3.0-2      2022-06-15 [1] CRAN (R 4.1.2)
# scales         1.3.0      2023-11-28 [1] CRAN (R 4.1.2)
# sessioninfo    1.2.2      2021-12-06 [1] CRAN (R 4.1.2)
# statmod        1.5.0      2023-01-06 [1] CRAN (R 4.1.2)
# stringi        1.8.4      2024-05-06 [1] CRAN (R 4.1.2)
# stringr      * 1.5.1      2023-11-14 [1] CRAN (R 4.1.2)
# survival       3.2-13     2021-08-24 [4] CRAN (R 4.1.1)
# systemfonts    1.0.4      2022-02-11 [1] CRAN (R 4.1.2)
# textshaping    0.3.6      2021-10-13 [1] CRAN (R 4.1.2)
# TH.data        1.1-2      2023-04-17 [1] CRAN (R 4.1.2)
# tibble       * 3.2.1      2023-03-20 [1] CRAN (R 4.1.2)
# tidyr        * 1.3.1      2024-01-24 [1] CRAN (R 4.1.2)
# tidyselect     1.2.1      2024-03-11 [1] CRAN (R 4.1.2)
# tidyverse    * 2.0.0      2023-02-22 [1] CRAN (R 4.1.2)
# timechange     0.2.0      2023-01-11 [1] CRAN (R 4.1.2)
# TMB            1.9.11     2024-04-03 [1] CRAN (R 4.1.2)
# tweenr         2.0.2      2022-09-06 [1] CRAN (R 4.1.2)
# tzdb           0.3.0      2022-03-28 [1] CRAN (R 4.1.2)
# utf8           1.2.4      2023-10-22 [1] CRAN (R 4.1.2)
# vctrs          0.6.5      2023-12-01 [1] CRAN (R 4.1.2)
# withr          3.0.0      2024-01-16 [1] CRAN (R 4.1.2)
# xtable         1.8-4      2019-04-21 [1] CRAN (R 4.1.2)
# zoo            1.8-11     2022-09-17 [1] CRAN (R 4.1.2)
# 
# [1] /home/fkch/R/x86_64-pc-linux-gnu-library/4.1
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
# 
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────