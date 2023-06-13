library(tidyverse)
library(magrittr)
library(future)
library(forecast)

source("R/mip_l0.R")
source("R/reconcile.R")
source("R/mip_l0_diag.R") # large hierarchy
source("R/reconcile_diag.R") # large hierarchy
source("R/qp_l0.R") # QP with output z

# Utility function
reconcile_forecast <- function(index, fits, train, basefc, resids, test, S, large = FALSE,
                               deteriorate = FALSE, deteriorate_series, deteriorate_rate,
                               G_bench, nlambda = 20, unbiased = TRUE,
                               parallel = FALSE, workers = 8, .progress = FALSE){
  
  n <- NCOL(fits)
  fitted_values <- fits[fits$Index == index, -n] |> as.matrix()
  residuals <- resids[resids$Index == index, -n] |> as.matrix()
  train_data <- train[train$Index == index, -n] |> as.matrix()
  base_forecasts <- basefc[basefc$Index == index, -n] |> as.matrix()
  
  if (deteriorate){
    fitted_values[, deteriorate_series] <- fitted_values[, deteriorate_series]*deteriorate_rate
    base_forecasts[, deteriorate_series] <- base_forecasts[, deteriorate_series]*deteriorate_rate
    residuals <- train_data - fitted_values
  }
  
  if (large){
    f.reconcile <- reconcile_diag
  } else{
    f.reconcile <- reconcile
  }
  
  BU <- f.reconcile(base_forecasts = base_forecasts, S = S, method = "bu")
  OLS <- f.reconcile(base_forecasts = base_forecasts, S = S, method = "ols")
  OLS_subset <- f.reconcile(base_forecasts = base_forecasts, S = S,
                            method = "ols", residuals = residuals,
                            fitted_values = fitted_values, train_data = train_data,
                            subset = TRUE, G_bench = G_bench, nlambda = nlambda,
                            unbiased = unbiased, lambda_0 = c(0,0.1,1,10,100),
                            parallel = parallel, workers = workers,
                            .progress = .progress)
  print("OLS_subset Finished")
  # OLS_lasso <- f.reconcile(base_forecasts = base_forecasts, S = S, 
  #                           method = "ols", residuals = residuals, 
  #                           fitted_values = fitted_values, train_data = train_data,
  #                           lasso = TRUE, ELasso = FALSE, nlambda = nlambda,
  #                           parallel = parallel, workers = workers, 
  #                           .progress = .progress)
  # WLSs <- f.reconcile(base_forecasts = base_forecasts, S = S, method = "wls_struct")
  # WLSs_subset <- f.reconcile(base_forecasts = base_forecasts, S = S,
  #                            method = "wls_struct", residuals = residuals,
  #                            fitted_values = fitted_values, train_data = train_data,
  #                            subset = TRUE, G_bench = G_bench, nlambda = nlambda,
  #                            unbiased = unbiased,
  #                            parallel = parallel, workers = workers,
  #                            .progress = .progress)
  # print("WLSs_subset Finished")
  # WLSs_lasso <- f.reconcile(base_forecasts = base_forecasts, S = S, 
  #                          method = "wls_struct", residuals = residuals, 
  #                          fitted_values = fitted_values, train_data = train_data,
  #                          lasso = TRUE, ELasso = FALSE, nlambda = nlambda,
  #                          parallel = parallel, workers = workers, 
  #                          .progress = .progress)
  # WLSv <- f.reconcile(base_forecasts = base_forecasts, S = S, method = "wls_var",
  #                     residuals = residuals)
  # WLSv_subset <- f.reconcile(base_forecasts = base_forecasts, S = S,
  #                            method = "wls_var", residuals = residuals,
  #                            fitted_values = fitted_values, train_data = train_data,
  #                            subset = TRUE, G_bench = G_bench, nlambda = nlambda,
  #                            unbiased = unbiased,
  #                            parallel = parallel, workers = workers,
  #                            .progress = .progress)
  # print("WLSv_subset Finished")
  # WLSv_lasso <- f.reconcile(base_forecasts = base_forecasts, S = S, 
  #                           method = "wls_var", residuals = residuals, 
  #                           fitted_values = fitted_values, train_data = train_data,
  #                           lasso = TRUE, ELasso = FALSE, nlambda = nlambda,
  #                           parallel = parallel, workers = workers, 
  #                           .progress = .progress)
  # MinT <- f.reconcile(base_forecasts = base_forecasts, S = S, method = "mint_cov",
  #                     residuals = residuals)
  # MinT_subset <- f.reconcile(base_forecasts = base_forecasts, S = S,
  #                            method = "mint_cov", residuals = residuals,
  #                            fitted_values = fitted_values, train_data = train_data,
  #                            subset = TRUE, G_bench = G_bench, nlambda = nlambda,
  #                            unbiased = unbiased,
  #                            parallel = parallel, workers = workers,
  #                            .progress = .progress)
  # MinT_lasso <- f.reconcile(base_forecasts = base_forecasts, S = S, 
  #                           method = "mint_cov", residuals = residuals, 
  #                           fitted_values = fitted_values, train_data = train_data,
  #                           lasso = TRUE, ELasso = FALSE, nlambda = nlambda,
  #                           parallel = parallel, workers = workers, 
  #                           .progress = .progress)
  MinTs <- f.reconcile(base_forecasts = base_forecasts, S = S, method = "mint_shrink",
                       residuals = residuals)
  MinTs_subset <- f.reconcile(base_forecasts = base_forecasts, S = S,
                              method = "mint_shrink", residuals = residuals,
                              fitted_values = fitted_values, train_data = train_data,
                              subset = TRUE, G_bench = G_bench, nlambda = nlambda,
                              unbiased = unbiased, lambda_0 = c(0,0.1,1,10,100),
                              parallel = parallel, workers = workers,
                              .progress = .progress)
  print("MinTs_subset Finished")
  # MinTs_lasso <- f.reconcile(base_forecasts = base_forecasts, S = S, 
  #                           method = "mint_shrink", residuals = residuals, 
  #                           fitted_values = fitted_values, train_data = train_data,
  #                           lasso = TRUE, ELasso = FALSE, nlambda = nlambda,
  #                           parallel = parallel, workers = workers, 
  #                           .progress = .progress)
  # Emp_lasso <- f.reconcile(base_forecasts = base_forecasts, S = S, 
  #                            method = "ols", residuals = residuals, 
  #                            fitted_values = fitted_values, train_data = train_data,
  #                            lasso = TRUE, ELasso = TRUE, nlambda = nlambda,
  #                            parallel = parallel, workers = workers, 
  #                            .progress = .progress)
  
  list(Base = list(y_tilde = base_forecasts, G = NA, z = NA, lambda_report = NA), 
       BU = BU,
       OLS = OLS, OLS_subset = OLS_subset, 
       #OLS_lasso = OLS_lasso,
       #WLSs = WLSs, WLSs_subset = WLSs_subset, 
       #WLSs_lasso = WLSs_lasso,
       #WLSv = WLSv, WLSv_subset = WLSv_subset, 
       #WLSv_lasso = WLSv_lasso,
       #MinT = MinT, MinT_subset = MinT_subset, MinT_lasso = MinT_lasso,
       MinTs = MinTs, MinTs_subset = MinTs_subset
       #MinTs_lasso = MinTs_lasso,
       #Emp_lasso = Emp_lasso
       )
}

#################################################
# Data
#################################################
#----------------------------------------------------------------------
# Simulation data
## Total/Middle/Bottom: 3 levels, n = 7
## Training set:  1978Q1-2018Q4
## Test set:      2019Q1-2022Q4
#----------------------------------------------------------------------
data_label <- "simulation" # used to name the results to be imported and saved
workers <- 8 # number of workers used to run in parallel
large <- FALSE

#----------------------------------------------------------------------
# Australian domestic tourism (only considering hierarchical structure)
##
## Monthly series from 1998Jan-2017Dec: 240 months (20 years) for each series
##
## Total/State/Zone/Region: 4 levels, n = 111 series in total
##
## Training set:  1998Jan-2016Dec
## Test set:      2017Jan-2017Dec
#----------------------------------------------------------------------
data_label <- "tourism" # used to name the results to be imported and saved
workers <- 8 # number of workers used to run in parallel
large <- FALSE # use reconcile_diag
unbiased <- TRUE # only work when using reconcile_diag
# S <- rbind(rep(1, 7), diag(1, 7))
# S <- rbind(rep(1, 27),
#            c(rep(1, 6), rep(0, 21)),
#            c(rep(0, 6), rep(1, 5), rep(0, 16)),
#            c(rep(0, 11), rep(1, 4), rep(0, 12)),
#            c(rep(0, 15), rep(1, 4), rep(0, 8)),
#            c(rep(0, 19), rep(1, 3), rep(0, 5)),
#            c(rep(0, 22), rep(1, 3), rep(0, 2)),
#            c(rep(0, 25), rep(1, 2)),
#            diag(1, 27))

#################################################
# Import base forecast results
#################################################
for (i in c("S", "fits", "resids", "train", "basefc","test")){
  assign(i, readRDS(file = paste0("data/", data_label, "_", i, ".rds")))
}

indices <- unique(fits$Index)
if (length(indices) == 1){
  map_fun <- purrr::map
} else {
  future::plan(multisession, workers = workers)
  map_fun <- furrr::future_map
}

#################################################
# Reconcile forecasts
#################################################
reconsf <- indices |>
  map_fun(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, S, 
                                      large = large, deteriorate = FALSE,
                                      G_bench = "Zero", nlambda = 5,
                                      parallel = FALSE,
                                      unbiased = unbiased,
                                      workers = workers, .progress = TRUE),
                    .progress = ifelse(length(indices) == 1, FALSE, TRUE))
saveRDS(reconsf, file = paste0("data/", data_label, "_reconsf_lasso.rds"))
rm(reconsf)

#################################################
# Reconcile forecasts - deteriorate base forecasts for some time series
#################################################
deteriorate_series <- c("AA") # s1  deteriorate_series <- c("WA", "QLD")
# deteriorate_series <- c("Total")
deteriorate_rate <- 1.5
scenario <- "s1"
# scenario <- "s2"

plan(multisession, workers = workers)
reconsf_s <- indices |>
  map_fun(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, S,
                                      large = large, deteriorate = TRUE, 
                                      deteriorate_series = deteriorate_series,
                                      deteriorate_rate = deteriorate_rate,
                                      G_bench = "Zero", nlambda = 20,
                                      parallel = ifelse(length(indices) == 1, TRUE, FALSE),
                                      workers = workers, .progress = ifelse(length(indices) == 1, TRUE, FALSE)),
          .progress = ifelse(length(indices) == 1, FALSE, TRUE))
saveRDS(reconsf_s, file = paste0("data/", data_label, "_reconsf_lasso_", scenario, ".rds"))


