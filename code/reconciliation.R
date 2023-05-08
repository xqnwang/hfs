library(tidyverse)
library(magrittr)
library(future)
library(forecast)

source("R/mip_l0.R")
source("R/reconcile.R")

# Utility function
reconcile_forecast <- function(index, fits, train, basefc, resids, test, S,
                               deteriorate = FALSE, deteriorate_series, deteriorate_rate,
                               G_bench, nlambda_0 = 20, 
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
  
  BU <- reconcile(base_forecasts = base_forecasts, S = S, method = "bu")
  OLS <- reconcile(base_forecasts = base_forecasts, S = S, method = "ols")
  OLS_subset <- reconcile(base_forecasts = base_forecasts, S = S, method = "ols",
                          residuals = residuals, 
                          fitted_values = fitted_values, train_data = train_data,
                          subset = TRUE, G_bench = G_bench, nlambda_0 = nlambda_0,
                          parallel = parallel, workers = workers, .progress = .progress)
  WLSs <- reconcile(base_forecasts = base_forecasts, S = S, method = "wls_struct")
  WLSs_subset <- reconcile(base_forecasts = base_forecasts, S = S, method = "wls_struct",
                           residuals = residuals, 
                           fitted_values = fitted_values, train_data = train_data,
                           subset = TRUE, G_bench = G_bench, nlambda_0 = nlambda_0,
                           parallel = parallel, workers = workers, .progress = .progress)
  WLSv <- reconcile(base_forecasts = base_forecasts, S = S, method = "wls_var",
                    residuals = residuals)
  WLSv_subset <- reconcile(base_forecasts = base_forecasts, S = S, method = "wls_var",
                           residuals = residuals, 
                           fitted_values = fitted_values, train_data = train_data,
                           subset = TRUE, G_bench = G_bench, nlambda_0 = nlambda_0,
                           parallel = parallel, workers = workers, .progress = .progress)
  MinT <- reconcile(base_forecasts = base_forecasts, S = S, method = "mint_cov",
                    residuals = residuals)
  MinT_subset <- reconcile(base_forecasts = base_forecasts, S = S, method = "mint_cov",
                           residuals = residuals,
                           fitted_values = fitted_values, train_data = train_data,
                           subset = TRUE, G_bench = G_bench, nlambda_0 = nlambda_0,
                           parallel = parallel, workers = workers, .progress = .progress)
  MinTs <- reconcile(base_forecasts = base_forecasts, S = S, method = "mint_shrink",
                     residuals = residuals)
  MinTs_subset <- reconcile(base_forecasts = base_forecasts, S = S, method = "mint_shrink",
                            residuals = residuals, 
                            fitted_values = fitted_values, train_data = train_data,
                            subset = TRUE, G_bench = G_bench, nlambda_0 = nlambda_0,
                            parallel = parallel, workers = workers, .progress = .progress)
  
  list(Base = list(y_tilde = base_forecasts, G = NA, z = NA, lambda0_report = NA), 
       BU = BU, 
       OLS = OLS, OLS_subset = OLS_subset,
       WLSs = WLSs, WLSs_subset = WLSs_subset,
       WLSv = WLSv, WLSv_subset = WLSv_subset,
       MinT = MinT, MinT_subset = MinT_subset,
       MinTs = MinTs, MinTs_subset = MinTs_subset)
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
S <- rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
           diag(rep(1, 4))) # S matrix

#----------------------------------------------------------------------
# Australian domestic tourism (only considering hierarchical structure)
## Quarterly series from 1998Q1-2017Q4: 80 observations for each series
##
## Total/State: 2 levels, n = 9
## Training set:  1998Q1-2015Q4
## Test set:      2016Q1-2017Q4
#----------------------------------------------------------------------
data_label <- "tourism" # used to name the results to be imported and saved
workers <- 8 # number of workers used to run in parallel
S <- rbind(rep(1, 8), diag(rep(1, 8)))

#################################################
# Import base forecast results
#################################################
for (i in c("fits", "resids", "train", "basefc","test")){
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
plan(multisession, workers = workers)
reconsf <- indices |>
  map_fun(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, S,
                                      deteriorate = FALSE,
                                      G_bench = "Zero", nlambda_0 = 20,
                                      parallel = ifelse(length(indices) == 1, TRUE, FALSE),
                                      workers = workers, .progress = ifelse(length(indices) == 1, TRUE, FALSE)),
                    .progress = ifelse(length(indices) == 1, FALSE, TRUE))
saveRDS(reconsf, file = paste0("data/", data_label, "_reconsf.rds"))
rm(reconsf)

#################################################
# Reconcile forecasts - deteriorate base forecasts for some time series
#################################################
deteriorate_series <- c("AA", "B") # s1
# deteriorate_series <- c("Total") # s2
deteriorate_rate <- 1.5
scenario <- "s1"

plan(multisession, workers = workers)
reconsf_s <- indices |>
  map_fun(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, S,
                                      deteriorate = TRUE, 
                                      deteriorate_series = deteriorate_series,
                                      deteriorate_rate = deteriorate_rate,
                                      G_bench = "Zero", nlambda_0 = 20,
                                      parallel = ifelse(length(indices) == 1, TRUE, FALSE),
                                      workers = workers, .progress = ifelse(length(indices) == 1, TRUE, FALSE)),
          .progress = ifelse(length(indices) == 1, FALSE, TRUE))
saveRDS(reconsf_s, file = paste0("data/", data_label, "_reconsf_", scenario, ".rds"))


