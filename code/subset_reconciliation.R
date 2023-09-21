library(tidyverse)
library(magrittr)
library(future)
library(forecast)
library(doParallel)
library(foreach)
library(parallel)

# Setup
data_label <- commandArgs(trailingOnly = TRUE)
# data_label <- "simulation"
# data_label <- "tourism"
nlambda <- 20
MonARCH <- TRUE
workers <- parallelly::availableCores()
source("R/subset_reconcile.R")

# Utility function
reconcile_forecast <- function(index, fits, train, basefc, resids, test, S, nvalid,
                               method, method_name, nlambda,
                               deteriorate = FALSE, deteriorate_series, deteriorate_rate,
                               MIPFocus, Cuts, TimeLimit,
                               MIPVerbose, MonARCH, workers){
  
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
  
  Base <- list(y_tilde = base_forecasts, G = NA, z = NA, lambda_report = NA)
  BU <- subset.reconcile(base_forecasts = base_forecasts, S = S, method = "bu",
                         MonARCH = MonARCH, workers = workers)
  for(i in 1:length(method)){
    assign(method_name[i], 
           subset.reconcile(base_forecasts = base_forecasts, S = S, method = method[i], 
                            residuals = residuals, MonARCH = MonARCH, workers = workers))
    assign(paste0(method_name[i], "_subset"), 
           subset.reconcile(base_forecasts = base_forecasts, S = S,
                            method = method[i], residuals = residuals,
                            fitted_values = fitted_values, train_data = train_data, nvalid = nvalid,
                            subset = TRUE, ridge = TRUE, nlambda = nlambda,
                            MIPFocus = MIPFocus, Cuts = Cuts, TimeLimit = TimeLimit,
                            MIPVerbose = MIPVerbose, MonARCH = MonARCH, workers = workers))
    print(paste("===", method_name[i], "finished!"))
  }
  
  mget(c("Base", "BU", method_name, paste0(method_name, "_subset")))
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
if (data_label == "simulation"){
  nvalid = 16; MIPFocus = 0; Cuts = -1; TimeLimit = 600; MIPVerbose = FALSE
  method <- c("ols", "wls_struct", "wls_var", "mint_cov", "mint_shrink")
  method_name <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")
}

#----------------------------------------------------------------------
# Simulation data - correlation
## Total/Middle/Bottom: 3 levels, n = 7
## Training set:  1-100
## Test set:      1
#----------------------------------------------------------------------
if (grepl("corr", data_label)){
  nvalid = NULL; MIPFocus = 0; Cuts = -1; TimeLimit = 600; MIPVerbose = FALSE
  method <- c("ols", "wls_struct", "wls_var", "mint_cov", "mint_shrink")
  method_name <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")
}

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
if (data_label == "tourism"){
  nvalid = 12; MIPFocus = 3; Cuts = 2; TimeLimit = 600; MIPVerbose = FALSE
  method <- c("ols", "wls_struct", "wls_var", "mint_shrink")
  method_name <- c("OLS", "WLSs", "WLSv", "MinTs")
}

#----------------------------------------------------------------------
# ABS - Unemployed persons by Duration of job search, State and Territory
##
## 6291.0.55.001 - UM2 - Unemployed persons by Duration of job search, State and Territory, January 1991 onwards
## 
## Monthly series
## Duration of job search (Duration, 6) * State and territory (STT, 8): n = 63 series in total, nb = 48 series at the bottom level
##
## Training set:  2010Jan-2022Jul
## Test set:      2022Aug-2023Jul
#----------------------------------------------------------------------
if (data_label == "labour"){
  nvalid = 12; MIPFocus = 3; Cuts = 2; TimeLimit = 600; MIPVerbose = FALSE
  method <- c("ols", "wls_struct", "wls_var", "mint_shrink")
  method_name <- c("OLS", "WLSs", "WLSv", "MinTs")
}

#################################################
# Import base forecast results
#################################################
for (i in c("S", "fits", "resids", "train", "basefc","test")){
  assign(i, readRDS(file = paste0("data/", data_label, "_", i, ".rds")))
}
indices <- unique(fits$Index)

#################################################
# Reconcile forecasts
#################################################
reconsf <- indices |>
  purrr::map(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, S, nvalid,
                                         method, method_name, nlambda,
                                         deteriorate = FALSE, 
                                         MIPFocus = MIPFocus, Cuts = Cuts, TimeLimit = TimeLimit,
                                         MIPVerbose = MIPVerbose, MonARCH = MonARCH, workers = workers))
saveRDS(reconsf, file = paste0("data_new/", data_label, "_subset_reconsf.rds"))
rm(reconsf)
print("Reconciliation finished!")

#################################################
# Reconcile forecasts - deteriorate base forecasts for some time series
#################################################
if (data_label == "simulation"){
  scenario <- c("s1", "s2", "s3")
  deteriorate_series <- c("AA", "A", "Total")
  deteriorate_rate <- rep(1.5, 3)
  for (i in 1:length(scenario)){
    reconsf_s <- indices |>
      purrr::map(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, S, nvalid,
                                             method, method_name, nlambda,
                                             deteriorate = TRUE, 
                                             deteriorate_series = deteriorate_series[i],
                                             deteriorate_rate = deteriorate_rate[i], 
                                             MIPFocus = MIPFocus, Cuts = Cuts, TimeLimit = TimeLimit,
                                             MIPVerbose = MIPVerbose, MonARCH = MonARCH, workers = workers))
    saveRDS(reconsf_s, file = paste0("data_new/", data_label, "_subset_reconsf_", scenario[i], ".rds"))
    rm(reconsf_s)
    print(paste("Scenario", i, "finished!"))
  }
}


