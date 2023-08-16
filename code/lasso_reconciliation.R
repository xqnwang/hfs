library(tidyverse)
library(magrittr)
library(future)
library(forecast)

# Setup
data_label <- "simulation"
# data_label <- "tourism"
nlambda <- 20
MonARCH <- TRUE
workers <- parallel::detectCores()
if (MonARCH){
  path <- "~/.local/share/r-miniconda/envs/r-reticulate/bin/python3.8"
  setwd(Sys.glob(file.path("~/wm15/", "*", "hfs")))
} else{
  path <- "~/Library/r-miniconda-arm64/bin/python3.10"
}
reticulate::use_python(path, required = T)
reticulate::source_python("Python/glasso.py")
reticulate::source_python("Python/eglasso.py")
source("R/lasso_reconcile.R")

# Utility function
reconcile_forecast <- function(index, fits, train, basefc, resids, test, S,
                               method, method_name, nlambda,
                               deteriorate = FALSE, deteriorate_series, deteriorate_rate,
                               SearchVerbose){
  
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
  BU <- lasso.reconcile(base_forecasts = base_forecasts, S = S, method = "bu", lasso = NULL)
  for(i in 1:length(method)){
    assign(method_name[i],
           lasso.reconcile(base_forecasts = base_forecasts, S = S, method = method[i],
                           residuals = residuals, lasso = NULL))
    assign(paste0(method_name[i], "_Lasso"),
           lasso.reconcile(base_forecasts = base_forecasts, S = S,
                           method = method[i], residuals = residuals,
                           fitted_values = fitted_values, train_data = train_data,
                           lasso = "Lasso", nlambda = nlambda,
                           SearchVerbose = SearchVerbose))
  }
  ELasso <- lasso.reconcile(base_forecasts = base_forecasts, S = S,
                            method = "ols", residuals = residuals,
                            fitted_values = fitted_values, train_data = train_data,
                            lasso = "ELasso", nlambda = nlambda,
                            SearchVerbose = SearchVerbose)
  mget(c("Base", "BU", method_name, paste0(method_name, "_Lasso"), "ELasso"))
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
  SearchVerbose = FALSE
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
  SearchVerbose = TRUE
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
  purrr::map(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, S,
                                         method, method_name, nlambda,
                                         SearchVerbose = SearchVerbose),
          .progress = !SearchVerbose)
saveRDS(reconsf, file = paste0("data_new/", data_label, "_lasso_reconsf.rds"))

#################################################
# Reconcile forecasts - deteriorate base forecasts for some time series
#################################################
if (data_label == "simulation"){
  scenario <- c("s1", "s2", "s3")
  deteriorate_series <- c("AA", "A", "Total")
  deteriorate_rate <- rep(1.5, 3)
  for (i in 1:3){
    reconsf_s <- indices |>
      purrr::map(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, S,
                                             method, method_name, nlambda,
                                             deteriorate = TRUE, 
                                             deteriorate_series = deteriorate_series[i],
                                             deteriorate_rate = deteriorate_rate[i], 
                                             SearchVerbose),
                 .progress = !SearchVerbose
      )
    saveRDS(reconsf_s, file = paste0("data_new/", data_label, "_lasso_reconsf_", scenario[i], ".rds"))
    rm(reconsf_s)
    print(paste0("Scenario s", i, " finished!"))
  } 
}


