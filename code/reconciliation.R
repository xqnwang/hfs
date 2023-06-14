library(tidyverse)
library(magrittr)
library(future)
library(forecast)

reticulate::use_python("/Users/xwan0362/Library/r-miniconda-arm64/bin/python3.10", required = T)
reticulate::source_python("Python/subset.py")
source("R/reconcile.R")

# Utility function
reconcile_forecast <- function(index, fits, train, basefc, resids, test, S,
                               method, method_name,
                               deteriorate = FALSE, deteriorate_series, deteriorate_rate,
                               MIPFocus, Cuts, TimeLimit,
                               MIPVerbose, SearchVerbose){
  
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
  BU <- reconcile(base_forecasts = base_forecasts, S = S, method = "bu")
  for(i in 1:length(method)){
    assign(method_name[i], 
           reconcile(base_forecasts = base_forecasts, S = S, method = method[i], 
                     residuals = residuals))
    assign(paste0(method_name[i], "_subset"), 
           reconcile(base_forecasts = base_forecasts, S = S,
                     method = method[i], residuals = residuals,
                     fitted_values = fitted_values, train_data = train_data,
                     subset = TRUE, ridge = TRUE,
                     MIPFocus = MIPFocus, Cuts = Cuts, TimeLimit = TimeLimit,
                     MIPVerbose = MIPVerbose, SearchVerbose = SearchVerbose))
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
data_label <- "simulation" # used to name the results to be imported and saved
MIPFocus = 0; Cuts = -1; TimeLimit = 600; MIPVerbose = FALSE; SearchVerbose = FALSE
method <- c("ols", "wls_struct", "wls_var", "mint_cov", "mint_shrink")
method_name <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")

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
MIPFocus = 3; Cuts = 2; TimeLimit = 600; MIPVerbose = FALSE; SearchVerbose = TRUE
method <- c("ols", "wls_struct", "wls_var", "mint_shrink")
method_name <- c("OLS", "WLSs", "WLSv", "MinTs")

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
                                         method, method_name,
                                         deteriorate = FALSE, 
                                         MIPFocus = MIPFocus, Cuts = Cuts, TimeLimit = TimeLimit,
                                         MIPVerbose = MIPVerbose, SearchVerbose = SearchVerbose),
          .progress = !SearchVerbose)
saveRDS(reconsf, file = paste0("data_new/", data_label, "_reconsf.rds"))
rm(reconsf)

#################################################
# Reconcile forecasts - deteriorate base forecasts for some time series
#################################################
scenario <- c("s1", "s2", "s3")
deteriorate_series <- c("AA", "A", "Total")
deteriorate_rate <- rep(1.5, 3)
for (i in 1:3){
  reconsf_s <- indices |>
    purrr::map(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, S,
                                           method, method_name,
                                           deteriorate = TRUE, 
                                           deteriorate_series = deteriorate_series[i],
                                           deteriorate_rate = deteriorate_rate[i], 
                                           MIPFocus = MIPFocus, Cuts = Cuts, TimeLimit = TimeLimit,
                                           MIPVerbose = MIPVerbose, SearchVerbose = SearchVerbose),
               .progress = !SearchVerbose)
  saveRDS(reconsf_s, file = paste0("data_new/", data_label, "_reconsf_", scenario[i], ".rds"))
  rm(reconsf_s)
  print(paste("i =", i, "finished!"))
}

