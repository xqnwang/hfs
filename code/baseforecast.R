library(tidyverse)
library(magrittr)
library(future)
library(forecast)

data_label <- commandArgs(trailingOnly = TRUE) # "simulation", "corr_i", "tourism_i", or "labour_i"

# Utility function
base_forecast <- function(hts, method, h, special = NULL){
  out <- list()
  for(i in 1:NCOL(hts)){
    series <- hts[, i]
    seriesname <- colnames(hts)[i]
    if (method == 'arima'){
      model <- auto.arima(series)
    }
    if (method == 'ets'){
      model <- ets(series)
    }
    if (method == "ar"){
      if (is.null(special)){
        model <- auto.arima(series, allowmean = TRUE, stationary = TRUE, stepwise = FALSE)
      } else{
        model <- auto.arima(series, allowmean = !(seriesname %in% special), stationary = TRUE, stepwise = FALSE)
      }
    }
    fc <- forecast(model, h)$mean
    resid <- residuals(model, type = 'response') # use regular residuals rather than innovation residuals
    out[[i]] <- list(fitted = model$fitted, residuals = resid, 
                     train = model$x, forecast = fc)
  }
  names(out) <- colnames(hts)
  out
}

#################################################
# Import data
#################################################
#----------------------------------------------------------------------
# Simulation data
## Total/Middle/Bottom: 3 levels, n = 7
## Training set:  1978Q1-2018Q4
## Test set:      2019Q1-2022Q4
#----------------------------------------------------------------------
if (data_label == "simulation"){
  # Formalize data (Index, Time, Series1, ..., Seriesn)
  dat <- readr::read_csv(paste0("data/", data_label, "_data.csv")) |>
    mutate(Time = tsibble::yearquarter(Quarter),
           A = AA + AB,
           B = BA + BB,
           Total = AA + AB + BA + BB) |>
    select(c(Index, Time, Total, A, B, AA, AB, BA, BB)) 
  
  # Data details
  freq <- 4
  start_train <- c(1978, 1)
  end_train <- c(2018, 4)
  start_test <- c(2019, 1)
  end_test <- c(2022, 4)
  
  # Forecasting method
  fmethod <- "ets"
  special <- NULL
  
  # S matrix
  S <- rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
             diag(rep(1, 4)))
}

#----------------------------------------------------------------------
# Simulation data - correlation
## Total/Middle/Bottom: 3 levels, n = 7
## Training set:  1-100
## Test set:      1
#----------------------------------------------------------------------
if (grepl("corr", data_label)){
  freq <- 1
  start_train <- 1
  end_train <- 100
  start_test <- 101
  end_test <- 101
  dat <- readr::read_csv(paste0("data/", data_label, "_data.csv")) |>
    mutate(A = AA + AB,
           B = BA + BB,
           Total = AA + AB + BA + BB) |>
    select(c(Index, Time, Total, A, B, AA, AB, BA, BB))
  
  # Forecasting method
  fmethod <- "ar"
  special <- c("Total", "A", "BA")
  
  # S matrix
  S <- rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
             diag(rep(1, 4)))
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
if (grepl("tourism", data_label)){
  # Formalize data (Index, Time, Series1, ..., Seriesn)
  dat <- readRDS("data/tourism_data.rds")
  h <- 12
  
  # Data details
  k <- sub('.*_', '', data_label) |> as.numeric() # 1:12
  freq <- 12
  start_train <- c(1998, 1)
  data_eg <- dat |>
    filter(Index == 1) |>
    select(!c(Index, Time)) |>
    ts(frequency = freq, start = start_train)
  
  length_train <- 240 - h - k + 1
  end_train <- time(data_eg)[length_train]
  start_test <- time(data_eg)[length_train + 1]
  end_test <- time(data_eg)[length_train + h]
  
  # Forecasting method
  fmethod <- "ets"
  special <- NULL
  
  # S matrix
  S <- readRDS("data/tourism_S.rds")
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
if (grepl("labour", data_label)){
  # Formalize data (Index, Time, Series1, ..., Seriesn)
  dat <- readRDS("data/labour_data.rds")
  h <- 12
  
  # Data details
  k <- sub('.*_', '', data_label) |> as.numeric() # 1:12
  freq <- 12
  start_train <- c(2010, 1)
  data_eg <- dat |>
    filter(Index == 1) |>
    select(!c(Index, Time)) |>
    ts(frequency = freq, start = start_train)
  
  length_train <- 163  - h - k + 1
  end_train <- time(data_eg)[length_train]
  start_test <- time(data_eg)[length_train + 1]
  end_test <- time(data_eg)[length_train + h]
  
  # Forecasting method
  fmethod <- "ets"
  special <- NULL
  
  # S matrix
  S <- readRDS("data/labour_S.rds")
}

#################################################
# Generate base forecasts
#################################################
indices <- dat |>
  distinct(Index) |>
  pull(Index) # all indices

fits <- resids <- train <- basefc <- test <- data.frame()
pb <- lazybar::lazyProgressBar(length(indices))
for (index in indices){
  # Hierarchical time series
  data_index <- dat |>
    filter(Index == index) |>
    select(!c(Index, Time)) |>
    ts(frequency = freq, start = start_train)
  train_index <- window(data_index, start = start_train, end = end_train)
  test_index <- window(data_index, start = start_test, end = end_test)
  
  # Base forecasts
  basef <- base_forecast(train_index, method = fmethod, h = NROW(test_index), special = special)
  
  # Results
  fits <- rbind(fits, 
                sapply(basef, function(l) l$fitted) |> 
                  data.frame(Index = index))
  resids <- rbind(resids, 
                  sapply(basef, function(l) l$residuals) |> 
                    data.frame(Index = index))
  train <- rbind(train, 
                 sapply(basef, function(l) l$train) |> 
                   data.frame(Index = index))
  basefc_index <- sapply(basef, function(l) l$forecast)
  if (is.vector(basefc_index)){
    basefc_index <- t(basefc_index)
  }
  basefc <- rbind(basefc, 
                  basefc_index |> 
                    data.frame(Index = index))
  test <- rbind(test, 
                test_index |> data.frame(Index = index))
  
  pb$tick()$print()
}

#################################################
# Save base forecast results
#################################################
for (i in c("S", "fits", "resids", "train", "basefc","test")){
  saveRDS(get(i), file = paste0("data/", data_label, "_", i, ".rds"))
}

