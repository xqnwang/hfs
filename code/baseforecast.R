library(tidyverse)
library(magrittr)
library(future)
library(forecast)

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
# Formalize data (Index, Time, Series1, ..., Seriesn)
dat <- readr::read_csv("data/simulation_data.csv") |>
  mutate(Time = tsibble::yearquarter(Quarter),
         A = AA + AB,
         B = BA + BB,
         Total = AA + AB + BA + BB) |>
  select(c(Index, Time, Total, A, B, AA, AB, BA, BB)) 

# Data details
data_label <- "simulation" # used to name the saved results
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

#----------------------------------------------------------------------
# Simulation data - correlation
## Total/Middle/Bottom: 3 levels, n = 7
## Training set:  1-100
## Test set:      1
#----------------------------------------------------------------------
k <- 1 # 1:9
data_label <- paste0("corr_", k)
freq <- 1
start_train <- 1
end_train <- 100
start_test <- 101
end_test <- 101
dat <- readr::read_csv(paste0("data/corr_", k, "_data.csv")) |>
  mutate(A = AA + AB,
         B = BA + BB,
         Total = AA + AB + BA + BB) |>
  select(c(Index, Time, Total, A, B, AA, AB, BA, BB))

# Forecasting method
fmethod <- "ar"
special <- c("A", "BA")

# S matrix
S <- rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
           diag(rep(1, 4)))

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
# Formalize data (Index, Time, Series1, ..., Seriesn)
dat <- readRDS("data/tourism_data.rds")

# Data details
data_label <- "tourism" # used to name the saved results
freq <- 12
start_train <- c(1998, 1)
end_train <- c(2016, 12)
start_test <- c(2017, 1)
end_test <- c(2017, 12)

# Forecasting method
fmethod <- "ets"
special <- NULL

# S matrix
S <- readRDS("data/tourism_S.rds")

#----------------------------------------------------------------------
# Australian prison population
##
## Total number of prisoners in Australia over the period 2005Q1–2016Q4
## Quarterly series: 48 quarters (12 years) for each series
##
## Gender * Legal * State: n = 81 series in total, nb = 32 series at the bottom level
##
## Training set:  2005Q1-2014Q4
## Test set:      2015Q1-2016Q4
#----------------------------------------------------------------------
# Formalize data (Index, Time, Series1, ..., Seriesn)
dat <- readRDS("data/prison_data.rds")

# Data details
data_label <- "prison" # used to name the saved results
freq <- 4
start_train <- c(2005, 1)
end_train <- c(2014, 4)
start_test <- c(2015, 1)
end_test <- c(2016, 4)

# Forecasting method
fmethod <- "ets"
special <- NULL

# S matrix
S <- readRDS("data/prison_S.rds")

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
# Formalize data (Index, Time, Series1, ..., Seriesn)
dat <- readRDS("data/labour_data.rds")

# Data details
data_label <- "labour" # used to name the saved results
freq <- 12
start_train <- c(2010, 1)
end_train <- c(2022, 7)
start_test <- c(2022, 8)
end_test <- c(2023, 7)

# Forecasting method
fmethod <- "ets"
special <- NULL

# S matrix
S <- readRDS("data/labour_S.rds")

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
  basefc <- rbind(basefc, 
                  sapply(basef, function(l) l$forecast) |> 
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

