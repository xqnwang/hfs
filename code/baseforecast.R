library(tidyverse)
library(magrittr)
library(future)
library(forecast)

# Utility function
base_forecast <- function(hts, method, h){
  out <- list()
  for(i in 1:NCOL(hts)){
    series <- hts[, i]
    if (method == 'arima'){
      model <- auto.arima(series)
    }
    if (method == 'ets'){
      model <- ets(series)
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

# S matrix
S <- rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
           diag(rep(1, 4)))

#----------------------------------------------------------------------
# Australian domestic tourism (only considering hierarchical structure)
## Quarterly series from 1998Q1-2017Q4: 80 observations for each series
##
## Total/State: 2 levels, n = 9
## Training set:  1998Q1-2015Q4
## Test set:      2016Q1-2017Q4
#----------------------------------------------------------------------
# Formalize data (Index, Time, Series1, ..., Seriesn)
dat <- readRDS("data/tourism_data.rds")

# Data details
data_label <- "tourism" # used to name the saved results
freq <- 4
start_train <- c(1998, 1)
end_train <- c(2015, 4)
start_test <- c(2016, 1)
end_test <- c(2017, 4)

# Forecasting method
fmethod <- "ets"

# S matrix
S <- rbind(rep(1, 8), diag(rep(1, 8)))

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
  basef <- base_forecast(train_index, method = fmethod, h = NROW(test_index))
  
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
for (i in c("fits", "resids", "train", "basefc","test")){
  saveRDS(get(i), file = paste0("data/", data_label, "_", i, ".rds"))
}

