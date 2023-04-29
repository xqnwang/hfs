library(tidyverse)
library(magrittr)
library(tsibble)
library(forecast)

#################################################
# Import data
#################################################

data_sim <- readr::read_csv("data/simulation_data.csv") |>
  mutate(Time = yearquarter(Quarter),
         A = AA + AB,
         B = BA + BB,
         Total = AA + AB + BA + BB) |>
  select(c(Index, Time, AA, AB, BA, BB, A, B, Total)) 

S <- rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
           diag(rep(1, 4)))

indices <- data_sim |>
  distinct(Index) |>
  pull(Index)

base_forecast <- function(hts, method, h){
  hts |> 
    apply(2, function(series){
      if (method == 'arima'){
        model <- auto.arima(series)
      }
      if (method == 'ets'){
        model <- ets(series)
      }
      fc <- forecast(model, h)$mean
      list(fitted = model$fitted, residuals = model$residuals, 
           train = model$x, forecast = fc)
    })
}

fits <- resids <- train <- basefc <- test <- data.frame()
pb <- lazybar::lazyProgressBar(length(indices))
for (index in indices){
  # Hierarchical time series
  data_index <- data_sim |>
    filter(Index == index) |>
    select(!c(Index, Time)) |>
    ts(frequency = 4, start = c(1978, 1))
  train_index <- window(data_index, start = c(1978, 1), end = c(2018, 4))
  test_index <- window(data_index, start = c(2019, 1), end = c(2022, 4))
  
  # Base forecasts
  basef <- base_forecast(train_index, method = "ets", h = NROW(test_index))
  
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

saveRDS(fits, file = "data/simulation_fits.rds")
saveRDS(resids, file = "data/simulation_resids.rds")
saveRDS(train, file = "data/simulation_train.rds")
saveRDS(basefc, file = "data/simulation_basefc.rds")
saveRDS(test, file = "data/simulation_test.rds")
