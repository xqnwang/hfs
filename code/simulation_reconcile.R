library(tidyverse)
library(magrittr)
library(future)
library(forecast)

source('R/sourceDir.R')
sourceDir("R", recursive = TRUE)

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
    out[[i]] <- list(fitted = model$fitted, residuals = model$residuals, 
         train = model$x, forecast = fc)
  }
  names(out) <- colnames(hts)
  out
}

reconcile_forecast <- function(index, fits, train, basefc, resids, test, 
                               deteriorate = FALSE, deteriorate_series, deteriorate_rate,
                               G_bench, nlambda_0 = 20, parallel = FALSE, workers = 8){
  
  fitted_values <- fits[fits$Index == index, 1:7] |> as.matrix()
  residuals <- resids[resids$Index == index, 1:7] |> as.matrix()
  train_data <- train[train$Index == index, 1:7] |> as.matrix()
  base_forecasts <- basefc[basefc$Index == index, 1:7] |> as.matrix()
  
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
                          parallel = parallel, workers = workers)
  WLSs <- reconcile(base_forecasts = base_forecasts, S = S, method = "wls_struct")
  WLSs_subset <- reconcile(base_forecasts = base_forecasts, S = S, method = "wls_struct",
                           residuals = residuals, 
                           fitted_values = fitted_values, train_data = train_data,
                           subset = TRUE, G_bench = G_bench, nlambda_0 = nlambda_0,
                           parallel = parallel, workers = workers)
  WLSv <- reconcile(base_forecasts = base_forecasts, S = S, method = "wls_var",
                    residuals = residuals)
  WLSv_subset <- reconcile(base_forecasts = base_forecasts, S = S, method = "wls_var",
                           residuals = residuals, 
                           fitted_values = fitted_values, train_data = train_data,
                           subset = TRUE, G_bench = G_bench, nlambda_0 = nlambda_0,
                           parallel = parallel, workers = workers)
  MinT <- reconcile(base_forecasts = base_forecasts, S = S, method = "mint_cov",
                    residuals = residuals)
  MinT_subset <- reconcile(base_forecasts = base_forecasts, S = S, method = "mint_cov",
                           residuals = residuals,
                           fitted_values = fitted_values, train_data = train_data,
                           subset = TRUE, G_bench = G_bench, nlambda_0 = nlambda_0,
                           parallel = parallel, workers = workers)
  MinTs <- reconcile(base_forecasts = base_forecasts, S = S, method = "mint_shrink",
                     residuals = residuals)
  MinTs_subset <- reconcile(base_forecasts = base_forecasts, S = S, method = "mint_shrink",
                            residuals = residuals, 
                            fitted_values = fitted_values, train_data = train_data,
                            subset = TRUE, G_bench = G_bench, nlambda_0 = nlambda_0,
                            parallel = parallel, workers = workers)
  
  list(Base = list(y_tilde = base_forecasts, G = NA, z = NA, lambda0_report = NA), 
       BU = BU, 
       OLS = OLS, OLS_subset = OLS_subset,
       WLSs = WLSs, WLSs_subset = WLSs_subset,
       WLSv = WLSv, WLSv_subset = WLSv_subset,
       MinT = MinT, MinT_subset = MinT_subset,
       MinTs = MinTs, MinTs_subset = MinTs_subset)
}

extract_element <- function(data, index, method, element){
  out <- data[[index]][[method]][[element]]
  if(length(out) == 1){
    if(is.na(out)){
      out <- NULL 
    }
  } else{
    if(is.vector(out)){
      out <- cbind(matrix(out, nrow = 1), Index = index) 
    }else{
      out <- cbind(out, Index = index) 
    }
  }
  out
}

calc_rmse <- function(fc, data, h){
  err <- subset(data, select=-Index) - subset(fc, select=-Index)
  err <- cbind(err, Index = subset(data, select=Index))
  rmse <- err |> 
    as_tibble() |> 
    group_by(Index) |> 
    mutate(Horizon = row_number()) |> 
    filter(Horizon <= h) |>
    summarise_at(vars(Total:BB), function(x) sqrt(mean(x^2))) |>
    ungroup() |>
    summarise_at(vars(Total:BB), mean)
  return(rmse)
}

#################################################
# Import data
#################################################
data_sim <- readr::read_csv("data/simulation_data_noisy.csv") |>
  mutate(Time = tsibble::yearquarter(Quarter),
         A = AA + AB,
         B = BA + BB,
         Total = AA + AB + BA + BB) |>
  select(c(Index, Time, Total, A, B, AA, AB, BA, BB)) 

S <- rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
           diag(rep(1, 4)))

indices <- data_sim |>
  distinct(Index) |>
  pull(Index)

#################################################
# Generate base forecasts
#################################################
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

#################################################
# Save base forecast results
#################################################
saveRDS(fits, file = "data/simulation_noisy_fits.rds")
saveRDS(resids, file = "data/simulation_noisy_resids.rds")
saveRDS(train, file = "data/simulation_noisy_train.rds")
saveRDS(basefc, file = "data/simulation_noisy_basefc.rds")
saveRDS(test, file = "data/simulation_noisy_test.rds")

#################################################
# Import base forecast results
#################################################
fits <- readRDS("data/simulation_noisy_fits.rds")
resids <- readRDS("data/simulation_noisy_resids.rds")
train <- readRDS("data/simulation_noisy_train.rds")
basefc <- readRDS("data/simulation_noisy_basefc.rds")
test <- readRDS("data/simulation_noisy_test.rds")

#################################################
# Reconcile ets forecasts
#################################################
plan(multisession, workers = 8)
reconsf <- indices |>
  furrr::future_map(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, 
                                                deteriorate = FALSE,
                                                G_bench = "Zero", nlambda_0 = 20,
                                                parallel = FALSE, workers = 8),
                    .progress = TRUE)
reconsf <- saveRDS(reconsf, file = "data/simulation_noisy_reconsf.rds")
rm(reconsf)

#################################################
# Reconcile ets forecasts - Scenario 1: base forecasts of AA and B are poor
#################################################
plan(multisession, workers = 8)
reconsf_s1 <- indices |>
  furrr::future_map(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, 
                                                deteriorate = TRUE, 
                                                deteriorate_series = c("AA", "B"),
                                                deteriorate_rate = 1.5,
                                                G_bench = "Zero", nlambda_0 = 20, 
                                                parallel = FALSE, workers = 8),
                    .progress = TRUE)
reconsf_s1 <- saveRDS(reconsf_s1, file = "data/simulation_noisy_reconsf_s1.rds")
rm(reconsf_s1)

#################################################
# Reconcile ets forecasts - Scenario 2: base forecasts of Total are poor
#################################################
plan(multisession, workers = 8)
reconsf_s2 <- indices |>
  furrr::future_map(\(index) reconcile_forecast(index, fits, train, basefc, resids, test, 
                                                deteriorate = TRUE, 
                                                deteriorate_series = c("Total"), 
                                                deteriorate_rate = 1.5,
                                                G_bench = "Zero", nlambda_0 = 20, 
                                                parallel = FALSE, workers = 8),
                    .progress = TRUE)
reconsf_s2 <- saveRDS(reconsf_s2, file = "data/simulation_noisy_reconsf_s2.rds")
rm(reconsf_s2)

#################################################
# Evaluation
#################################################
reconsf <- readRDS("data/simulation_noisy_reconsf_s2.rds")

# Extract reconciled forecasts
methods <- c("Base", "BU", "OLS", "OLS_subset", 
             "WLSs", "WLSs_subset", "WLSv", "WLSv_subset", 
             "MinT", "MinT_subset", 
             "MinTs", "MinTs_subset")
for(method in methods) { 
  out <- indices |> 
    purrr::map(\(index) 
               extract_element(data = reconsf, index = index, 
                               method = method, element = "y_tilde")) %>% 
    do.call(rbind, .)
  assign(tolower(method), out)
}

# Calculate RMSE values
for(h in c(1, 8, 16)){
  out <- bind_rows(Base = calc_rmse(base, test, h = h),
                   BU = calc_rmse(bu, test, h = h),
                   OLS = calc_rmse(ols, test, h = h),
                   OLS_subset = calc_rmse(ols_subset, test, h = h),
                   WLSs = calc_rmse(wlss, test, h = h),
                   WLSs_subset = calc_rmse(wlss_subset, test, h = h),
                   WLSv = calc_rmse(wlsv, test, h = h),
                   WLSv_subset = calc_rmse(wlsv_subset, test, h = h),
                   MinT = calc_rmse(mint, test, h = h),
                   MinT_subset = calc_rmse(mint_subset, test, h = h),
                   MinTs = calc_rmse(mints, test, h = h),
                   MinTs_subset = calc_rmse(mints_subset, test, h = h),
                   .id = "Method") |>
    mutate(Top = Total,
           Middle = (A + B)/2,
           Bottom = (AA + AB + BA + BB)/4,
           Average = (Total + A + B + AA + AB + BA + BB)/7) |>
    select(Method, Top, Middle, Bottom, Average)
  assign(paste0("rmse_h", h), out)
}

# Extract z
methods <- c("OLS_subset", "WLSs_subset", "WLSv_subset", 
             "MinT_subset", "MinTs_subset")
z_summary <- NULL
for(method in methods) { 
  out <- indices |> 
    purrr::map(\(index) 
               extract_element(data = reconsf, index = index, 
                               method = method, element = "z")) %>% 
    do.call(rbind, .)
  out <- ifelse(out < 1e-3, 0, out) |> cbind(Method = method)
  z_summary <- rbind(z_summary, out)
}
colnames(z_summary) <- c("Total", "A", "B", "AA", "AB", "BA", "BB", "Index", "Method")

z_summary |>
  as_tibble() |>
  group_by(Method) |>
  summarise_at(vars(Total:BB), function(x) sum(x==0)) |>
  pivot_longer(
    cols = Total:BB,
    names_to = "Series",
    values_to = "Frequency") |>
  ggplot(aes(x = Series, y = Frequency)) +
  geom_bar(stat="identity") +
  facet_grid(vars(Method), scales = "free_y") +
  labs(title = "Frequency of zero lambda_0",
       y= "")

# Extract lambda0_report
methods <- c("OLS_subset", "WLSs_subset", "WLSv_subset", 
             "MinT_subset", "MinTs_subset")
lambda0_summary <- NULL
for(method in methods) { 
  out <- indices |> 
    purrr::map(\(index) 
               extract_element(data = reconsf, index = index, 
                               method = method, element = "lambda0_report")) %>% 
    do.call(rbind, .)
  out <- cbind(out, Method = method)
  lambda0_summary <- rbind(lambda0_summary, out)
}

lambda0_summary |>
  group_by(Method, Index) |>
  summarise(sse_index = which.min(sse)) |>
  ungroup() |>
  group_by(Method) |>
  ggplot(aes(x = factor(sse_index))) +
  geom_bar(stat = "count") +
  facet_grid(vars(Method), scales = "free_y") +
  labs(title = "Frequency of zero lambda_0",
       y= "")

