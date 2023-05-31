library(tidyverse)
library(magrittr)
library(future)
library(forecast)
library(latex2exp)

# Utility functions
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
  err <- subset(data, select = -Index) - subset(fc, select = -Index)
  err <- cbind(err, Index = subset(data, select = Index))
  rmse <- err |> 
    as_tibble() |> 
    group_by(Index) |> 
    mutate(Horizon = row_number()) |> 
    filter(Horizon <= h) |>
    summarise_at(1:NCOL(err), function(x) sqrt(mean(x^2))) |>
    ungroup() |>
    select(!c("Index", "Horizon")) |>
    summarise_all(mean)
  return(rmse)
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
data_label <- "simulation"
scenario <- NULL
if (is.null(scenario)){
  reconsf <- readRDS(file = paste0("data/", data_label, "_reconsf.rds"))
} else{
  reconsf <- readRDS(file = paste0("data/", data_label, "_reconsf_", scenario, ".rds"))
}
test <- readRDS(file = paste0("data/", data_label, "_test.rds"))

# Structure information used to calculate RMSE across levels
top <- 1
middle <- 2:3
bottom <- 4:7
avg <- 1:7
horizon <- c(1, 4, 8, 16)

#----------------------------------------------------------------------
# Australian domestic tourism (only considering hierarchical structure)
## Quarterly series from 1998Q1-2017Q4: 80 observations for each series
##
## Total/State: 2 levels, n = 9
## Training set:  1998Q1-2015Q4
## Test set:      2016Q1-2017Q4
#----------------------------------------------------------------------
data_label <- "tourism"
scenario <- NULL
reconsf <- readRDS(file = paste0("data/", data_label, "_reconsf.rds"))
test <- readRDS(file = paste0("data/", data_label, "_test.rds"))

# Structure information used to calculate RMSE across levels
top <- 1
state <- 2:8
zone <- 9:35
region <- 36:111
avg <- 1:111
horizon <- c(1, 4, 8, 12)

#################################################
# Extract reconciled forecasts
#################################################
methods <- c("Base", "BU", 
             "OLS", "OLS_subset", "OLS_lasso",
             "WLSs", "WLSs_subset", "WLSs_lasso",
             "WLSv", "WLSv_subset", "WLSv_lasso",
             "MinT", "MinT_subset", "MinT_lasso",
             "MinTs", "MinTs_subset", "MinTs_lasso",
             "Emp_lasso")
subset_methods <- c("OLS_subset", "OLS_lasso", 
                    "WLSs_subset", "WLSs_lasso", 
                    "WLSv_subset", "WLSv_lasso",
                    "MinT_subset", "MinT_lasso",
                    "MinTs_subset", "MinTs_lasso",
                    "Emp_lasso")
indices <- unique(test$Index)

for(method in methods) { 
  out <- indices |> 
    purrr::map(\(index) 
               extract_element(data = reconsf, index = index, 
                               method = method, element = "y_tilde")) %>% 
    do.call(rbind, .)
  assign(tolower(method), out)
}

#################################################
# Calculate RMSE values
#################################################
for(h in horizon){
  out <- bind_rows(Base = calc_rmse(base, test, h = h),
                   BU = calc_rmse(bu, test, h = h),
                   OLS = calc_rmse(ols, test, h = h),
                   OLS_subset = calc_rmse(ols_subset, test, h = h),
                   OLS_lasso = calc_rmse(ols_lasso, test, h = h),
                   WLSs = calc_rmse(wlss, test, h = h),
                   WLSs_subset = calc_rmse(wlss_subset, test, h = h),
                   WLSs_lasso = calc_rmse(wlss_lasso, test, h = h),
                   WLSv = calc_rmse(wlsv, test, h = h),
                   WLSv_subset = calc_rmse(wlsv_subset, test, h = h),
                   WLSv_lasso = calc_rmse(wlsv_lasso, test, h = h),
                   MinT = calc_rmse(mint, test, h = h),
                   MinT_subset = calc_rmse(mint_subset, test, h = h),
                   MinT_lasso = calc_rmse(mint_lasso, test, h = h),
                   MinTs = calc_rmse(mints, test, h = h),
                   MinTs_subset = calc_rmse(mints_subset, test, h = h),
                   MinTs_lasso = calc_rmse(mints_lasso, test, h = h),
                   Emp_lasso = calc_rmse(emp_lasso, test, h = h),
                   .id = "Method") |>
    rowwise() |>
    mutate(Top = mean(c_across(top + 1)),
           # State = mean(c_across(state + 1)),
           # Zone = mean(c_across(zone + 1)),
           # Region = mean(c_across(region + 1)),
           Middle = mean(c_across((middle + 1))),
           Bottom = mean(c_across((bottom + 1))),
           Average = mean(c_across(avg + 1))) |>
    select(Method, Top, Middle, Bottom, Average)
    # select(Method, Top, State, Zone, Region, Average)
  assign(paste0("rmse_h", h), out)
  # if (is.null(scenario)){
  #   saveRDS(out, file = paste0("data/", data_label, "_reconsf_rmse_", h, ".rds"))
  # } else{
  #   saveRDS(out, file = paste0("data/", data_label, "_reconsf_", scenario, "_rmse_", h, ".rds"))
  # }
}

#################################################
# Extract z
#################################################
z_summary <- NULL
for(method in subset_methods) { 
  out <- indices |> 
    purrr::map(\(index) 
               extract_element(data = reconsf, index = index, 
                               method = method, element = "z")) %>% 
    do.call(rbind, .)
  out <- ifelse(out < 1e-3, 0, out) |> cbind(Method = method)
  if (length(out) > 1){
    z_summary <- rbind(z_summary, out) 
  }
}
series_name <- c(colnames(test), "Method")
colnames(z_summary) <- series_name
if (is.null(scenario)){
  saveRDS(z_summary, file = paste0("data/", data_label, "_reconsf_z_summary.rds"))
} else{
  saveRDS(z_summary, file = paste0("data/", data_label, "_reconsf_", scenario, "_z_summary.rds"))
}

z_summary |>
  as_tibble() |>
  group_by(Method) |>
  summarise_at(1:(NCOL(test)-1), function(x) sum(x==0)) |>
  pivot_longer(
    cols = 2:NCOL(test),
    names_to = "Series",
    values_to = "Frequency") |>
  ggplot(aes(x = factor(Series, levels = series_name), y = Frequency)) +
  geom_bar(stat = "identity") +
  facet_grid(vars(Method), scales = "free_y") +
  labs(title = "Frequency of being zeroed out",
       x = "",
       y= "")

#################################################
# Extract lambda_report
#################################################
lambda_summary <- NULL
for(method in subset_methods) { 
  out <- indices |> 
    purrr::map(\(index) 
               extract_element(data = reconsf, index = index, 
                               method = method, element = "lambda_report")) %>% 
    do.call(rbind, .)
  out <- cbind(out, Method = method)
  lambda_summary <- rbind(lambda_summary, out)
}
if (is.null(scenario)){
  saveRDS(lambda_summary, file = paste0("data/", data_label, "_reconsf_lambda_summary.rds"))
} else{
  saveRDS(lambda_summary, file = paste0("data/", data_label, "_reconsf_", scenario, "_lambda_summary.rds"))
}

lambda_summary |>
  group_by(Method, Index) |>
  summarise(sse_index = which.min(sse)) |>
  ungroup() |>
  group_by(Method) |>
  ggplot(aes(x = factor(sse_index))) +
  geom_bar(stat = "count") +
  facet_grid(vars(Method), scales = "free_y") +
  labs(title = TeX(r"(Frequency of being selected as the optimal $\lambda_0$)"),
       x = TeX(r"(index of $\lambda_0$, $\lambda_{0\min} = 0$)"),
       y= "")

