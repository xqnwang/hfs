library(fable)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(purrr)
library(tsibble)
library(rlang)
library(ggplot2)

source('R/sourceDir.R')
sourceDir("R", recursive = TRUE)

# Import simulation data
data_sim <- readr::read_csv("data/simulation_data.csv") |>
  dplyr::mutate(Quarter = yearquarter(Quarter)) |>
  tidyr::pivot_longer(cols = c("AA", "AB", "BA", "BB"),
               names_to = "Bottom",
               values_to = "Value") |>
  dplyr::mutate(Middle = str_sub(Bottom, 1, 1)) |>
  select(c(Index, Quarter, Middle, Bottom, Value)) |>
  as_tsibble(key = c(Index, Middle, Bottom),
             index = Quarter) |>
  relocate(Quarter)

data_repeats <- purrr::map(distinct(data_sim, Index) |> pull(),
                           function(i){
                             return(data_sim |> 
                                      filter(Index == i))
                           })

# Scenario 1: models for some bottom-level series (AA and BB) are mis-specified
forecast_repeat_scenario1 <- function(data_r, subset, lasso, ridge, 
                            lambda_0 = NULL, lambda_1 = 0, lambda_2 = 0,
                            nlambda_0 = 20){
  # hts
  data_hts <- data_r |>
    aggregate_key(Middle/Bottom, Value = sum(Value))
  
  # visualization
  # data_hts |>
  #   autoplot(Value) +
  #   labs(title = "Simulation data",
  #        x = "Quarter") +
  #   facet_wrap(vars(Bottom), scales = "free_y", ncol = 3) +
  #   theme(legend.position = "none")
  
  # base modeling
  fit_base_drift <- data_hts |>
    filter(year(Quarter) <= 2018) |>
    filter(Bottom == "AA" | Bottom == "BB" ) |>
    model(base = RW(Value ~ drift()))
  fit_base_ets <- data_hts |>
    filter(year(Quarter) <= 2018) |>
    filter(Bottom != "AA" & Bottom != "BB" ) |>
    model(base = ETS(Value))
  fit_base <- dplyr::bind_rows(fit_base_drift, fit_base_ets)
  
  # reconciliation
  fit_rec <- fit_base |>
    reconcile(
      BU = bottom_up(base),
      OLS = min_trace(base, method = "ols"),
      OLS_subset = min_trace_subset(base, method = "ols", 
                                    subset, lasso, ridge, 
                                    lambda_0, lambda_1, lambda_2, 
                                    nlambda_0, sparse = FALSE),
      WLSv = min_trace(base, method = "wls_var"),
      WLSv_subset = min_trace_subset(base, method = "wls_var", 
                                     subset, lasso, ridge, 
                                     lambda_0, lambda_1, lambda_2, 
                                     nlambda_0, sparse = FALSE),
      WLSs = min_trace(base, method = "wls_struct"),
      WLSs_subset = min_trace_subset(base, method = "wls_struct", 
                                     subset, lasso, ridge, 
                                     lambda_0, lambda_1, lambda_2, 
                                     nlambda_0, sparse = FALSE),
      MinTc = min_trace(base, method = "mint_cov"),
      MinTc_subset = min_trace_subset(base, method = "mint_cov", 
                                      subset, lasso, ridge, 
                                      lambda_0, lambda_1, lambda_2, 
                                      nlambda_0, sparse = FALSE),
      MinTs = min_trace(base, method = "mint_shrink"),
      MinTs_subset = min_trace_subset(base, method = "mint_shrink", 
                                     subset, lasso, ridge, 
                                     lambda_0, lambda_1, lambda_2, 
                                     nlambda_0, sparse = FALSE)
    )
  
  # forecast
  fc <- fit_rec |>
    forecast(h = "4 years") |>
    group_by(Middle, Bottom, .model) |>
    mutate(h = row_number()) |>
    ungroup()
  
  # accuracy
  error <- fc |>
    as_fable(response = "Value", distribution = Value) |>
    accuracy(data = data_hts, by = c("h", ".model", "Middle", "Bottom")) |>
    select(-ACF1)
  
  return(error)
}

all_accuracy <- purrr::map(data_repeats, forecast_repeat_scenario1, 
                           subset = TRUE, lasso = FALSE, ridge = FALSE, 
                           nlambda_0 = 20, .progress = TRUE)
saveRDS(all_accuracy, file = "data/simulation_accuracy_scenario1.rds")
# all_accuracy <- readRDS("data/simulation_accuracy_scenario1.rds")
rm(all_accuracy)

# Scenario 2: models for top-level series (Total) are mis-specified
forecast_repeat_scenario2 <- function(data_r, subset, lasso, ridge, 
                            lambda_0 = NULL, lambda_1 = 0, lambda_2 = 0,
                            nlambda_0 = 20){
  # hts
  data_hts <- data_r |>
    aggregate_key(Middle/Bottom, Value = sum(Value))
  
  # visualization
  # data_hts |>
  #   autoplot(Value) +
  #   labs(title = "Simulation data",
  #        x = "Quarter") +
  #   facet_wrap(vars(Bottom), scales = "free_y", ncol = 3) +
  #   theme(legend.position = "none")
  
  # base modeling
  fit_base_drift <- data_hts |>
    filter(year(Quarter) <= 2018) |>
    filter(is_aggregated(Middle) & is_aggregated(Bottom)) |>
    model(base = RW(Value ~ drift()))
  fit_base_ets <- data_hts |>
    filter(year(Quarter) <= 2018) |>
    filter(!is_aggregated(Middle) | !is_aggregated(Bottom)) |>
    model(base = ETS(Value))
  fit_base <- dplyr::bind_rows(fit_base_drift, fit_base_ets)
  
  # reconciliation
  fit_rec <- fit_base |>
    reconcile(
      BU = bottom_up(base),
      OLS = min_trace(base, method = "ols"),
      OLS_subset = min_trace_subset(base, method = "ols",
                                    subset, lasso, ridge,
                                    lambda_0, lambda_1, lambda_2,
                                    nlambda_0, sparse = FALSE),
      WLSv = min_trace(base, method = "wls_var"),
      WLSv_subset = min_trace_subset(base, method = "wls_var",
                                     subset, lasso, ridge,
                                     lambda_0, lambda_1, lambda_2,
                                     nlambda_0, sparse = FALSE),
      WLSs = min_trace(base, method = "wls_struct"),
      WLSs_subset = min_trace_subset(base, method = "wls_struct",
                                     subset, lasso, ridge,
                                     lambda_0, lambda_1, lambda_2,
                                     nlambda_0, sparse = FALSE),
      MinTc = min_trace(base, method = "mint_cov"),
      MinTc_subset = min_trace_subset(base, method = "mint_cov",
                                      subset, lasso, ridge,
                                      lambda_0, lambda_1, lambda_2,
                                      nlambda_0, sparse = FALSE),
      MinTs = min_trace(base, method = "mint_shrink"),
      MinTs_subset = min_trace_subset(base, method = "mint_shrink", 
                                      subset, lasso, ridge, 
                                      lambda_0, lambda_1, lambda_2, 
                                      nlambda_0, sparse = FALSE)
    )
  
  # forecast
  fc <- fit_rec |>
    forecast(h = "4 years") |>
    group_by(Middle, Bottom, .model) |>
    mutate(h = row_number()) |>
    ungroup()
  
  # accuracy
  error <- fc |>
    as_fable(response = "Value", distribution = Value) |>
    accuracy(data = data_hts, by = c("h", ".model", "Middle", "Bottom")) |>
    select(-ACF1)
  
  return(error)
}

all_accuracy <- purrr::map(data_repeats, forecast_repeat_scenario2, 
                           subset = TRUE, lasso = FALSE, ridge = FALSE, 
                           nlambda_0 = 20, .progress = TRUE)
saveRDS(all_accuracy, file = "data/simulation_accuracy_scenario2.rds")
# all_accuracy <- readRDS("data/simulation_accuracy_scenario2.rds")

# Calculate accuracy
sim_accuracy <- bind_rows(all_accuracy, .id = "id")

model_name <- sim_accuracy |> distinct(.model) |> pull(.model)
filter_tab <- matrix(NA, ncol = 1, nrow = 4)
filter_tab[1] <- "sim_accuracy |> filter(is_aggregated(Middle) & is_aggregated(Bottom))"
filter_tab[2] <- "sim_accuracy |> filter(!is_aggregated(Middle) & is_aggregated(Bottom))"
filter_tab[3] <- "sim_accuracy |> filter(!is_aggregated(Middle) & !is_aggregated(Bottom))"
filter_tab[4] <- "sim_accuracy"
accuracy_summary <- lapply(1:4, function(i){
  err <- eval(parse(text = filter_tab[i])) |>
    group_by(.model) |>
    summarise_at(c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "RMSSE"), mean)
  err <- err |> arrange(factor(.model, levels = model_name))
})
names(accuracy_summary) <- c("Total", "Middle", "Bottom", "All")

accuracy_summary_1 <- lapply(1:4, function(i){
  err <- eval(parse(text = filter_tab[i])) |>
    filter(h == 1) |>
    group_by(.model) |>
    summarise_at(c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "RMSSE"), mean)
  err <- err |> arrange(factor(.model, levels = model_name))
})
names(accuracy_summary_1) <- c("Total", "Middle", "Bottom", "All")

accuracy_summary_1_8 <- lapply(1:4, function(i){
  err <- eval(parse(text = filter_tab[i])) |>
    filter(h <= 8) |>
    group_by(.model) |>
    summarise_at(c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "RMSSE"), mean)
  err <- err |> arrange(factor(.model, levels = model_name))
})
names(accuracy_summary_1_8) <- c("Total", "Middle", "Bottom", "All")

## RMSE
RMSE_total <- cbind(accuracy_summary_1$Total |> pull(RMSE), 
                    accuracy_summary_1_8$Total |> pull(RMSE),
                    accuracy_summary$Total |> pull(RMSE))
RMSE_middle <- cbind(accuracy_summary_1$Middle |> pull(RMSE), 
                     accuracy_summary_1_8$Middle |> pull(RMSE),
                     accuracy_summary$Middle |> pull(RMSE))
RMSE_bottom <- cbind(accuracy_summary_1$Bottom |> pull(RMSE), 
                     accuracy_summary_1_8$Bottom |> pull(RMSE),
                     accuracy_summary$Bottom |> pull(RMSE))
RMSE_all <- cbind(accuracy_summary_1$All |> pull(RMSE), 
                  accuracy_summary_1_8$All |> pull(RMSE),
                  accuracy_summary$All |> pull(RMSE))

## MASE
MASE_total <- cbind(accuracy_summary_1$Total |> pull(MASE), 
                    accuracy_summary_1_8$Total |> pull(MASE),
                    accuracy_summary$Total |> pull(MASE))
MASE_middle <- cbind(accuracy_summary_1$Middle |> pull(MASE), 
                     accuracy_summary_1_8$Middle |> pull(MASE),
                     accuracy_summary$Middle |> pull(MASE))
MASE_bottom <- cbind(accuracy_summary_1$Bottom |> pull(MASE), 
                     accuracy_summary_1_8$Bottom |> pull(MASE),
                     accuracy_summary$Bottom |> pull(MASE))
MASE_all <- cbind(accuracy_summary_1$All |> pull(MASE), 
                  accuracy_summary_1_8$All |> pull(MASE),
                  accuracy_summary$All |> pull(MASE))

total <- data.frame(Method = model_name, cbind(RMSE_total, MASE_total) |> round(3))
middle <- data.frame(Method = model_name, cbind(RMSE_middle, MASE_middle) |> round(3))
bottom <- data.frame(Method = model_name, cbind(RMSE_bottom, MASE_bottom) |> round(3))
all <- data.frame(Method = model_name, cbind(RMSE_all, MASE_all) |> round(3))

library(flextable)
table_total <- flextable(total) |>
  add_header_row(
    top = TRUE,        # New header goes on top of existing header row
    values = c("Total",     # Header values for each column below
               "RMSE", 
               "",     # This will be the top-level header for this and two next columns
               "",
               "MASE",
               "", # This will be the top-level header for this and two next columns
               "")) |>
  set_header_labels(         # Rename the columns in original header row
    X1 = "h=1", 
    X2 = "1-8",                  
    X3 = "1-16",
    X4 = "h=1",
    X5 = "1-8",
    X6 = "1-16") |>
  bg(part = "body", bg = "gray95") |>
  bg(j = 2, i = ~ X1 == min(X1), part = "body", bg = "orange") |>
  bg(j = 3, i = ~ X2 == min(X2), part = "body", bg = "orange") |>
  bg(j = 4, i = ~ X3 == min(X3), part = "body", bg = "orange") |>
  bg(j = 5, i = ~ X4 == min(X4), part = "body", bg = "orange") |>
  bg(j = 6, i = ~ X5 == min(X5), part = "body", bg = "orange") |>
  bg(j = 7, i = ~ X6 == min(X6), part = "body", bg = "orange")

table_middle <- flextable(middle) |>
  add_header_row(
    top = TRUE,        # New header goes on top of existing header row
    values = c("Middle",     # Header values for each column below
               "RMSE", 
               "",     # This will be the top-level header for this and two next columns
               "",
               "MASE",
               "", # This will be the top-level header for this and two next columns
               "")) |>
  set_header_labels(         # Rename the columns in original header row
    X1 = "h=1", 
    X2 = "1-8",                  
    X3 = "1-16",
    X4 = "h=1",
    X5 = "1-8",
    X6 = "1-16") |>
  bg(part = "body", bg = "gray95") |>
  bg(j = 2, i = ~ X1 == min(X1), part = "body", bg = "orange") |>
  bg(j = 3, i = ~ X2 == min(X2), part = "body", bg = "orange") |>
  bg(j = 4, i = ~ X3 == min(X3), part = "body", bg = "orange") |>
  bg(j = 5, i = ~ X4 == min(X4), part = "body", bg = "orange") |>
  bg(j = 6, i = ~ X5 == min(X5), part = "body", bg = "orange") |>
  bg(j = 7, i = ~ X6 == min(X6), part = "body", bg = "orange")

table_bottom <- flextable(bottom) |>
  add_header_row(
    top = TRUE,        # New header goes on top of existing header row
    values = c("Bottom",     # Header values for each column below
               "RMSE", 
               "",     # This will be the top-level header for this and two next columns
               "",
               "MASE",
               "", # This will be the top-level header for this and two next columns
               "")) |>
  set_header_labels(         # Rename the columns in original header row
    X1 = "h=1", 
    X2 = "1-8",                  
    X3 = "1-16",
    X4 = "h=1",
    X5 = "1-8",
    X6 = "1-16") |>
  bg(part = "body", bg = "gray95") |>
  bg(j = 2, i = ~ X1 == min(X1), part = "body", bg = "orange") |>
  bg(j = 3, i = ~ X2 == min(X2), part = "body", bg = "orange") |>
  bg(j = 4, i = ~ X3 == min(X3), part = "body", bg = "orange") |>
  bg(j = 5, i = ~ X4 == min(X4), part = "body", bg = "orange") |>
  bg(j = 6, i = ~ X5 == min(X5), part = "body", bg = "orange") |>
  bg(j = 7, i = ~ X6 == min(X6), part = "body", bg = "orange")

table_all <- flextable(all) |>
  add_header_row(
    top = TRUE,        # New header goes on top of existing header row
    values = c("All",     # Header values for each column below
               "RMSE", 
               "",     # This will be the top-level header for this and two next columns
               "",
               "MASE",
               "", # This will be the top-level header for this and two next columns
               "")) |>
  set_header_labels(         # Rename the columns in original header row
    X1 = "h=1", 
    X2 = "1-8",                  
    X3 = "1-16",
    X4 = "h=1",
    X5 = "1-8",
    X6 = "1-16") |>
  bg(part = "body", bg = "gray95") |>
  bg(j = 2, i = ~ X1 == min(X1), part = "body", bg = "orange") |>
  bg(j = 3, i = ~ X2 == min(X2), part = "body", bg = "orange") |>
  bg(j = 4, i = ~ X3 == min(X3), part = "body", bg = "orange") |>
  bg(j = 5, i = ~ X4 == min(X4), part = "body", bg = "orange") |>
  bg(j = 6, i = ~ X5 == min(X5), part = "body", bg = "orange") |>
  bg(j = 7, i = ~ X6 == min(X6), part = "body", bg = "orange")

table_total; table_middle; table_bottom; table_all

sim_accuracy |> 
  filter(.model %in% c("OLS", "OLS_subset")) |>
  arrange(Middle, Bottom) |>
  select(!.type) |>
  group_by(id, Middle, Bottom) |>
  mutate(subset_better = c(0, diff(MASE))) |>
  ungroup() |>
  mutate(better = ifelse(subset_better < 0, TRUE, FALSE)) |>
  mutate(worse = ifelse(subset_better > 0, TRUE, FALSE)) |>
  select(better, worse) |>
  summarise(better = mean(better), worse = mean(worse))

sim_accuracy |> 
  filter(.model %in% c("WLSs", "WLSs_subset")) |>
  arrange(Middle, Bottom) |>
  select(!.type) |>
  group_by(id, Middle, Bottom) |>
  mutate(subset_better = c(0, diff(MASE))) |>
  ungroup() |>
  mutate(better = ifelse(subset_better < 0, TRUE, FALSE)) |>
  mutate(worse = ifelse(subset_better > 0, TRUE, FALSE)) |>
  select(better, worse) |>
  summarise(better = mean(better), worse = mean(worse))

sim_accuracy |> 
  filter(.model %in% c("WLSv", "WLSv_subset")) |>
  arrange(Middle, Bottom) |>
  select(!.type) |>
  group_by(id, Middle, Bottom) |>
  mutate(subset_better = c(0, diff(MASE))) |>
  ungroup() |>
  mutate(better = ifelse(subset_better < 0, TRUE, FALSE)) |>
  mutate(worse = ifelse(subset_better > 0, TRUE, FALSE)) |>
  select(better, worse) |>
  summarise(better = mean(better), worse = mean(worse))

sim_accuracy |> 
  filter(.model %in% c("MinTc", "MinTc_subset")) |>
  arrange(Middle, Bottom) |>
  select(!.type) |>
  group_by(id, Middle, Bottom) |>
  mutate(subset_better = c(0, diff(MASE))) |>
  ungroup() |>
  mutate(better = ifelse(subset_better < 0, TRUE, FALSE)) |>
  mutate(worse = ifelse(subset_better > 0, TRUE, FALSE)) |>
  select(better, worse) |>
  summarise(better = mean(better), worse = mean(worse))

sim_accuracy |> 
  filter(.model %in% c("MinTs", "MinTs_subset")) |>
  arrange(Middle, Bottom) |>
  select(!.type) |>
  group_by(id, Middle, Bottom) |>
  mutate(subset_better = c(0, diff(MASE))) |>
  ungroup() |>
  mutate(better = ifelse(subset_better < 0, TRUE, FALSE)) |>
  mutate(worse = ifelse(subset_better > 0, TRUE, FALSE)) |>
  select(better, worse) |>
  summarise(better = mean(better), worse = mean(worse))

