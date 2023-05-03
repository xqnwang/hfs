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
    summarise_at(1:NCOL(err), mean)
  return(rmse)
}

#################################################
# Evaluation
#################################################
data_label <- "simulation"
reconsf <- readRDS(file = paste0("data/", data_label, "_reconsf.rds"))
# scenario <- "s2"
# reconsf <- readRDS(file = paste0("data/", data_label, "_reconsf_", scenario, ".rds"))
test <- readRDS(file = paste0("data/", data_label, "_test.rds"))

methods <- c("Base", "BU", "OLS", "OLS_subset", 
             "WLSs", "WLSs_subset", "WLSv", "WLSv_subset", 
             "MinT", "MinT_subset", 
             "MinTs", "MinTs_subset")
subset_methods <- c("OLS_subset", "WLSs_subset", "WLSv_subset", "MinT_subset", "MinTs_subset")
indices <- unique(test$Index)

#################################################
# Extract reconciled forecasts
#################################################
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
  z_summary <- rbind(z_summary, out)
}
series_name <- c(colnames(test), "Method")
colnames(z_summary) <- series_name

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
       y= "")

#################################################
# Extract lambda0_report
#################################################
lambda0_summary <- NULL
for(method in subset_methods) { 
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
  labs(title = TeX(r"(Frequency of being selected as the optimal $\lambda_0$)"),
       x = TeX(r"(index of $\lambda_0$, $\lambda_{0\min} = 0$)"),
       y= "")
