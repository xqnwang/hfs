library(tidyverse)
library(magrittr)
library(future)
library(forecast)
library(latex2exp)

# Setup
input <- commandArgs(trailingOnly = TRUE)

method_label <- input[1] # "subset", "lasso", "intuitive"
data_label <- input[2] # "simulation", "tourism"
if(is.na(input[3])){
  scenario <- NULL
} else{
  scenario <- input[3]
} # NULL, "s1", "s2", "s3"

# Utility functions
source("R/analysis.R")

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
  # Import results
  freq <- 4
  if (is.null(scenario)){
    reconsf <- readRDS(file = paste0("data_new/", data_label, "_", method_label, "_reconsf.rds"))
  } else{
    reconsf <- readRDS(file = paste0("data_new/", data_label, "_", method_label, "_reconsf_", scenario, ".rds"))
  }
  train <- readRDS(file = paste0("data/", data_label, "_train.rds"))
  test <- readRDS(file = paste0("data/", data_label, "_test.rds"))
  method <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")
  
  # Structure information used to calculate RMSE across levels
  top <- 1
  middle <- 2:3
  bottom <- 4:7
  avg <- 1:7
  horizon <- c(1, 4, 8, 16)
  
  # Reconciliation methods considered
  methods <- c("Base", "BU", 
               sapply(method, function(l) c(l, paste0(l, "_", method_label))) |> as.character())
  if (method_label == "lasso") methods <- c(methods, "Elasso")
  reconcile_methods <- grep(method_label, methods, value = TRUE)
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
  # Import results
  freq <- 12
  scenario <- NULL
  reconsf <- readRDS(file = paste0("data_new/", data_label, "_", method_label, "_reconsf.rds"))
  train <- readRDS(file = paste0("data/", data_label, "_train.rds"))
  test <- readRDS(file = paste0("data/", data_label, "_test.rds"))
  method <- c("OLS", "WLSs", "WLSv", "MinTs")
  
  # Structure information used to calculate RMSE across levels
  top <- 1
  state <- 2:8
  zone <- 9:35
  region <- 36:111
  avg <- 1:111
  horizon <- c(1, 4, 8, 12)
  
  # Reconciliation methods considered
  methods <- c("Base", "BU", 
               sapply(method, function(l) c(l, paste0(l, "_", method_label))) |> as.character())
  if (method_label == "lasso") methods <- c(methods, "Elasso")
  reconcile_methods <- grep(method_label, methods, value = TRUE)
}

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
if (data_label == "prison"){
  # Import results
  freq <- 4
  scenario <- NULL
  reconsf <- readRDS(file = paste0("data_new/", data_label, "_", method_label, "_reconsf.rds"))
  train <- readRDS(file = paste0("data/", data_label, "_train.rds"))
  test <- readRDS(file = paste0("data/", data_label, "_test.rds"))
  method <- c("OLS", "WLSs", "WLSv", "MinTs")
  
  # Structure information used to calculate RMSE across levels
  top <- 1
  gender <- 2:3
  legal <- 4:5
  state <- 6:13
  gender_legal <- 14:17
  gender_state <- 18:33
  legal_state <- 34:49
  gender_legal_state <- 50:81
  avg <- 1:81
  horizon <- c(1, 2, 4, 8)
  
  # Reconciliation methods considered
  methods <- c("Base", "BU", 
               sapply(method, function(l) c(l, paste0(l, "_", method_label))) |> as.character())
  if (method_label == "lasso") methods <- c(methods, "Elasso")
  reconcile_methods <- grep(method_label, methods, value = TRUE)
}

#################################################
# Extract reconciled forecasts
#################################################
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
# Calculate RMSE & MASE values
#################################################
for(h in horizon){
  rmse <- lapply(methods, function(lmethod){
    assign(lmethod, calc_rmse(fc = get(tolower(lmethod)), test = test, h = h))
  })
  names(rmse) <- methods
  
  if (data_label == "simulation"){
    out <- bind_rows(rmse, .id = "Method") |>
      rowwise() |>
      mutate(Top = mean(c_across(top + 1)),
             Middle = mean(c_across((middle + 1))),
             Bottom = mean(c_across((bottom + 1))),
             Average = mean(c_across(avg + 1))) |>
      select(Method, Top, Middle, Bottom, Average)
  } else if (data_label == "tourism"){
    out <- bind_rows(rmse, .id = "Method") |>
      rowwise() |>
      mutate(Top = mean(c_across(top + 1)),
             State = mean(c_across(state + 1)),
             Zone = mean(c_across(zone + 1)),
             Region = mean(c_across(region + 1)),
             Average = mean(c_across(avg + 1))) |>
      select(Method, Top, State, Zone, Region, Average)
  } else if (data_label == "prison"){
    out <- bind_rows(rmse, .id = "Method") |>
      rowwise() |>
      mutate(Top = mean(c_across(top + 1)),
             Gender = mean(c_across(gender + 1)),
             Legal = mean(c_across(legal + 1)),
             State = mean(c_across(state + 1)),
             Gender_Legal = mean(c_across(gender_legal + 1)),
             Gender_State = mean(c_across(gender_state + 1)),
             Legal_State = mean(c_across(legal_state + 1)),
             Gender_Legal_State = mean(c_across(gender_legal_state + 1)),
             Average = mean(c_across(avg + 1))) |>
      select(Method, Top, Gender, Legal, State, Gender_Legal, Gender_State, Legal_State, Gender_Legal_State, Average)
  }
  assign(paste0("rmse_h", h), out)
  
  if (is.null(scenario)){
    saveRDS(out, file = paste0("data_new/", data_label, "_", method_label, "_reconsf_rmse_", h, ".rds"))
  } else{
    saveRDS(out, file = paste0("data_new/", data_label, "_", method_label, "_reconsf_", scenario, "_rmse_", h, ".rds"))
  }
}

for(h in horizon){
  mase <- lapply(methods, function(lmethod){
    assign(lmethod, calc_mase(fc = get(tolower(lmethod)), 
                              train = train, test = test, 
                              freq = freq, h = h))
  })
  names(mase) <- methods
  
  if (data_label == "simulation"){
    out <- bind_rows(mase, .id = "Method") |>
      rowwise() |>
      mutate(Top = mean(c_across(top + 1)),
             Middle = mean(c_across((middle + 1))),
             Bottom = mean(c_across((bottom + 1))),
             Average = mean(c_across(avg + 1))) |>
      select(Method, Top, Middle, Bottom, Average)
  } else if (data_label == "tourism"){
    out <- bind_rows(mase, .id = "Method") |>
      rowwise() |>
      mutate(Top = mean(c_across(top + 1)),
             State = mean(c_across(state + 1)),
             Zone = mean(c_across(zone + 1)),
             Region = mean(c_across(region + 1)),
             Average = mean(c_across(avg + 1))) |>
      select(Method, Top, State, Zone, Region, Average)
  } else if (data_label == "prison"){
    out <- bind_rows(mase, .id = "Method") |>
      rowwise() |>
      mutate(Top = mean(c_across(top + 1)),
             Gender = mean(c_across(gender + 1)),
             Legal = mean(c_across(legal + 1)),
             State = mean(c_across(state + 1)),
             Gender_Legal = mean(c_across(gender_legal + 1)),
             Gender_State = mean(c_across(gender_state + 1)),
             Legal_State = mean(c_across(legal_state + 1)),
             Gender_Legal_State = mean(c_across(gender_legal_state + 1)),
             Average = mean(c_across(avg + 1))) |>
      select(Method, Top, Gender, Legal, State, Gender_Legal, Gender_State, Legal_State, Gender_Legal_State, Average)
  }
  assign(paste0("mase_h", h), out)
  
  if (is.null(scenario)){
    saveRDS(out, file = paste0("data_new/", data_label, "_", method_label, "_reconsf_mase_", h, ".rds"))
  } else{
    saveRDS(out, file = paste0("data_new/", data_label, "_", method_label, "_reconsf_", scenario, "_mase_", h, ".rds"))
  }
}

#################################################
# Extract z
#################################################
z_summary <- NULL
for(method in reconcile_methods) { 
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
  saveRDS(z_summary, file = paste0("data_new/", data_label, "_", method_label, "_reconsf_z_summary.rds"))
} else{
  saveRDS(z_summary, file = paste0("data_new/", data_label, "_", method_label, "_reconsf_", scenario, "_z_summary.rds"))
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


