library(dplyr)
library(reshape)
library(ggplot2)
library(tidyverse)
library(knitr)
library(kableExtra)
library(latex2exp)

## MCB test
source("R/nemenyi.R")

## Other functions used for analysis
source("R/analysis.R")

#----------------------------------------------------------------------
# Simulation setup 1
# Exploring the effect of model misspecification
#----------------------------------------------------------------------
# RMSE table
measure <- "rmse"
data_label <- "simulation"
horizons <- c(1, 4, 8, 16)
methods <- c("subset", "intuitive", "lasso")
for (i in 1:3){
  scenario <- paste0("s", i)
  out_all <- combine_table(data_label, methods, measure, scenario = scenario, horizons)
  saveRDS(out_all, file = paste0("paper/results/sim_rmse_", scenario, ".rds"))
}

# Selection ratio table
scenarios <- c("s1", "s2", "s3")
series_name <- c("Top", "A", "B", "AA", "AB", "BA", "BB")
simulation_info <- combine_z(data_label, methods, scenarios, series_name)
saveRDS(simulation_info, file = "paper/results/sim_selection.rds")



