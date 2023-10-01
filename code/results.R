library(dplyr)
library(reshape)
library(ggplot2)
library(fabletools)
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

#----------------------------------------------------------------------
# Tourism application
#----------------------------------------------------------------------
tourism_data <- readRDS("data/tourism_data.rds")

# Time series plot
tourism_ts <- tourism_data |>
  select(-1) |>
  as_tibble() |>
  pivot_longer(cols = Total:GBD, names_to = "Series", values_to = "Value") |>
  mutate(Time = yearmonth(Time)) |>
  as_tsibble(index = Time, key = Series)
tourism_ts <- tourism_ts |>
  filter(Series %in% c("Total", LETTERS[1:7])) |>
  mutate(Series = recode(Series,
                        `A` = "NSW",
                        `B` = "VIC",
                        `C` = "QLD",
                        `D` = "SA",
                        `E` = "WA",
                        `F` = "TAS",
                        `G` = "NT"
  )) |>
  mutate(Series = factor(Series, 
                         levels = c("Total", "NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT")))
saveRDS(tourism_ts, "paper/results/tourism_hts.rds")

# RMSE table
measure <- "rmse"
data_label <- "tourism"
scenario <- NULL
horizons <- c(1, 4, 8, 12)
methods <- c("subset", "intuitive", "lasso")

out_all <- combine_table(data_label, methods, measure, scenario, horizons)
saveRDS(out_all, file = paste0("paper/results/tourism_rmse.rds"))

# Series retained table
tourism_subset_reconsf <- readRDS(file = "data_new/tourism_subset_reconsf.rds")
tourism_lasso_reconsf <- readRDS(file = "data_new/tourism_lasso_reconsf.rds")
Top <- 1
State <- 2:8
Zone <- 9:35
Region <- 36:111
Total <- 1:111

tourism_subsetonly_reconsf <- c(tourism_subset_reconsf[[1]][grepl("_subset", names(tourism_subset_reconsf[[1]]))],
                                Elasso = list(tourism_lasso_reconsf[[1]]$Elasso))
tourism_subset_z <- sapply(tourism_subsetonly_reconsf, function(lentry) lentry$z) |> t()
tourism_subset_lambda <- sapply(tourism_subsetonly_reconsf, function(lentry) lentry$lambda_opt[c("l0", "l1", "l2")] |> unname()) |> t()
tourism_subset_info <- data.frame(tourism_subset_z, tourism_subset_lambda)

None <- c(Top = length(Top), State = length(State), Zone = length(Zone), Region = length(Region),
          Total = length(Total), "$\\lambda$" = NA, "$\\lambda_0$" = NA, "$\\lambda_2$" = NA)
tourism_subset_info <- apply(tourism_subset_info, 1, function(lentry){
  c(Top = sum(lentry[Top]), State = sum(lentry[State]), 
    Zone = sum(lentry[Zone]), Region = sum(lentry[Region]), 
    Total = sum(lentry[Total]), 
    "$\\lambda$" = round(tail(lentry, 3), 2)[2],
    "$\\lambda_0$" = round(tail(lentry, 3), 2)[1],
    "$\\lambda_2$" = round(tail(lentry, 3), 2)[3])
}) %>% t() %>% rbind(None, .)
rownames(tourism_subset_info) <- sub("_", "-", rownames(tourism_subset_info))
saveRDS(tourism_subset_info, "paper/results/tourism_info.rds")

# Heatmap
data_label <- "tourism"
method_label <- "subset"
h <- 12

test <- readRDS(file = paste0("data/", data_label, "_test.rds"))
reconsf <- readRDS(file = paste0("data_new/", data_label, "_", method_label, "_reconsf.rds"))
lasso_reconsf <- readRDS(file = "data_new/tourism_lasso_reconsf.rds")
reconsf[[1]] <- c(reconsf[[1]], Elasso = list(lasso_reconsf[[1]]$Elasso))

methods <- names(reconsf[[1]])
for(method in methods) {
  out <- unique(test$Index) |> 
    purrr::map(\(index) 
               extract_element(data = reconsf, index = index, 
                               method = method, element = "y_tilde")) %>% 
    do.call(rbind, .)
  assign(tolower(method), out)
}

RMSE <- sapply(methods, function(lmethod){
  assign(lmethod, calc_rmse(fc = get(tolower(lmethod)), test = test, h = h))
}, USE.NAMES = TRUE) |> as.data.frame() |> data.matrix() # h=1-12 average RMSE (111 series)
colnames(RMSE) <- sub("_", "-", colnames(RMSE))

RMSE_heatmap <- RMSE
rownames(RMSE_heatmap) <- 1:NROW(RMSE_heatmap)
RMSE_melt <- melt(RMSE_heatmap)  
colnames(RMSE_melt) <- c("Series", "Method", "RMSE")
saveRDS(RMSE_melt, file = "paper/results/tourism_heatmap.rds")

ggplot(RMSE_melt, aes(x = Series, y = Method, fill = RMSE)) +
  geom_tile() +
  geom_vline(xintercept = c(1.5, 8.5, 35.5), linetype = "dashed", linewidth = 0.5) +
  scale_fill_gradientn(colors = c("#f7d9a6", 
                                  rev(hcl.colors(100, "Purples"))[c(seq(1, 50, 10), seq(51, 100, 1))])) +
  labs(x = "Time series", y = "") +
  scale_x_continuous(expand = c(0, 0), 
                     breaks =  seq(20, 100, 20),
                     sec.axis = dup_axis(name = "",
                                         breaks = c(4.5, 11.5, 39),
                                         labels = c("States", "Zones", "Regions"))) +
  scale_y_discrete(expand = c(0, 0), 
                   limits = rev(c("Base", "BU", "OLS", "OLS-subset", 
                                  "WLSs", "WLSs-subset", "WLSv", "WLSv-subset", 
                                  "MinTs", "MinTs-subset", "EMinT", "Elasso"))) +
  theme(
    plot.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    axis.title.x = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    axis.text = element_text(face = "bold", size = 10),
    axis.ticks.x.top = element_blank()
  ) +
  guides(fill = guide_colourbar(barwidth = 10,
                                barheight = 1.5))

#----------------------------------------------------------------------
# Labour application
#----------------------------------------------------------------------
labour_data <- readRDS("data/labour_data.rds")

# Time series plot
labour_ts <- labour_data |>
  select(-1) |>
  as_tibble() |>
  pivot_longer(cols = Total:D6WAS, names_to = "Series", values_to = "Value") |>
  mutate(Time = yearmonth(Time))
labour_ts <- labour_ts |>
  filter(grepl("Duration", Series)|grepl("Total", Series)|grepl("STT", Series)) |>
  mutate(Series = recode(Series,
                         `Duration/D1` = "Under 1 month",
                         `Duration/D2` = "1-3 months",
                         `Duration/D3` = "3-6 months",
                         `Duration/D4` = "6-12 months",
                         `Duration/D5` = "1-2 years",
                         `Duration/D6` = "2 years and over",
                         `STT/NSW` = "NSW",
                         `STT/VIC` = "VIC",
                         `STT/QLD` = "QLD",
                         `STT/SAS` = "SA",
                         `STT/WAS` = "WA",
                         `STT/TAS` = "TAS",
                         `STT/NTT` = "NT",
                         `STT/ACT` = "ACT"
  )) |>
  mutate(Series = factor(Series, 
                         levels = c("Total", 
                                    "NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT",
                                    "Under 1 month", "1-3 months", "3-6 months", "6-12 months", 
                                    "1-2 years", "2 years and over"))) |>
  as_tsibble(index = Time, key = Series)

saveRDS(labour_ts, "paper/results/labour_gts.rds")

# RMSE table
measure <- "rmse"
data_label <- "labour"
scenario <- NULL
horizons <- c(1, 4, 8, 12)
methods <- c("subset", "intuitive", "lasso")

out_all <- combine_table(data_label, methods, measure, scenario, horizons)
saveRDS(out_all, file = paste0("paper/results/labour_rmse.rds"))

# Series retained table
labour_subset_reconsf <- readRDS(file = "data_new/labour_subset_reconsf.rds")
labour_lasso_reconsf <- readRDS(file = "data_new/labour_lasso_reconsf.rds")
Top <- 1
Duration <- 2:7
STT <- 8:15
Duration_STT <- 16:63
Total <- 1:63

labour_subsetonly_reconsf <- c(labour_subset_reconsf[[1]][grepl("_subset", names(labour_subset_reconsf[[1]]))],
                               Elasso = list(labour_lasso_reconsf[[1]]$Elasso))
labour_subset_z <- sapply(labour_subsetonly_reconsf, function(lentry) lentry$z) |> t()
labour_subset_lambda <- sapply(labour_subsetonly_reconsf, function(lentry) lentry$lambda_opt[c("l0", "l1", "l2")] |> unname()) |> t()
labour_subset_info <- data.frame(labour_subset_z, labour_subset_lambda)

None <- c(Top = length(Top), Duration = length(Duration), STT = length(STT), Duration_STT = length(Duration_STT),
          Total = length(Total), "$\\lambda$" = NA, "$\\lambda_0$" = NA, "$\\lambda_2$" = NA)
labour_subset_info <- apply(labour_subset_info, 1, function(lentry){
  c(Top = sum(lentry[Top]), Duration = sum(lentry[Duration]), 
    STT = sum(lentry[STT]), Duration_STT = sum(lentry[Duration_STT]), 
    Total = sum(lentry[Total]), 
    "$\\lambda$" = round(tail(lentry, 3), 2)[2],
    "$\\lambda_0$" = round(tail(lentry, 3), 2)[1],
    "$\\lambda_2$" = round(tail(lentry, 3), 2)[3])
}) %>% t() %>% rbind(None, .)
rownames(labour_subset_info) <- sub("_", "-", rownames(labour_subset_info))
colnames(labour_subset_info) <- sub("Duration_STT", "Duration x STT", colnames(labour_subset_info))
saveRDS(labour_subset_info, "paper/results/labour_info.rds")

