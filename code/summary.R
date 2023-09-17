# Load packages
library(dplyr)
library(reshape)
library(ggplot2)
library(tidyverse)
library(knitr)
library(kableExtra)
library(latex2exp)

# Functions
## MCB test
source("R/nemenyi.R")

## Other functions used for analysis
source("R/analysis.R")

#----------------------------------------------------------------------
# Table: Out-of-sample forecast performance (average RMSE/MASE)
#----------------------------------------------------------------------
measure <- "rmse"
data_label <- "tourism"
scenario <- NULL
horizons <- c(1, 4, 8, 12)
methods <- c("subset", "intuitive", "lasso")

out_all <- combine_table(data_label, methods, measure, scenario, horizons)
latex_table(out_all)

#----------------------------------------------------------------------
# Table: Number of time series retained in the structure with different methods for the tourism data
#----------------------------------------------------------------------
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
          Total = length(Total), "$\\lambda_0$" = 0, "$\\lambda_1$" = 0, "$\\lambda_2$" = 0)
tourism_subset_info <- apply(tourism_subset_info, 1, function(lentry){
  c(Top = sum(lentry[Top]), State = sum(lentry[State]), 
    Zone = sum(lentry[Zone]), Region = sum(lentry[Region]), 
    Total = sum(lentry[Total]), 
    "$\\lambda_0$" = round(tail(lentry, 3), 2)[1],
    "$\\lambda_1$" = round(tail(lentry, 3), 2)[2],
    "$\\lambda_2$" = round(tail(lentry, 3), 2)[3])
}) %>% t() %>% rbind(None, .)
rownames(tourism_subset_info) <- sub("_", "-", rownames(tourism_subset_info))

options(knitr.kable.NA = '')
tourism_subset_info |>
  kable(format = "latex",
        booktabs = TRUE,
        digits = 2,
        align = rep("r", 7),
        escape = FALSE,
        linesep = "") |>
  kable_styling(latex_options = c("hold_position", "repeat_header", "scale_down")) |>
  add_header_above(c("", "Number of time series retained" = 5, "Optimal parameters" = 2), align = "c")

#----------------------------------------------------------------------
# Plots: MCB test conducted on the methods examined using the tourism data, the ranks are computed considering all the time series of the hierarchy included in the data set, i.e. a total of 111 series in the hierarchy.
#----------------------------------------------------------------------
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

# Plot 1: h=1-12 average RMSE (111 series)
RMSE <- sapply(methods, function(lmethod){
  assign(lmethod, calc_rmse(fc = get(tolower(lmethod)), test = test, h = h))
}, USE.NAMES = TRUE) |> as.data.frame() |> data.matrix()
colnames(RMSE) <- sub("_", "-", colnames(RMSE))
# nemenyi(RMSE, conf.level = 0.95, plottype = "vmcb", Title = "RMSE - Tourism data") # select = match("WLSs-subset", colnames(RMSE))

highlight <- c("OLS", "WLSs", "WLSv", "MinTs")
target <- sapply(highlight, function(len){
  c(paste0(len, c("", "-subset")))
}) |> as.vector()
target <- c("Base", "BU", target, "Elasso") |> rev()
par(mfrow=c(1, 1),
    mar=c(5,0.1,0.1,0.3)
)
nemenyi(RMSE[,target], conf.level = 0.95, plottype = "vmcb",
        sort = FALSE,
        shadow = FALSE,
        group = list(1, 2:3, 4:5, 6:7, 8:9, 10:11),
        Title = "",
        Xlab = "Mean ranks",
        Ylab = "")

# Plot 2: Absolute error across all horizons (111*12 instances)
abserr_ts_horizon <- sapply(methods, function(lmethod){
  assign(lmethod, calc_abserror(fc = get(tolower(lmethod)), test = test))
}, USE.NAMES = TRUE)
nemenyi(abserr_ts_horizon, conf.level = 0.95, plottype = "vmcb", Title = "Absolute error across all horizons - Tourism data")

#----------------------------------------------------------------------
# Plot: Average reconciliation errors in terms of RMSE (1- to 12-step-ahead) after forecast reconciliation, for a single series, between disaggregate and aggregate views of the tourism data, for the different reconciliation methods. Time series are ordered in the horizontal axis.
#----------------------------------------------------------------------
RMSE_heatmap <- RMSE
rownames(RMSE_heatmap) <- 1:NROW(RMSE_heatmap)
RMSE_melt <- melt(RMSE_heatmap)  
colnames(RMSE_melt) <- c("Series", "Method", "RMSE")

# # Purples palette
# ggplot(RMSE_melt, aes(x = Series, y = Method, fill = RMSE)) +
#   geom_tile() +
#   scale_fill_gradientn(colors = rev(hcl.colors(50, "Purples"))[c(seq(1, 20, 10), seq(21, 50, 1))]) +
#   labs(x = "Time series", y = "") +
#   scale_x_discrete(expand = c(0, 0),
#                    breaks = seq(20, 100, 20)) +
#   scale_y_discrete(limits=rev(c("Base","BU","OLS","OLS-subset", "WLSs", "WLSs-subset", "WLSv", "WLSv-subset", "MinTs", "MinTs-subset"))) +
#   theme(
#     legend.text = element_text(face = "bold"),
#     plot.background = element_blank(),
#     panel.border = element_blank()
#   )

# 10*4
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
                                  "MinTs", "MinTs-subset", "Elasso"))) +
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
# Table: Number of time series retained in the structure with different methods for the simulation data
#----------------------------------------------------------------------
data_label <- "simulation"
methods <- c("subset", "intuitive", "lasso")
scenarios <- c("s0", "s1", "s2", "s3")
series_name <- c("Top", "A", "B", "AA", "AB", "BA", "BB")
simulation_info <- combine_z(data_label, methods, scenarios, series_name)
latex_sim_nos_table(simulation_info$out_s0$z,
                    simulation_info$out_s0$n)
latex_sim_nos_table(simulation_info$out_s1$z,
                    simulation_info$out_s1$n)
latex_sim_nos_table(simulation_info$out_s2$z,
                    simulation_info$out_s2$n)
latex_sim_nos_table(simulation_info$out_s3$z,
                    simulation_info$out_s3$n)

#----------------------------------------------------------------------
# Figures: MCB test for the simulation data
#----------------------------------------------------------------------
data_label <- "simulation"
methods <- c("subset", "intuitive", "lasso")
scenarios <- c("s0", "s1", "s2", "s3")
highlight <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")
target <- sapply(highlight, function(len){
  c(paste0(len, c("", "-subset", "-intuitive", "-lasso")))
}) |> as.vector()
target <- c("Base", "BU", target, "Elasso") |> rev()
h <- 16

RMSE_sim_s0 <- RMSE_MCB_sim("simulation", methods=methods, scenario="s0", h=16)[,target]
RMSE_sim_s1 <- RMSE_MCB_sim("simulation", methods=methods, scenario="s1", h=16)[,target]
RMSE_sim_s2 <- RMSE_MCB_sim("simulation", methods=methods, scenario="s2", h=16)[,target]
RMSE_sim_s3 <- RMSE_MCB_sim("simulation", methods=methods, scenario="s3", h=16)[,target]

MASE_sim_s0 <- MASE_MCB_sim("simulation", methods=methods, scenario="s0", h=16)[,target]
MASE_sim_s1 <- MASE_MCB_sim("simulation", methods=methods, scenario="s1", h=16)[,target]
MASE_sim_s2 <- MASE_MCB_sim("simulation", methods=methods, scenario="s2", h=16)[,target]
MASE_sim_s3 <- MASE_MCB_sim("simulation", methods=methods, scenario="s3", h=16)[,target]

par(mfrow=c(1, length(scenarios)),
    mar=c(2,0.1,0.1,0.3)
    #mar=c(5.1,2.1,4.1,2.1)
)

measure <- "RMSE"
scenarios <- c("s0", "s1", "s2", "s3")
for(j in scenarios){
  char <- paste0(measure, "_sim_", j)
  x <- get(char)
  
  if (j == "s0"){
    title <- "Simulation"
  } else if (j == "s1"){
    title <- "Scenario I"
  } else if (j == "s2"){
    title <- "Scenario II"
  } else if (j == "s3"){
    title <- "Scenario III"
  }
  
  nemenyi(x, conf.level = 0.95, plottype = "vmcb",
          sort = FALSE, 
          shadow = FALSE,
          group = list(1, 2:5, 6:9, 10:13, 14:17, 18:21, 22:23),
          Title = title,
          Xlab = "Mean ranks",
          Ylab = "")
}


