library(dplyr)
library(reshape)
library(ggplot2)
library(fabletools)
library(tidyverse)
library(knitr)
library(kableExtra)
library(latex2exp)
library(lubridate)
library(fable)
library(tsibble)
library(patchwork)

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
scenarios <- c("s1", "s2", "s3")
horizons <- c(1, 4, 8, 16)
methods <- c("subset", "intuitive", "lasso")
sim_rmse <- lapply(scenarios, function(scenario){
  out_all <- combine_table(data_label, methods, measure, scenario = scenario, horizons)
  out_all$table_out$Method <- gsub("intuitive", "parsim", out_all$table_out$Method)
  return(out_all)
})
names(sim_rmse) <- scenarios
saveRDS(sim_rmse, file = "paper/results/sim_rmse.rds")

# MCB tests for each scenario
h <- 16
highlight <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")
target <- sapply(highlight, function(len){
  c(paste0(len, c("", "-subset", "-intuitive", "-lasso")))
}) |> as.vector()
target <- c("Base", "BU", target, "EMinT", "Elasso") |> rev()
sim_rmse_mcb <- lapply(scenarios, function(scenario){
  for (method in methods){
    readRDS(paste0("data_new/", data_label, "_", method, "_reconsf_", scenario, "_rmse_hts_", h, ".rds")) |>
      assign(method, value = _)
  }
  rmse_hfs <- cbind(subset,
                    intuitive[, grepl("intuitive", colnames(intuitive))],
                    lasso[, grepl("lasso", colnames(lasso))])
  colnames(rmse_hfs) <- gsub("_", "-", colnames(rmse_hfs))
  colnames(rmse_hfs) <- gsub("intuitive", "parsim", colnames(rmse_hfs))
  rmse_hfs <- rmse_hfs[, gsub("intuitive", "parsim", target)]
  return(rmse_hfs)
})
names(sim_rmse_mcb) <- scenarios
saveRDS(sim_rmse_mcb, file = "paper/results/sim_rmse_mcb.rds")

# Selection ratio table
series_name <- c("Top", "A", "B", "AA", "AB", "BA", "BB")
simulation_info <- combine_z(data_label, methods, scenarios, series_name)
simulation_info <- lapply(simulation_info, function(len){
  len$n$Method <- gsub("intuitive", "parsim", len$n$Method)
  rownames(len$z) <- gsub("intuitive", "parsim", rownames(len$z))
  return(len)
})
names(simulation_info) <- scenarios
saveRDS(simulation_info, file = "paper/results/sim_selection.rds")

# Forecast variance (RMSSE)
freq <- 4; h <- 16
train <- readRDS(file = paste0("data/", data_label, "_train.rds"))
test <- readRDS(file = paste0("data/", data_label, "_test.rds"))
indices <- unique(test$Index)
simulation_series_rmsse <- lapply(scenarios, function(scenario) {
  reconsf <- readRDS(file = paste0("data_new/", data_label, "_subset_reconsf_", scenario, ".rds"))
  fc <- indices |>
    purrr::map(\(index)
               extract_element(data = reconsf, index = index,
                               method = "Base", element = "y_tilde")) %>% 
    do.call(rbind, .)
  return(calc_rmsse(fc, train, test, freq, h))
})
names(simulation_series_rmsse) <- scenarios
saveRDS(simulation_series_rmsse, file = "paper/results/sim_series_rmsse.rds")

# Forecast correlation
h <- 1
train <- readRDS(file = paste0("data/", data_label, "_train.rds"))
test <- readRDS(file = paste0("data/", data_label, "_test.rds"))
indices <- unique(test$Index)
simulation_cormat <- lapply(scenarios, function(scenario) {
  reconsf <- readRDS(file = paste0("data_new/", data_label, "_subset_reconsf_", scenario, ".rds"))
  fc <- indices |>
    purrr::map(\(index)
               extract_element(data = reconsf, index = index,
                               method = "Base", element = "y_tilde")) %>% 
    do.call(rbind, .)
  er_h <- calc_error(fc, test, h)
  cormat <- cor(er_h)
  cormat[upper.tri(cormat)]<- NA
  melted_cormat <- reshape2::melt(cormat, na.rm = TRUE)
  return(melted_cormat)
})
names(simulation_cormat) <- scenarios
saveRDS(simulation_cormat, file = "paper/results/sim_cormat.rds")

#----------------------------------------------------------------------
# Simulation setup 2
# Exploring the effect of correlation
#----------------------------------------------------------------------
# Time plot
dat <- readr::read_csv("data/corr_1_data.csv")
data <- dat |>
  filter(Index == 1) |>
  select(!Index) |>
  mutate(A = AA+AB,
         B = BA+BB,
         Total = A+B) |>
  pivot_longer(cols = !Time, names_to = "Series", values_to = "Value") |>
  as_tsibble(index = Time, key = Series) |>
  mutate(Series = factor(Series, 
                         levels = c("Total", "A", "B", "AA", "AB", "BA", "BB")))
fc <- readRDS("data/corr_1_basefc.rds")
fit <- readRDS("data/corr_1_fits.rds")
resid <- readRDS("data/corr_1_resids.rds")
fc <- fc[fc$Index == 1, ]
fit <- fit[fit$Index == 1, ]
resid <- resid[resid$Index == 1, ]
resid <- resid |>
  as_tibble() |>
  select(!Index) |>
  mutate(Time = 1:100) |>
  pivot_longer(cols = !Time, names_to = "Series", values_to = "Value") |> 
  as_tsibble(index = Time, key = Series) |>
  mutate(Series = factor(Series, 
                         levels = c("Total", "A", "B", "AA", "AB", "BA", "BB")))

saveRDS(data, file = "paper/results/corr_data_neg.rds")
saveRDS(resid, file = "paper/results/corr_resid_neg.rds")

# RMSE table
measure <- "rmse"
data_label <- "corr"
corr <- seq(-0.8, 0.8, 0.2)
index <- c(1, 3, 5, 7, 9)
methods <- c("subset", "intuitive", "lasso")

out_all <- combine_corr_table(data_label, methods, corr, index, measure)
out_all$table_out$Method <- gsub("intuitive", "parsim", out_all$table_out$Method)
saveRDS(out_all, file = "paper/results/corr_rmse.rds")

# Selection ratio table
series_name <- c("Top", "A", "B", "AA", "AB", "BA", "BB")
corr_selection <- lapply(index, function(i){
  corr_info <- combine_z(paste0("corr_", i), methods, scenarios = "s0", series_name)
  corr_info_out <- lapply(corr_info, function(len){
    len$n$Method <- gsub("intuitive", "parsim", len$n$Method)
    rownames(len$z) <- gsub("intuitive", "parsim", rownames(len$z))
    return(len)
  })
  return(corr_info_out)
})
names(corr_selection) <- paste0("p", index)
saveRDS(corr_selection, file = "paper/results/corr_selection.rds")

# Forecast variance (RMSSE)
freq <- 1; h <- 1
corr_series_rmsse <- lapply(index, function(i){
  train <- readRDS(file = paste0("data/", data_label, "_", i, "_train.rds"))
  test <- readRDS(file = paste0("data/", data_label, "_", i, "_test.rds"))
  indices <- unique(test$Index)
  reconsf <- readRDS(file = paste0("data_new/", data_label, "_", i, "_subset_reconsf.rds"))
  fc <- indices |>
    purrr::map(\(index)
               extract_element(data = reconsf, index = index,
                               method = "Base", element = "y_tilde")) %>% 
    do.call(rbind, .)
  list(out_s0 = calc_rmsse(fc, train, test, freq, h))
})
names(corr_series_rmsse) <- paste0("p", index)
saveRDS(corr_series_rmsse, file = "paper/results/corr_series_rmsse.rds")

# MCB tests for each error correlation
h <- 1
highlight <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")
target <- sapply(highlight, function(len){
  c(paste0(len, c("", "-subset", "-intuitive", "-lasso")))
}) |> as.vector()
target <- c("Base", "BU", target, "EMinT", "Elasso") |> rev()
corr_rmse_mcb <- lapply(index, function(i){
  for (method in methods){
    readRDS(paste0("data_new/", data_label, "_", i, "_", method, "_reconsf_rmse_hts_", h, ".rds")) |>
      assign(method, value = _)
  }
  rmse_hfs <- cbind(subset,
                    intuitive[, grepl("intuitive", colnames(intuitive))],
                    lasso[, grepl("lasso", colnames(lasso))])
  colnames(rmse_hfs) <- gsub("_", "-", colnames(rmse_hfs))
  rmse_hfs <- rmse_hfs[, target]
  return(rmse_hfs)
})
names(corr_rmse_mcb) <- paste0("p", index)
saveRDS(corr_rmse_mcb, file = paste0("paper/results/corr_rmse_mcb.rds"))
  
# Forecast correlation
data_label <- "corr"
h <- 1
corr_cormat <- lapply(index, function(len) {
  train <- readRDS(file = paste0("data/", data_label, "_", len, "_train.rds"))
  test <- readRDS(file = paste0("data/", data_label, "_", len, "_test.rds"))
  indices <- unique(test$Index)
  reconsf <- readRDS(file = paste0("data_new/", data_label, "_", len, "_subset_reconsf.rds"))
  fc <- indices |>
    purrr::map(\(x)
               extract_element(data = reconsf, index = x,
                               method = "Base", element = "y_tilde")) %>% 
    do.call(rbind, .)
  er_h <- calc_error(fc, test, h)
  cormat <- cor(er_h)
  cormat[upper.tri(cormat)]<- NA
  melted_cormat <- reshape2::melt(cormat, na.rm = TRUE)
  return(melted_cormat)
})
names(corr_cormat) <- paste0("p", index)
saveRDS(corr_cormat, file = "paper/results/corr_cormat.rds")

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

# RMSE table - average results for rolling windows
measure <- "rmse"
data_label <- "tourism"
scenario <- NULL
horizons <- c(1, 4, 8, 12)
methods <- c("subset", "intuitive", "lasso")
test_indices <- 1:12

for (method_label in methods){
  for (h in horizons){
    for (test_lable in test_indices){
      assign(paste0(data_label, "_", test_lable, "_", method_label, "_", measure, "_h", h), 
             readRDS(file = paste0("data_new/", data_label, "_", test_lable, "_", method_label, "_reconsf", 
                                   "_", measure, "_", h, ".rds"))
      )
    }
    Method <- get(paste0(data_label, "_", test_lable, "_", method_label, "_", measure, "_h", h))[,1]
    out <- mget(paste0(data_label, "_", test_indices, "_", method_label, "_", measure, "_h", h), inherits = TRUE) |> 
      lapply(function(lentry) lentry[, -1])
    out_h <- data.frame(Method = Method, Reduce('+', out)/length(out))
    saveRDS(out_h, file = paste0("data_new/", data_label, "_", method_label, "_reconsf_rmse_", h, ".rds"))
  }
}

out_all <- combine_table(data_label, methods, measure, scenario, horizons)
out_all$table_out$Method <- gsub("intuitive", "parsim", out_all$table_out$Method)
saveRDS(out_all, file = paste0("paper/results/tourism_rmse.rds"))

# RMSE table - the last window
measure <- "rmse"
data_label <- "tourism_1"
scenario <- NULL
horizons <- c(1, 4, 8, 12)
methods <- c("subset", "intuitive", "lasso")

out_1_all <- combine_table(data_label, methods, measure, scenario, horizons)
out_1_all$table_out$Method <- gsub("intuitive", "parsim", out_1_all$table_out$Method)
saveRDS(out_1_all, file = "paper/results/tourism_1_rmse.rds")

# Series retained table
tourism_subset_reconsf <- readRDS(file = "data_new/tourism_1_subset_reconsf.rds")
tourism_lasso_reconsf <- readRDS(file = "data_new/tourism_1_lasso_reconsf.rds")
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

# Forecast variance (RMSSE) - the last window
freq <- 12; h <- 12
data_label <- "tourism_1"
train <- readRDS(file = paste0("data/", data_label, "_train.rds"))
test <- readRDS(file = paste0("data/", data_label, "_test.rds"))
indices <- unique(test$Index)
reconsf <- readRDS(file = paste0("data_new/", data_label, "_subset_reconsf.rds"))
fc <- indices |>
  purrr::map(\(index)
             extract_element(data = reconsf, index = index,
                             method = "Base", element = "y_tilde")) %>% 
  do.call(rbind, .)
tourism_series_rmsse <- calc_rmsse(fc, train, test, freq, h) |> as.matrix()
tourism_select_index <- rbind(None = rep(1, ncol(tourism_subset_z)), tourism_subset_z)
tourism_select_index <- tourism_select_index == TRUE
tourism_hlevel <- list(Top, State, Zone, Region)
tourism_selected_rmsse <- sapply(1:nrow(tourism_select_index), function(i) {
  ind <- tourism_select_index[i, ]
  avg <- sapply(tourism_hlevel, function(j) {
    mean(tourism_series_rmsse[, j][ind[j]])
  })
  return(avg)
}) |> t()
tourism_selected_rmsse[is.na(tourism_selected_rmsse)] <- 0
colnames(tourism_selected_rmsse) <- c("Top", "State", "Zone", "Region")
rownames(tourism_selected_rmsse) <- rownames(tourism_select_index)
saveRDS(tourism_selected_rmsse, file = "paper/results/tourism_selected_rmsse.rds")

# Heatmap
data_label <- "tourism_1"
method_label <- "subset"
h <- 12

test <- readRDS(file = paste0("data/", data_label, "_test.rds"))
reconsf <- readRDS(file = paste0("data_new/", data_label, "_", method_label, "_reconsf.rds"))
lasso_reconsf <- readRDS(file = paste0("data_new/", data_label, "_lasso", "_reconsf.rds"))
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

# RMSE table - average results for rolling windows
measure <- "rmse"
data_label <- "labour"
scenario <- NULL
horizons <- c(1, 4, 8, 12)
methods <- c("subset", "intuitive", "lasso")
test_indices <- 1:12

for (method_label in methods){
  for (h in horizons){
    for (test_lable in test_indices){
      assign(paste0(data_label, "_", test_lable, "_", method_label, "_", measure, "_h", h), 
             readRDS(file = paste0("data_new/", data_label, "_", test_lable, "_", method_label, "_reconsf", 
                                   "_", measure, "_", h, ".rds"))
      )
    }
    Method <- get(paste0(data_label, "_", test_lable, "_", method_label, "_", measure, "_h", h))[,1]
    out <- mget(paste0(data_label, "_", test_indices, "_", method_label, "_", measure, "_h", h), inherits = TRUE) |> 
      lapply(function(lentry) lentry[, -1])
    out_h <- data.frame(Method = Method, Reduce('+', out)/length(out))
    saveRDS(out_h, file = paste0("data_new/", data_label, "_", method_label, "_reconsf_rmse_", h, ".rds"))
  }
}

out_all <- combine_table(data_label, methods, measure, scenario, horizons)
out_all$table_out$Method <- gsub("intuitive", "parsim", out_all$table_out$Method)
saveRDS(out_all, file = paste0("paper/results/labour_rmse.rds"))

# RMSE table - the last window
measure <- "rmse"
data_label <- "labour_1"
scenario <- NULL
horizons <- c(1, 4, 8, 12)
methods <- c("subset", "intuitive", "lasso")

out_1_all <- combine_table(data_label, methods, measure, scenario, horizons)
out_1_all$table_out$Method <- gsub("intuitive", "parsim", out_1_all$table_out$Method)
saveRDS(out_1_all, file = "paper/results/labour_1_rmse.rds")

# Series retained table
data_label <- "labour_1"
labour_subset_reconsf <- readRDS(file = paste0("data_new/", data_label, "_subset_reconsf.rds"))
labour_lasso_reconsf <- readRDS(file = paste0("data_new/", data_label, "_lasso_reconsf.rds"))
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

# Forecast variance (RMSSE) - the last window
freq <- 12; h <- 12
data_label <- "labour_1"
train <- readRDS(file = paste0("data/", data_label, "_train.rds"))
test <- readRDS(file = paste0("data/", data_label, "_test.rds"))
indices <- unique(test$Index)
reconsf <- readRDS(file = paste0("data_new/", data_label, "_subset_reconsf.rds"))
fc <- indices |>
  purrr::map(\(index)
             extract_element(data = reconsf, index = index,
                             method = "Base", element = "y_tilde")) %>% 
  do.call(rbind, .)
labour_series_rmsse <- calc_rmsse(fc, train, test, freq, h) |> as.matrix()
labour_select_index <- rbind(None = rep(1, ncol(labour_subset_z)), labour_subset_z)
labour_select_index <- labour_select_index == TRUE
labour_hlevel <- list(Top, Duration, STT, Duration_STT)
labour_selected_rmsse <- sapply(1:nrow(labour_select_index), function(i) {
  ind <- labour_select_index[i, ]
  avg <- sapply(labour_hlevel, function(j) {
    mean(labour_series_rmsse[, j][ind[j]])
  })
  return(avg)
}) |> t()
labour_selected_rmsse[is.na(labour_selected_rmsse)] <- 0
colnames(labour_selected_rmsse) <- c("Top", "Duration", "STT", "Duration_STT")
rownames(labour_selected_rmsse) <- rownames(labour_select_index)
saveRDS(labour_selected_rmsse, file = "paper/results/labour_selected_rmsse.rds")
