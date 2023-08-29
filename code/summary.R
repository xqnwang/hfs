# Load packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(knitr)
library(kableExtra)
library(latex2exp)

# Functions
## MCB test
source("R/nemenyi.R")

## Combine results of different methods
combine_table <- function(data_label, methods, measure, scenario = NULL, horizons){
  for(method_label in methods){
    data_method <- paste0(data_label, "_", method_label)
    for(h in horizons){
      assign(paste0(data_method, "_", measure, "_h", h), 
             readRDS(file = paste0("data_new/", data_method, "_reconsf", 
                                   ifelse(is.null(scenario), "", paste0("_", scenario)), "_", measure, "_", h, ".rds"))
      )
    }
    example <- get(paste0(data_method, "_", measure, "_h", horizons[1]))
    levels <- colnames(example)[-1]
    method <- sub("_", "-", example$Method)
    out <- lapply(levels, function(level){
      mget(paste0(data_method, "_", measure, "_h", horizons), inherits = TRUE) |> 
        sapply(function(lentry) as.numeric(lentry[[level]]))
    })
    out <- data.frame(method, do.call(cbind, out))
    colnames(out) <- c("Method", c(ifelse(horizons == 1, paste0("h=", 1), paste0("1-", horizons))) |> rep(length(levels)))
    assign(paste0(data_method, "_", measure), out)
  }
  
  out_com <- mget(paste0(data_label, "_", methods, "_", measure), inherits = TRUE) %>% do.call(rbind, .)
  rownames(out_com) <- NULL
  out_com <- out_com[!duplicated(out_com), ]
  if(data_label == "tourism"){
    candidates <- c("OLS", "WLSs", "WLSv", "MinTs")
  } else if (data_label == "tourism"){
    candidates <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")
  }
  target <- c("Base", "BU", 
              sapply(candidates, 
                     function(len) c(len, paste0(len, "-", c("subset", "intuitive", "lasso")))) |> as.vector(), 
              "Elasso")
  out_com <- out_com[match(target, out_com$Method), ]
  colnames_out <- colnames(out_com)
  table_out <- sapply(1:NROW(out_com), function(i){
    r1 <- out_com[1, -1] |> as.numeric()
    rt <- out_com[i, -1] |> as.numeric()
    if (i==1){
      return(r1)
    } else{
      return(100*(rt - r1)/r1)
    }
  }) |> t()
  table_out <- data.frame(out_com[, 1], table_out)
  colnames(table_out) <- colnames_out
  return(list(table_out = table_out, levels = levels, candidates = candidates))
}

## Output table latex
latex_table <- function(out_all){
  out <- head(out_all$table_out, -1)
  levels <- out_all$levels
  header <- c("", rep(as.character(length(horizons)), length(levels)))
  names(header) <- c("", levels)
  candidates <- out_all$candidates
  
  # Red entries identify the best performing approaches
  # Bold entries identify methods that perform better than the corresponding benchmark method
  out[, -1] <- lapply(out[, -1], function(x) {
    x <- round(x, 1)
    origin_x <- x
    min_x <- min(origin_x)
    comp_x <- c(origin_x[1:2], rep(origin_x[3+4*(0:(length(candidates)-1))], each = 4))
    out_x <- ifelse(x == min(x), cell_spec(format(x, nsmall = 1), bold = TRUE, color = "red"), ifelse(x < comp_x, cell_spec(format(x, nsmall = 1), bold = TRUE), format(x, nsmall = 1)))
  })
  
  out |>
    kable(format = "latex",
          booktabs = TRUE,
          digits = 1,
          align = c("l", rep("r", length(horizons)*length(levels))),
          escape = FALSE,
          linesep = "") |>
    row_spec(2+4*(0:(length(candidates)-1)), hline_after = TRUE) |>
    # kable_paper("striped", full_width = F) |>
    kable_styling(latex_options = c("hold_position", "repeat_header", "scale_down")) |>
    row_spec(grepl("-", out$Method) |> which(), 
             background = "#e6e3e3") |>
    add_header_above(header, align = "c") |> print()
}

## Extract element from forecast results
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

## Calculate RMSE for each series
calc_rmse <- function(fc, test, h){
  err <- subset(test, select = -Index) - subset(fc, select = -Index)
  err <- cbind(err, Index = subset(test, select = Index))
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

## Calculate absolute errors for each series across different horizons
calc_abserror <- function(fc, test){
  abserr <- abs(subset(test, select = -Index) - subset(fc, select = -Index))
  abserr |> as.data.frame() |> data.matrix() |> as.vector()
}

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
# Table: Number of time series retained in the structure with different methods
#----------------------------------------------------------------------
tourism_subset_reconsf <- readRDS(file = "data_new/tourism_subset_reconsf.rds")
Top <- 1
State <- 2:8
Zone <- 9:35
Region <- 36:111
Total <- 1:111

tourism_subsetonly_reconsf <- tourism_subset_reconsf[[1]][grepl("_subset", names(tourism_subset_reconsf[[1]]))]
tourism_subset_z <- sapply(tourism_subsetonly_reconsf, function(lentry) lentry$z) |> t()
tourism_subset_lambda <- sapply(tourism_subsetonly_reconsf, function(lentry) lentry$lambda_opt) |> t()
tourism_subset_info <- data.frame(tourism_subset_z, tourism_subset_lambda)

None <- c(Top = length(Top), State = length(State), Zone = length(Zone), Region = length(Region),
          Total = length(Total), "$\\lambda_0$" = 0, "$\\lambda_2$" = 0)
tourism_subset_info <- apply(tourism_subset_info, 1, function(lentry){
  c(Top = sum(lentry[Top]), State = sum(lentry[State]), 
    Zone = sum(lentry[Zone]), Region = sum(lentry[Region]), 
    Total = sum(lentry[Total]), 
    "$\\lambda_0$" = round(tail(lentry, 2), 2)[1],
    "$\\lambda_2$" = round(tail(lentry, 2), 2)[2])
}) %>% t() %>% rbind(None, .)
rownames(tourism_subset_info) <- sub("_", "-", rownames(tourism_subset_info))

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
nemenyi(RMSE, conf.level = 0.95, plottype = "vmcb", title = "RMSE - Tourism data")

# Plot 2: Absolute error across all horizons (111*12 instances)
abserr_ts_horizon <- sapply(methods, function(lmethod){
  assign(lmethod, calc_abserror(fc = get(tolower(lmethod)), test = test))
}, USE.NAMES = TRUE)
nemenyi(abserr_ts_horizon, conf.level = 0.95, plottype = "vmcb", title = "Absolute error across all horizons - Tourism data")



