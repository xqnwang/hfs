# Load packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(knitr)
library(kableExtra)
library(latex2exp)

# Functions
combine_table <- function(data_label, methods, scenario = NULL, horizons){
  for(method_label in methods){
    data_method <- paste0(data_label, "_", method_label)
    for(h in horizons){
      assign(paste0(data_method, "_rmse_h", h), 
             readRDS(file = paste0("data_new/", data_method, "_reconsf", ifelse(is.null(scenario), "", paste0("_", scenario)), "_rmse_", h, ".rds"))
      )
    }
    example <- get(paste0(data_method, "_rmse_h", horizons[1]))
    levels <- colnames(example)[-1]
    method <- sub("_", "-", example$Method)
    rmse <- lapply(levels, function(level){
      mget(paste0(data_method, "_rmse_h", horizons), inherits = TRUE) |> 
        sapply(function(lentry) as.numeric(lentry[[level]]))
    })
    rmse <- data.frame(method, do.call(cbind, rmse))
    colnames(rmse) <- c("Method", c(ifelse(horizons == 1, paste0("h=", 1), paste0("1-", horizons))) |> rep(length(levels)))
    assign(paste0(data_method, "_rmse"), rmse)
  }
  
  rmse_out <- mget(paste0(data_label, "_", methods, "_rmse"), inherits = TRUE) %>% do.call(rbind, .)
  rownames(rmse_out) <- NULL
  rmse_out <- rmse_out[!duplicated(rmse_out), ]
  if(data_label == "tourism"){
    candidates <- c("OLS", "WLSs", "WLSv", "MinTs")
  } else if (data_label == "tourism"){
    candidates <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")
  }
  target <- c("Base", "BU", 
              sapply(candidates, 
                     function(len) c(len, paste0(len, "-", c("subset", "intuitive", "lasso")))) |> as.vector(), 
              "Elasso")
  rmse_out <- rmse_out[match(target, rmse_out$Method), ]
  colnames_out <- colnames(rmse_out)
  table_out <- sapply(1:NROW(rmse_out), function(i){
    r1 <- rmse_out[1, -1] |> as.numeric()
    rt <- rmse_out[i, -1] |> as.numeric()
    if (i==1){
      return(r1)
    } else{
      return(100*(rt - r1)/r1)
    }
  }) |> t()
  table_out <- data.frame(rmse_out[, 1], table_out)
  colnames(table_out) <- colnames_out
  return(list(table_out = table_out, levels = levels, candidates = candidates))
}

# Target results to be processed
data_label <- "tourism"
scenario <- NULL
horizons <- c(1, 4, 8, 12)
methods <- c("subset", "intuitive", "lasso")

# Combine results of different methods
out_all <- combine_table(data_label, methods, scenario, horizons)
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

# Latex table output
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
  add_header_above(header, align = "c")
