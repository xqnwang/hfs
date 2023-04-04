library(fable)
library(dplyr)
library(lubridate)
library(purrr)
library(tsibble)
library(rlang)
library(ggplot2)

source('R/sourceDir.R')
sourceDir("R", recursive = TRUE)

#----------------------------------------------------------------------
# Australian domestic tourism (only considering hierarchical structure)
#
## Quarterly series from 1998Q1-2017Q4: 80 observations for each series
## Total/State/Region: 2 levels, n = 85
##    Training set:  1998Q1-2015Q4
##    Test set:      2016Q1-2017Q4
#----------------------------------------------------------------------
tourism <- tsibble::tourism |>
  mutate(State = recode(State,
                        `New South Wales` = "NSW",
                        `Northern Territory` = "NT",
                        `Queensland` = "QLD",
                        `South Australia` = "SA",
                        `Tasmania` = "TAS",
                        `Victoria` = "VIC",
                        `Western Australia` = "WA"
  ))

tourism_hts <- tourism |>
  aggregate_key(State, Trips = sum(Trips))
tourism_hts
# tourism_hts |> distinct(State)

tourism_hts |>
  autoplot(Trips) +
  labs(y = "Trips ('000)",
       title = "Australian tourism: national and states") +
  facet_wrap(vars(State), scales = "free_y", ncol = 3) +
  theme(legend.position = "none")

lambda <- c(0.5, 0, 0)
subset <- TRUE
lasso <- FALSE
ridge <- FALSE

## cross-validation
# tourism_hts_tr <- tourism_hts |>
#   stretch_tsibble(.init = 72, .step = 1) |>
#   relocate(Quarter, State, .id)

# base forecast
fit_base_VIC <- tourism_hts |>
  filter(year(Quarter) <= 2015) |>
  filter(State == "VIC") |>
  model(base = RW(Trips ~ drift()))
fit_base_others <- tourism_hts |>
  filter(year(Quarter) <= 2015) |>
  filter(State != "VIC") |>
  model(base = ETS(Trips))
fit_base <- dplyr::bind_rows(fit_base_VIC, fit_base_others)

# reconciliation
fit_rec <- fit_base |>
  reconcile(
    BU = bottom_up(base),
    OLS = min_trace(base, method = "ols"),
    OLS_s = min_trace_subset(base, method = "ols", subset, lasso, ridge, lambda = lambda, sparse = FALSE),
    WLSv = min_trace(base, method = "wls_var"),
    WLSv_s = min_trace_subset(base, method = "wls_var", subset, lasso, ridge, lambda = lambda, sparse = FALSE),
    WLSs = min_trace(base, method = "wls_struct"),
    WLSs_s = min_trace_subset(base, method = "wls_struct", subset, lasso, ridge, lambda = lambda, sparse = FALSE),
    MinT = min_trace(base, method = "mint_shrink"),
    MinT_s = min_trace_subset(base, method = "mint_shrink", subset, lasso, ridge, lambda = lambda, sparse = FALSE)
  )

fc <- fit_rec |>
  forecast(h = "2 years")

fc |>
  accuracy(
    data = tourism_hts,
    measures = list(rmse = RMSE, mase = MASE)
  ) |>
  group_by(.model) |>
  summarise(rmse = mean(rmse), mase = mean(mase))

model_name <- fc |> distinct(.model) |> pull(.model)
tab <- matrix(NA, ncol = length(model_name)*2, nrow = 3)
rownames(tab) <- c("Total", "State", "All series")
filter_tab <- matrix(NA, ncol = 1, nrow = 3)
filter_tab[1] <- "fc %>% filter(is_aggregated(State))"
filter_tab[2] <- "fc %>% filter(!is_aggregated(State))"
filter_tab[3] <- "fc"
for (i in 1:3) {
  eval(parse(text = filter_tab[i])) %>%
    accuracy(data = tourism_hts, measures = list(rmse = RMSE, mase = MASE)) %>%
    group_by(.model) %>%
    summarise(rmse = mean(rmse), mase = mean(mase)) -> err
  err <- err |> arrange(factor(.model, levels = model_name))
  tab[i, ] <- cbind(t(err[, 2]), t(err[, 3]))
}
colnames(tab) <- err |> pull(.model) |> rep(2)
out <- list(RMSE = tab[, 1:length(model_name)], MASE = tab[, (length(model_name)+1):(length(model_name)*2)])
out$RMSE |> round(2)
out$MASE |> round(2)




