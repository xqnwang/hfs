library(fable)
library(dplyr)
library(tidyr)
library(tsibble)
library(tsibbledata)
library(lubridate)
library(purrr)
library(ggplot2)

source('R/sourceDir.R')
sourceDir("R", recursive = TRUE)

get_structure <- function(model){
  key_data <- key_data(model)
  agg_data <- fabletools:::build_key_data_smat(key_data)
  S <- matrix(0L, 
              nrow = length(agg_data$agg), 
              ncol = max(vctrs::vec_c(!!!agg_data$agg)))
  S[length(agg_data$agg)*(vctrs::vec_c(!!!agg_data$agg)-1) + rep(seq_along(agg_data$agg), lengths(agg_data$agg))] <- 1L
  return(list(key_data = key_data, S = S))
}

#----------------------------------------------------------------------
# Australian domestic tourism (only considering hierarchical structure)
## Quarterly series from 1998Q1-2017Q4: T=80
## Total/State/Region: 2 levels, n=85
### training set: 1998Q1-2007Q4 (10Y*4)
### validation set: 2008Q1-2012Q4 (5Y*4)
### test set: 2013Q1-2017Q4 (5Y*4)
### one-step-ahead rolling-origin base forecasts
#----------------------------------------------------------------------

# tourism data
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
#tourism_hts <- tourism |>
#  aggregate_key(State / Region, Trips = sum(Trips))
tourism_hts |> distinct(State)


# setup
Tw <- 40 # length of window
T_validation <- T_test <- 20 # length of validation and test period
h <- 1 # forecast horizon


# make data windows
index <- tourism_hts |> 
  distinct(Quarter) |>
  pull(Quarter)
data_windows <- purrr::map(1:(length(index) - Tw + 1 - h),
                           function(i){
                             index_w <- index[i:(i + Tw - 1)]
                             return(tourism_hts |> 
                                      filter(Quarter %in% index_w))
                           })
data_windows[[1]] |> head(3)
data_windows[[20]] |> tail(3)
data_windows[[40]] |> tail(3)


# benchmark S and key_data
st_ben <- data_windows[[1]] |>
  model(base = ETS(Trips)) |>
  get_structure()
S_ben <- st_ben$S
key_ben <- st_ben$key_data


# base modeling and forecasting
baseforec_window <- function(data_w, key_ben){
  md <- data_w |>
    model(base = ETS(Trips))
  fc <- md |>
    forecast(h = 1)
  ky <- key_data(md)
  if(all.equal(ky, key_ben)){
    mean <- fc$.mean
  }else{
    ky |> 
      mutate(Row = seq.int(NROW(ky))) |>
      select(-.rows) -> ky
    left_join(key_ben, ky, by = join_by(State, Region)) |> 
      pull(Row) -> ky_index
    mean <- fc$.mean[ky_index]
  }
  return(mean)
}


# base forecasts on validation set
data_windows[1:T_validation] |> 
  purrr::map(\(x) baseforec_window(x, key_ben = key_ben)) -> basefc_validation
#> saveRDS(basefc_validation, file = "data/tourism_basefc_validation.rds")
#> basefc_validation <- readRDS("data/tourism_basefc_validation.rds")


# extract observations on validation set
tourism_validation <- tourism_hts |> 
  filter(year(Quarter) %in% 2008:2012) # 2008Q1-2012Q4
index_validation <- tourism_validation |> 
  distinct(Quarter) |>
  pull(Quarter)
obs_validation <- purrr::map(1:length(index_validation), function(i){
  tv_validation <- tourism_validation |> 
    filter(Quarter == index_validation[i])
  tv_validation <- left_join(key_ben, tv_validation, 
                             by = join_by(State)) |> 
    pull(Trips)
  #tv_validation <- left_join(key_ben, tv_validation, 
  #                           by = join_by(State, Region)) |> 
  #  pull(Trips)
  return(tv_validation)
  })


# forecast reconciliation approach (calculate matrix G)
obs_valid <- matrix(unlist(obs_validation), 
                    nrow = length(obs_validation), 
                    byrow = TRUE)
basefc_valid <- matrix(unlist(basefc_validation), 
                       nrow = length(basefc_validation), 
                       byrow = TRUE)
duplicated_index <- which(duplicated(obs_valid, MARGIN = 2))
if(!identical(duplicated_index, integer(0))){
  obs_valid <- obs_valid[, -duplicated_index]
  basefc_valid <- basefc_valid[, -duplicated_index]
  S <- S_ben[-duplicated_index, ]
}else{
  S <- S_ben
}

y <- as.vector(obs_valid)
X <- kronecker(S, basefc_valid)
group_indices <- lapply(seq.int(NROW(S)), function(j) {
  j + (0:(NCOL(S)-1)) * NROW(S)
  })

model <- mip_l0l2(x = X, y = y, group_indices = group_indices, 
                  lambda_0 = 500, lambda_2 = 200, M = NCOL(S)*10, 
                  solver = "gurobi")

b <- model$coefficients
z <- model$binary
z
G <- matrix(b, nrow = NROW(S), ncol = NCOL(S), byrow = FALSE) |>
  t()
G

# testing
## observation data on test set
tourism_test <- tourism_hts |> 
  filter(year(Quarter) %in% 2013:2017) # 2013Q1-2017Q4
index_test <- tourism_test |> 
  distinct(Quarter) |>
  pull(Quarter)
obs_test <- purrr::map(1:length(index_test), function(i){
  tv_test <- tourism_test |> 
    filter(Quarter == index_test[i])
  tv_test <- left_join(key_ben, tv_test, 
                             by = join_by(State)) |> 
    pull(Trips)
  # tv_test <- left_join(key_ben, tv_test, 
  #                      by = join_by(State, Region)) |> 
  #   pull(Trips)
  return(tv_test)
})
obs_test <- matrix(unlist(obs_test), 
                    nrow = length(obs_test), 
                    byrow = TRUE)
if(!identical(duplicated_index, integer(0))){
  obs_test <- obs_test[, -duplicated_index]
}

## base forecasts and other reconciled forecasts (bu, td, ols, mint) on test set
reconcile_window <- function(data_w, key_ben, S, G){
  md <- data_w |>
    model(base = ETS(Trips)) |>
    reconcile(
      td = top_down(base, method = "forecast_proportions"),
      bu = bottom_up(base),
      ols = min_trace(base, method = "ols"),
      mint = min_trace(base, method = "mint_shrink")
    )
  fc <- md |>
    forecast(h = 1)
  ky <- key_data(md)
  if(all.equal(ky, key_ben)){
    out <- fc |>
      pivot_wider(id_cols = c(State, Quarter), names_from = .model, values_from = .mean) |>
      mutate(l0l2 = as.numeric(S %*% G %*% base))
  }else{
    ky |> 
      mutate(Row = seq.int(NROW(ky))) |>
      select(-.rows) -> ky
    left_join(key_ben, ky, by = join_by(State)) |> 
      pull(Row) -> ky_index
    # left_join(key_ben, ky, by = join_by(State, Region)) |> 
    #   pull(Row) -> ky_index
    out <- fc[ky_index, ] |>
      pivot_wider(id_cols = c(State, Quarter), names_from = .model, values_from = .mean) |>
      mutate(l0l2 = as.numeric(S %*% G %*% base))
  }
  return(out)
}
reconcile_test <- data_windows[(1+T_validation):(T_validation + T_test)] |> 
  purrr::map(\(x) reconcile_window(x, key_ben = key_ben, S, G))
out_test <- reduce(reconcile_test, bind_rows) # Reduce a list of tsibble objects to a tsibble
out_test <- left_join(out_test, tourism_hts, by = join_by(State, Quarter))
out_test |> 
  pivot_longer(cols = 3:9, names_to = "Method", values_to = "Trip") |>
  ggplot(aes(x = Quarter, y = Trip, colour = Method, group = Method)) +
  stat_summary(fun = sum, geom = "line") + 
  facet_wrap(~ as.character(State),
             nrow = 1, scales = "free_y")


# accuracy
tab <- vector(mode = "list", length = 2)
tab[[1]] <- tab[[2]] <- matrix(NA, ncol = 6, nrow = 3)
rownames(tab[[1]]) <- rownames(tab[[2]]) <- c("Total", "State", "All series")
colnames(tab[[1]]) <- colnames(tab[[2]]) <- c("Base", "Top-down", "Bottom-up", "OLS", "MinT", "L0L2")
names(tab) <- c("RMSE", "MASE")

filter_tab <- matrix(NA, ncol = 1, nrow = 3)
filter_tab[1] <- "out_test %>% filter(is_aggregated(State))"
filter_tab[2] <- "out_test %>% filter(!is_aggregated(State))"
filter_tab[3] <- "out_test"

train_tab <- vector(mode = "list", length = 3)
train_tab[[1]] <- data_windows[(1+T_validation):(T_validation + T_test)] |>
  purrr::map(\(x) x |> filter(is_aggregated(State)) |> as_tibble()) |>
  reduce(bind_rows) |> unique()
train_tab[[2]] <- data_windows[(1+T_validation):(T_validation + T_test)] |>
  purrr::map(\(x) x |> filter(!is_aggregated(State)) |> as_tibble()) |>
  reduce(bind_rows) |> unique()
train_tab[[3]] <- data_windows[(1+T_validation):(T_validation + T_test)] |>
  purrr::map(\(x) x |> as_tibble()) |>
  reduce(bind_rows) |> unique()

for (i in 1:3) {
  filter_tab_i <- eval(parse(text = filter_tab[i]))
  train_tab_i <- lapply(seq.int(NROW(filter_tab_i)), function(r){
    State_r <- filter_tab_i[r, ] |> pull(State)
    Quarter_r <- filter_tab_i[r,] |> pull(Quarter)
    train_tab_i_r <- train_tab[[i]] |> 
      filter((State == State_r) & 
              (Quarter >= (Quarter_r - 40)) & 
              (Quarter <= (Quarter_r-1)))
  })
  err_tab_i <- filter_tab_i[, c("base", "td", "bu", "ols", "mint", "l0l2")] |>
    purrr::map(\(x) x- filter_tab_i[, "Trips"]) |>
    reduce(bind_cols) |>
    as.matrix()
  RMSE_tab_i <- apply(err_tab_i, c(1,2), RMSE) %>% colMeans()
  MASE_tab_i <- sapply(seq.int(NROW(err_tab_i)), function(r){
    err_r <- err_tab_i[r, ] |> as.matrix()
    train_r <- train_tab_i[[r]] |> pull(Trips)
    mase_r <- apply(err_r, 1, 
                    function(e) MASE(.resid = e, .train = train_r, .period = 4))
  }) |> t() |> colMeans()
  tab[[1]][i, ] <- cbind(RMSE_tab_i)
  tab[[2]][i, ] <- cbind(MASE_tab_i)
}
tab
