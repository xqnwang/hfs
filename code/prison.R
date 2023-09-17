library(tidyverse)
library(tsibble)
library(fabletools)
library(magrittr)
library(lubridate)
library(fable)
library(dplyr)
library(hts)

prison <- readr::read_csv("https://OTexts.com/fpp3/extrafiles/prison_population.csv") |>
  mutate(Quarter = yearquarter(Date)) |>
  select(-Date)  |>
  as_tsibble(key = c(Gender, Legal, State, Indigenous),
             index = Quarter) |>
  relocate(Quarter)

prison_bts <- prison |>
  aggregate_key(Gender * Legal * State, Count = sum(Count)) |>
  filter(!is_aggregated(Gender), !is_aggregated(Legal),
         !is_aggregated(State)) |>
  mutate(ID = paste0(substr(Gender, start = 1, stop = 1), 
                     substr(Legal, start = 1, stop = 1),
                     substr(State, start = 1, stop = 2))) |>
  arrange(ID, Quarter)
Quarter <- prison_bts$Quarter |> unique()

prison_data <- pull(prison_bts, Count) |>
  matrix(nrow = 48, ncol = 32, byrow = FALSE)
colnames(prison_data) <- unique(prison_bts$ID)

prison_gts <- gts(prison_data, character=c(1,1,2),
               gnames = c("Gender",
                          "Legal",
                          "State",
                          "Gender x Legal",
                          "Gender x State",
                          "Legal x State"))
S <- smatrix(prison_gts)
labels <- do.call(c, prison_gts$labels) |> 
  as.character() %>%
  c("Total", ., colnames(prison_gts$bts))

prison_gts <- (prison_gts$bts %*% t(S)) %>% 
  data.frame(1, Quarter, .)
colnames(prison_gts) <- c("Index", "Time", labels)

saveRDS(prison_gts, file = "data/prison_data.rds")
saveRDS(S, file = "data/prison_S.rds")

