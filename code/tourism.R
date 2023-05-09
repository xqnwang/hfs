library(tidyverse)
library(tsibble)
library(fabletools)
library(magrittr)
library(lubridate)
library(fable)
library(hts)

#----------------------------------------------------------------------
# Australian domestic tourism (only considering hierarchical structure)
##
## Spiliotis et al. (2021) on ASC
##
## Monthly series from 1998Jan-2017Dec: 240 months (20 years) for each series
##
## Total/State/Zone/Region: 4 levels, n = 111 series in total
##
## Training set:  1998Jan-2016Dec
## Test set:      2017Jan-2017Dec
#----------------------------------------------------------------------
# Import data, only include bottom-level series
tourism <- readr::read_csv("data/TourismData-asc2021.csv", skip = 3) |>
  select(-1) |>
  slice(-(1:2)) |>
  rename(Year = ...2, Month = ...3) |>
  fill(Year, .direction = "down") |>
  mutate(Month = str_sub(Month, 1, 3)) |>
  mutate(across(AAA:GBD, as.numeric)) |>
  unite("Time", Year:Month, sep = " ", remove = TRUE)


# Hierarchical time series
tourism_hts <- hts(tourism[-1] |> ts(start = c(1998, 1), end = c(2017, 12), frequency = 12), 
                   characters = c(1, 1, 1))
labels <- do.call(c, tourism_hts$labels) |> 
  as.character()
S <- smatrix(tourism_hts)
tourism_hts <- (tourism_hts$bts %*% t(S)) %>% 
  data.frame(1, tourism |> pull(Time), .)
colnames(tourism_hts) <- c("Index", "Time", labels)
saveRDS(tourism_hts, file = "data/tourism_data.rds")
saveRDS(S, file = "data/tourism_data_S.rds")

# Time series plot
tourism_ts <- tourism_hts |>
  select(-1) |>
  as_tibble() |>
  pivot_longer(cols = Total:GBD, names_to = "Series", values_to = "Value") |>
  mutate(Time = yearmonth(Time)) |>
  as_tsibble(index = Time, key = Series)
tourism_ts |>
  #filter(Series %in% c("Total", LETTERS[1:7])) |>
  autoplot(Value) +
  facet_wrap(vars(Series), scales = "free_y") +
  theme(legend.position = "none")
