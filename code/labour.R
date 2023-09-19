library(tidyverse)
library(tsibble)
library(fabletools)
library(magrittr)
library(lubridate)
library(fable)
library(dplyr)
library(hts)

#----------------------------------------------------------------------
# ABS - Unemployed persons by Duration of job search, State and Territory
##
## 6291.0.55.001 - UM2 - Unemployed persons by Duration of job search, State and Territory, January 1991 onwards
## 
## Monthly series
## Duration of job search (Duration, 6) * State and territory (STT, 8): n = 63 series in total, nb = 48 series at the bottom level
##
## Training set:  2010Jan-2022Jul
## Test set:      2022Aug-2023Jul
#----------------------------------------------------------------------
# Import data
labour <- readr::read_csv("data/labour.csv") |>
  mutate(Month = yearmonth(Month)) |>
  filter(year(Month) >= 2010) |>
  dplyr::rename(Duration = `Duration of job search`, 
                STT = `State and territory (STT): ASGS (2011)`,
                Unemployed = `Unemployed total ('000)`) |>
  as_tsibble(key = c(Duration, STT),
             index = Month) |>
  relocate(Month) |>
  fill_gaps()

labour_fill <- labour |>
  # use a random walk to give linear interpolation between points
  model(naive = ARIMA(Unemployed ~ -1 + pdq(0,1,0) + PDQ(0,0,0))) |>
  interpolate(labour)

labour_data <- labour_fill |>
  mutate(Duration=replace(Duration, Duration=="Under 4 weeks (under 1 month)", "D1"),
         Duration=replace(Duration, Duration=="4 weeks and under 13 weeks (1-3 months)", "D2"),
         Duration=replace(Duration, Duration=="13 weeks and under 26 weeks (3-6 months)", "D3"),
         Duration=replace(Duration, Duration=="26 weeks and under 52 weeks (6-12 months)", "D4"),
         Duration=replace(Duration, Duration=="52 weeks and under 104 weeks (1-2 years)", "D5"),
         Duration=replace(Duration, Duration=="104 weeks and over (2 years and over)", "D6")) |>
  mutate(STT=replace(STT, STT=="Australian Capital Territory", "ACT"),
         STT=replace(STT, STT=="New South Wales", "NSW"),
         STT=replace(STT, STT=="Northern Territory", "NTT"),
         STT=replace(STT, STT=="Queensland", "QLD"),
         STT=replace(STT, STT=="South Australia", "SAS"),
         STT=replace(STT, STT=="Tasmania", "TAS"),
         STT=replace(STT, STT=="Victoria", "VIC"),
         STT=replace(STT, STT=="Western Australia", "WAS")) |>
  mutate(ID = paste0(Duration, STT)) |>
  arrange(ID, Month)

labour_bts <- pull(labour_data, Unemployed) |>
  matrix(nrow = length(unique(labour_data$Month)), ncol = 48, byrow = FALSE)
colnames(labour_bts) <- unique(labour_data$ID)

# Grouped time series
labour_gts <- gts(labour_bts, character=c(2,3),
                  gnames = c("Duration",
                             "STT"))
S <- smatrix(labour_gts)
labels <- do.call(c, labour_gts$labels) |> 
  as.character() %>%
  c("Total", ., colnames(labour_gts$bts))

labour_gts <- (labour_gts$bts %*% t(S)) %>% 
  data.frame(1, unique(labour_data$Month), .)
colnames(labour_gts) <- c("Index", "Time", labels)

saveRDS(labour_gts, file = "data/labour_data.rds")
saveRDS(S, file = "data/labour_S.rds")
