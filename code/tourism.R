library(tidyverse)
library(fabletools)
library(magrittr)
library(lubridate)
library(fable)

#----------------------------------------------------------------------
# Australian domestic tourism (only considering hierarchical structure)
##
## Quarterly series from 1998Q1-2017Q4: 80 observations for each series
##
## Total/State: 2 levels, n = 9
## Total/State/Region: 3 levels, n = 85
##
## Training set:  1998Q1-2015Q4
## Test set:      2016Q1-2017Q4
#----------------------------------------------------------------------
# target data
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

# time series plot
tourism_hts |>
  autoplot(Trips) +
  labs(y = "Trips ('000)",
       title = "Australian tourism: national and states") +
  facet_wrap(vars(State), scales = "free_y", ncol = 3) +
  theme(legend.position = "none")

# data construction
tourism_hts <- tourism_hts |>
  as_tibble() |>
  pivot_wider(names_from = State,
              values_from = Trips
              ) |>
  mutate(Index = 1,
         Time = Quarter,
         Total = `<aggregated>`,
         .before = `<aggregated>`) |>
  select(!c(Quarter, `<aggregated>`))
saveRDS(tourism_hts, file = "data/tourism_data.rds")
