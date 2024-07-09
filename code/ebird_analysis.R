# load libraries
library(tidyverse)
library(auk)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# read in sample ebird data
auk_set_ebd_path("data/ebird_sample")

# read in ebird data
# ebird_raw <- auk_ebd("ebd_US-AL-101_202204_202204_relApr-2022.txt") %>%
#   # filter for USA
#   auk_country("United States") %>%
#   auk_filter(file = "ebird_filtered.txt") %>%
#   # read in data from text file
#   read_ebd()

# built in sample data
ebird_raw <- system.file("extdata/ebd-sample.txt", package = "auk") %>%
  auk_ebd() %>%
  # filter for USA
  auk_country("United States") %>%
  auk_filter(file = "test_filtered.txt", overwrite = TRUE) %>%
  # read in data from text file
  read_ebd()

# drop unnecessary columns
ebird_clean <- ebird_raw %>%
  select(c(state, county_code, observation_date, scientific_name, common_name)) %>%
  # format county codes
  mutate(county_code = substr(county_code, 4, length(county_code)))

# read in NOAA temperature data
temps_raw <- read_csv("data/ebird_sample/april_temps.csv")

# join temperature data to bird observations
ebird_environment <- ebird_clean %>%
  left_join(temps_raw, by = c("county_code" = "ID")) %>%
  select(c(state, Name, county_code, observation_date, scientific_name, common_name, Value)) %>%
  rename(county = Name,
         avg_monthly_temp = Value) %>%
  mutate(avg_monthly_temp = as.numeric(avg_monthly_temp)) %>%
  drop_na()

# plot response functions for each species
ggplot(data = ebird_environment, aes(x = avg_monthly_temp)) +
  geom_histogram() +
  facet_wrap(~scientific_name) +
  labs(x = "Average Temperature (F)", y = "Number of Observations") +
  theme_bw() +
  theme(panel.grid = element_blank())

