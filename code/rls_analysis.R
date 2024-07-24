# load libraries -------
library(tidyverse)

# set working directory --------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in data
rls_raw <- read_csv("../data/rls_raw.csv")

rls_clean <- rls_raw %>%
  group_by(survey_date, latitude, longitude, ecoregion, species_name) %>%
  summarise(total_biomass = sum(biomass))

# bring in sea surface temperature data