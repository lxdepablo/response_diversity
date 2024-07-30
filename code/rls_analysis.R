# load libraries -------
library(tidyverse)

# set working directory --------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in data
rls_raw <- read_csv("../data/rls_raw.csv")

rls_clean <- rls_raw %>%
  group_by(survey_date, latitude, longitude, ecoregion, species_name) %>%
  summarise(total_biomass = sum(biomass))

# export coords for each site
coords <- rls_clean %>%
  ungroup() %>%
  select(c(survey_date, latitude, longitude)) %>%
  unique()

write_csv(coords, "../data/rls_coords.csv")

# bring in sea surface temperature data
sst_raw <- read_csv("../data/sst.csv")

# extract coordinates with regex
extract_values <- str_extract_all(sst_raw$.geo, "\\[(.*?)\\]")
extracted_values <- gsub("\\[|\\]", "", unlist(extract_values)) # Remove the brackets
# separate out latitude and longitude
separated_coords <- str_split(extracted_values, ",")

coords_df <- bind_rows(lapply(separated_coords, function(x){
  data.frame(longitude = x[[1]], latitude = x[[2]])
}))

# clean up SST data
sst_clean <- sst_raw %>%
  cbind(coords_df) %>%
  filter(sst != -9999) %>%
  select(c(survey_date, latitude, longitude, sst)) %>%
  mutate(latitude = trunc(as.numeric(latitude)*100)/100,
         longitude = trunc(as.numeric(longitude)*100)/100) %>%
  unique()

# join to fish data
rls_temp <- rls_clean %>%
  left_join(sst_clean, by = c("survey_date", "latitude", "longitude")) %>%
  drop_na()

# plot abundance for a couple fish
# get all fish
all_fish_sp <- unique(rls_temp$species_name)
# pick some random fish to plot
rand_fish <- sample(all_fish_sp, 5)

rls_temp %>%
  filter(species_name %in% rand_fish) %>%
  ggplot(aes(x = sst, y = total_biomass, col = species_name)) +
  geom_point()

# fit abundance functions for every species of fish




