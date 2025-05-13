# load libraries
library(tidyverse)
library(janitor)
library(lme4)
library(fixest)
library(BBmisc)
library(sjPlot)
library(cowplot)
library(parallel)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# get some helper functions
source("cc_utils.R")

# read in data
climate_raw <- read_csv("../data/cedar_creek/climate.csv") %>%
  clean_names()
biomass_raw <- read_csv("../data/cedar_creek/main_biomass.csv") %>%
  clean_names()
# soil trait data
phosphorous_raw <- read_tsv("../data/cedar_creek/soil_phosphorous.txt") %>%
  clean_names()
nitrate_ammonium_raw <- read_tsv("../data/cedar_creek/soil_nitrate_ammonium.txt") %>%
  clean_names()
nitrogen_raw <- read_tsv("../data/cedar_creek/soil_nitrogen.txt") %>%
  clean_names

# wrangle soil trait data
phosphorous_clean <- phosphorous_raw %>%
  select(c(year, plot, depth, phosphorus)) %>%
  filter(depth == "0.20cm") %>%
  select(-depth) %>%
  group_by(plot) %>%
  summarize(phosphorus = mean(phosphorus, na.rm=TRUE))
nitrate_ammonium_clean <- nitrate_ammonium_raw %>%
  select(c(year, plot, depth, no2no3, nh4)) %>%
  filter(depth == "0-20") %>%
  select(-depth) %>%
  group_by(plot) %>%
  summarize(no2no3 = mean(no2no3, na.rm=TRUE),
            nh4 = mean(nh4, na.rm=TRUE))
nitrogen_clean <- nitrogen_raw %>%
  select(c(year, plot, depth, p_total_ns)) %>%
  filter(depth == "0-20 cm") %>%
  select(-depth) %>%
  group_by(plot) %>%
  summarize(nitrogen = mean(p_total_ns, na.rm=TRUE))

# bring together all 4 soil traits
all_soil_traits <- phosphorous_clean %>%
  left_join(nitrate_ammonium_clean, by = "plot") %>%
  left_join(nitrogen_clean, by = "plot")

# normalize soil traits to 0-1 scale
soil_traits_norm <- all_soil_traits %>%
  mutate(phosphorus = normalize(phosphorus, method = "range"),
         no2no3 = normalize(no2no3, method = "range"),
         nh4 = normalize(nh4, method = "range"),
         nitrogen = normalize(nitrogen, method = "range"))
  

# calculate start of growing season each year?
#  first consecutive 6 days with temp > 5 C
#  according to https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/eap.2641
climate_clean <- climate_raw %>%
  mutate(date = mdy(date),
         day = day(date),
         month = month(date),
         year = year(date)) %>%
  drop_na() %>%
  arrange(date)



# function to find start and end of growing season
growing_season_start <- bind_rows(lapply(unique(climate_clean$year), function(y){
  # get data for current year
  curr_year <- filter(climate_clean, year == y)
  
  consec_days <- 0
  end_date <- data.frame(date = NA, year = NA)
  
  # iterate over data and check temperatures
  for(i in 1:nrow(curr_year)){
    #print(curr_year[i, ]$max_temp_deg_c)
    if(curr_year[i, ]$max_temp_deg_c > 5){
      consec_days <- consec_days + 1
    }
    # if enough consecutive days over temp threshold, return growing season start date
    if(consec_days == 6){
      end_date <- select(curr_year[i, ], c(date, year))
      break
    }
  }
  end_date
})) %>% 
  drop_na() %>%
  rename(growing_season_start = date)

ggplot(data = climate_clean, aes(x = max_temp_deg_c)) +
  geom_histogram()

# drop unneeded columns from biomass data
biomass_clean <- biomass_raw %>%
  select(c(year, month, date, plot, strip, species, num_sp, sp_num, biomass_g_m2)) %>%
  mutate(planted = ifelse(species %in% c("Miscellaneous sp.",
                          "Miscellaneous grasses",
                          "Unsorted Biomass",
                          "Unsorted biomass",
                          "Miscellaneous litter",
                          "Miscellaneous herbs",
                          "Miscellaneous sedges",
                          "Mosses & lichens",
                          "Grasses",
                          "Mosses",
                          "miscellaneous seedhead",
                          "Miscellaneous Forb",
                          "Miscellaneous forb",
                          "unknown grass",
                          "unknown forb",
                          "Green matter (alive)",
                          "Miscellaneous grass",
                          "32 Species Weeds",
                          "Real Weeds",
                          "Green matter",
                          "Weeds"),
                0, 1)) %>%
  filter(biomass_g_m2 > 0,
         !is.na(date))

# find end of growing season
growing_season_by_plot <- biomass_clean %>%
  pivot_wider(names_from = species, values_from = biomass_g_m2) %>%
  mutate(date = mdy(date),
         growing_season_end = date-14) %>%
  select(c(plot, year, strip, growing_season_end)) %>%
  # drop 2001 since it was sampled in both june and july
  filter(year != 2001) %>%
  left_join(growing_season_start, by = "year")

# get mean temp and precip for every plot/year
# parallelized because this is slow
# create a cluster with 8 workers
numCores <- 8
cl <- makeCluster(numCores, type = "PSOCK")

# export objects to the cluster
clusterExport(
  cl,
  varlist = c("growing_season_by_plot", "climate_clean"),
  envir = environment()  # ensures the cluster can see these objects
)

# load required packages on each worker
clusterEvalQ(cl, {
  library(dplyr)
  library(lubridate)
})

# get environmental data for each growing season/plot
par_results <- parLapply(cl, 1:nrow(growing_season_by_plot), function(i){
  curr_sample <- growing_season_by_plot[i, ]

  start_date <- curr_sample$growing_season_start
  end_date   <- curr_sample$growing_season_end

  start_ind <- which(climate_clean$date == start_date)
  end_ind   <- which(climate_clean$date == end_date)

  # get climate data within date range
  gs_clim <- climate_clean[start_ind:end_ind, ]

  # summarize climate data for growing season
  gs_clim_sum <- gs_clim %>%
    summarize(
      mean_max_temp = mean(max_temp_deg_c, na.rm=TRUE),
      mean_min_temp = mean(min_temp_deg_c, na.rm=TRUE),
      mean_precip   = mean(precip_mm, na.rm=TRUE)
    ) %>%
    mutate(
      date = curr_sample$date,
      plot = curr_sample$plot,
      year = curr_sample$year
    )

  # Return the data frame for each iteration
  gs_clim_sum
})

# combine results into a single data frame
climate_growing_season <- bind_rows(par_results)
# stop the cluster
stopCluster(cl)


# # non parallel version
# # get environmental data for each growing season/plot
# env_results <- lapply(1:nrow(growing_season_by_plot), function(i){
#   curr_sample <- growing_season_by_plot[i, ]
# 
#   start_date <- curr_sample$growing_season_start
#   end_date   <- curr_sample$growing_season_end
# 
#   start_ind <- which(climate_clean$date == start_date)
#   end_ind   <- which(climate_clean$date == end_date)
#   print(paste0("start: ", start_ind, ", end: ", end_ind))
# 
#   # get climate data within date range
#   gs_clim <- climate_clean[start_ind:end_ind, ]
# 
#   # summarize climate data for growing season
#   gs_clim_sum <- gs_clim %>%
#     summarize(
#       mean_max_temp = mean(max_temp_deg_c, na.rm=TRUE),
#       mean_min_temp = mean(min_temp_deg_c, na.rm=TRUE),
#       mean_precip   = mean(precip_mm, na.rm=TRUE)
#     ) %>%
#     mutate(
#       plot = curr_sample$plot,
#       year = curr_sample$year
#     )
# 
#   # Return the data frame for each iteration
#   gs_clim_sum
# })
# 
# # combine results into a single data frame
# climate_growing_season <- bind_rows(env_results)

climate_growing_season_clean <- climate_growing_season %>%
  group_by(plot, year) %>%
  summarize(mean_max_temp = mean(mean_max_temp),
            mean_min_temp = mean(mean_min_temp),
            mean_precip = mean(mean_precip))

# visualize growing season climate data
ggplot(data = climate_growing_season_clean, aes(x = mean_max_temp)) +
  geom_histogram()

# join climate data to biomass
biomass_clim <- biomass_clean %>%
  left_join(climate_growing_season_clean, by = join_by("year", "plot")) %>%
  filter(year != 2001)


# plot climate time series
ggplot(data = filter(climate_clean, date %in% seq(mdy("01-01-1990"), mdy("01-01-2025"), by="day")), aes(x = date, y = max_temp_deg_c)) +
  geom_line() +
  geom_point(data = biomass_clim, aes(x = mdy(date), y = mean_max_temp, col = plot), col = "red") +
  scale_color_viridis_c() +
  theme_classic()


# calculate temp and precip optimum for each species
optima <- biomass_clim %>%
  filter(num_sp == 1,
         planted == 1) %>%
  group_by(species) %>%
  summarize(max_temp_optimum = weighted.mean(mean_max_temp, biomass_g_m2, na.rm = T),
            precip_optimum = weighted.mean(mean_precip, biomass_g_m2, na.rm = T)) %>%
  drop_na()

# prep data for fitting response curves
response_curve_df <- biomass_clim %>%
  # use only monoculture data
  filter(num_sp == 1,
         # only fit curves for planted species
         planted == 1) %>%
  # add column for day of year
  mutate(day_of_year = yday(mdy(date)))

# fit response curves for each species
temperature_fits <- bind_rows(lapply(unique(response_curve_df$species),
                                     fit_response_curve,
                                     df = response_curve_df))

temperature_fits <- response_curve_df %>%
  # make one tibble per species
  nest(data = -species) %>% 
  mutate(curve = map2(data, species, fit_response_curve_doy, center_doy=TRUE)) %>%
  filter(!map_lgl(curve, is.null)) %>%
  mutate(curve = map(curve, ~ select(.x, -species))) %>%
  unnest(curve)

# visualize response curves
temperature_fits %>%
  ggplot(aes(mean_max_temp, predicted, colour = species)) +
  geom_line() +
  theme_bw()
# get all species
all_sp <- unique(temperature_fits$species)
# pick some random species to plot
rand_sp <- sample(all_sp, 5)

curr_sp <- biomass_clim %>%
  filter(species %in% rand_sp)

curr_fits <- temperature_fits %>%
  filter(species %in% rand_sp)

# plot random subset of species
ggplot(data = curr_sp, aes(x = mean_max_temp, y = biomass_g_m2)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_vline(data = curr_fits, aes(xintercept = b), col = "blue", linetype = "dashed") +
  geom_hline(data = curr_fits, aes(yintercept = biomass_g_m2), col = "blue", linetype = "dashed") +
  geom_line(data = curr_fits, aes(x = mean_max_temp, y = predicted), col = "red") +
  labs(x = "Temperature", y = "Biomass (g/m2)") +
  facet_grid(~species) +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())

# check how many pancake-shaped curves there are
diagnostics <- temperature_fits %>% 
  distinct(species, a, b, c) %>% 
  mutate(span      = map_dbl(species, ~ diff(range(response_curve_df$mean_max_temp[response_curve_df$species == .x]))),
         flat_flag = c > 0.5 * span)
table(diagnostics$flat_flag)

# visualize optima weighted means
ggplot(data = optima, aes(x = max_temp_optimum)) +
  geom_histogram()
ggplot(data = optima, aes(x = precip_optimum)) +
  geom_histogram()

# summarize fits
fits_summary <- temperature_fits %>%
  group_by(species) %>%
  summarize(max_temp_optimum = mean(b))

# visualize optima from curves next to weighted means
ggplot() +
  geom_histogram(data = fits_summary, aes(x = max_temp_optimum), fill = "red", alpha = 0.5) +
  geom_histogram(data = optima, aes(x = max_temp_optimum), fill = "blue", alpha=0.5)

# calculate functional stability across environments
stability <- biomass_clim %>%
  # exclude monocultures
  filter(num_sp > 1) %>%
  group_by(year, month, plot, species) %>%
  # get mean function for each species for each plot
  summarize(mean_biomass_g_m2 = mean(biomass_g_m2),
            # use realized richness
            richness = mean(sp_num)) %>%
  # sum across species
  group_by(year, month, plot) %>%
  summarize(total_function = sum(mean_biomass_g_m2),
            richness = mean(richness)) %>%
  # calculate stability from function/m2
  group_by(plot) %>%
  summarize(stability = 1/var(total_function)/mean(total_function),
            richness = mean(richness))

# calculate response diversity
response_diversity <- biomass_clim %>%
  # exclude monocultures
  filter(num_sp > 1) %>%
  # get mean function for each species each year for each plot
  group_by(year, month, plot, species) %>%
  summarize(mean_biomass_g_m2 = mean(biomass_g_m2),
            richness = mean(sp_num)) %>%
  # get mean function across all years for each species for each plot
  group_by(plot, species) %>%
  summarize(mean_biomass_g_m2 = mean(mean_biomass_g_m2), richness = mean(richness)) %>%
  # bring in optima for each species
  left_join(fits_summary, by = "species") %>%
  # calculate response diversity as weighted variance in optima, where weights are mean biomass of each species
  group_by(plot) %>%
  summarize(response_diversity = Hmisc::wtd.var(max_temp_optimum, mean_biomass_g_m2)/mean(max_temp_optimum, na.rm=TRUE),
            total_function = sum(mean_biomass_g_m2))

# bring together stability and RD for each plot
stab_rd <- stability %>%
  left_join(response_diversity, by = "plot") %>%
  # bring in soil traits
  left_join(soil_traits_norm, by = "plot") %>%
  # normalize rd and stability
  mutate(stability_norm = normalize(stability, method = "range") + 0.0000000001,
         rd_norm = normalize(response_diversity, method = "range"),
         richness_norm = normalize(richness, method = "range"),
         tf_norm = normalize(total_function, method = "range")) %>%
  drop_na()


# visualize data
ggplot(data = stab_rd, aes(x = response_diversity, y = log(stability), col = total_function, size = richness)) +
  geom_point(alpha=0.5) +
  #geom_smooth(method = "lm") +
  scale_color_viridis_c() +
  labs(x = "Response Diversity (Temperature)", y = "log(Functional Stability)", col = "Total Function (g/m^2)", size = "Richness") +
  theme_minimal()

# build models
model <- lm(log(stability) ~ response_diversity*total_function
            + richness
            + phosphorus
            + no2no3
            + nh4
            + nitrogen,
            data = stab_rd)
summary(model)

model_norm <- lm(log(stability_norm) ~ rd_norm*tf_norm
                 + richness_norm
                 + phosphorus
                 + no2no3
                 + nh4
                 + nitrogen,
                 data = stab_rd)
summary(model_norm)

simple_model <- lm(log(stability_norm) ~ rd_norm + richness_norm,
                   data = stab_rd)
summary(simple_model)

partials <- plot_model(model, type = "pred")
plot_grid(partials[[1]], partials[[2]], partials[[7]])

# test for heteroskedasticity
library(lmtest)

# run Breusch-Pagan test, p < 0.05 means there's heteroskedasticity
bptest(simple_model)

# calculate more robust standard erros
library(sandwich)

# robust standard errors (Huber-White)
coeftest(model_norm, vcov = vcovHC(model_norm, type = "HC1"))

# check for unobserved confounders
library(sensemakr)

sense <- sensemakr(model = model_norm, 
                   treatment = "rd_norm", 
                   benchmark_covariates = NULL, 
                   kd = 1)

summary(sense)
plot(sense)


# evaluate effect of environment on performance
biomass_clim_summary <- biomass_clim %>%
  group_by(year, month, plot, species) %>%
  summarize(mean_biomass_g_m2 = mean(biomass_g_m2),
            avg_max_temp_c = mean(mean_max_temp),
            avg_precip_mm = mean(mean_precip),
            richness = mean(num_sp)) %>%
  group_by(year, month, plot) %>%
  summarize(total_function = sum(mean_biomass_g_m2),
            avg_max_temp_c = mean(avg_max_temp_c),
            avg_precip_mm = mean(avg_precip_mm),
            richness = mean(richness))

# visualize
ggplot(data = biomass_clim_summary, aes(x = avg_max_temp_c, y = total_function, col = richness)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(data = biomass_clim_summary, aes(x = avg_precip_mm, y = total_function, col = richness)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(data = biomass_clim_summary, aes(x = richness, y = total_function, col = year)) +
  geom_point() +
  geom_smooth(method = "lm")



e_model <- lm(total_function ~ avg_precip_mm + avg_max_temp_c + richness, data = biomass_clim_summary)
summary(e_model)

# partial effects plots
partial_effects <- plot_model(model = e_model, type = "pred")
cowplot::plot_grid(partial_effects[[1]], partial_effects[[2]], partial_effects[[3]], ncol = 3, nrow = 1)

# check if environment drive function of each species
each_species_summary <- biomass_clim %>%
  group_by(year, month, plot, species) %>%
  summarize(richness = mean(num_sp),
            mean_biomass_g_m2 = mean(biomass_g_m2),
            avg_max_temp_c = mean(mean_max_temp),
            avg_precip_mm = mean(mean_precip))

# normalize values
each_species_norm <- each_species_summary %>%
  ungroup() %>%
  mutate(mean_biomass_g_m2 = (mean_biomass_g_m2 - min(mean_biomass_g_m2)) / (max(mean_biomass_g_m2) - min(mean_biomass_g_m2)),
         avg_max_temp_c = (avg_max_temp_c - min(avg_max_temp_c)) / (max(avg_max_temp_c) - min(avg_max_temp_c)),
         avg_precip_mm = (avg_precip_mm - min(avg_precip_mm)) / (max(avg_precip_mm) - min(avg_precip_mm)),
         richness = (richness - min(richness)) / (max(richness) - min(richness)))

# visualize
ggplot(data = each_species_summary, aes(x = avg_max_temp_c, y = mean_biomass_g_m2, col = species, size = richness)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  theme(legend.position = "none")

biomass_temp_model <- lmer(mean_biomass_g_m2 ~ avg_max_temp_c + avg_precip_mm + richness + (1 + avg_max_temp_c + avg_precip_mm + richness | species),
                           data = each_species_summary)
summary(biomass_temp_model)

model_fe <- feols(
  mean_biomass_g_m2 ~ avg_max_temp_c + avg_precip_mm + richness | species,
  cluster = ~ species,
  data = each_species_summary
)

summary(model_fe)





  
