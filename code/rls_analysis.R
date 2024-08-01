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

# discard outliers

# fit abundance functions for every species of fish

# get lower and upper bounds for environment vector
lower <- min(rls_temp$sst)
upper <- max(rls_temp$sst)

# Define the Gaussian function
gaussian <- function(x, a, b, c) {
  a * exp(-((x - b)^2) / (2 * c^2))
}

gaussian_fits <- bind_rows(lapply(unique(rls_temp$species_name), function(x){
  #print(x)
  curr_species <- rls_temp %>%
    filter(species_name ==  x) %>%
    ungroup()
  
  # check if there's more than 5 observations
  if(nrow(curr_species) >= 5){
    # drop outliers
    # Calculate Q1 and Q3
    Q1 <- quantile(curr_species$total_biomass, 0.25)
    Q3 <- quantile(curr_species$total_biomass, 0.75)
    
    # Compute the IQR
    IQR <- Q3 - Q1
    
    # Define the bounds
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    
    # Filter the vector to remove outliers
    curr_species_filtered <- curr_species %>%
      filter(total_biomass >= lower_bound & total_biomass <= upper_bound)
    
    # Initial parameter estimates
    start_params <- list(a = Q3, b = weighted.mean(curr_species_filtered$sst, curr_species_filtered$total_biomass), c = sd(curr_species_filtered$sst))
    
    fit_result <- tryCatch(
      {
        fit <- minpack.lm::nlsLM(
          total_biomass ~ gaussian(sst, a, b, c), 
          data = curr_species_filtered, 
          start = start_params,
          lower = c(-Inf, lower, -Inf),
          upper = c(Inf, upper, Inf),
          control = nls.control(maxiter = 1000, minFactor = 1/1024)
        )
        list(success = TRUE, fit = fit)
      },
      error = function(e) {
        warning("An error occurred during model fitting: ", e$message)
        list(success = FALSE, error = e$message)
      }
    )
    
    if (!fit_result$success) {
      # Handle the error as needed
      print(fit_result$error)
      fitted_curve <- data.frame(
        species_name = NULL,
        sst = NULL,
        predicted = NULL
      )
    } else {
      # Proceed with the fit object
      fit <- fit_result$fit
      # Generate fitted values for the original data points
      curr_species_filtered$predicted <- predict(fit, newdata = curr_species_filtered)
      
      # Generate fitted values for a smooth curve
      temp_seq <- seq(lower, upper, length.out = 100)
      fitted_curve <- data.frame(
        species_name = x,
        sst = temp_seq,
        predicted = predict(fit, newdata = data.frame(sst = temp_seq))
      )
    }
  } else {
    fitted_curve <- data.frame(
      species_name = NULL,
      sst = NULL,
      predicted = NULL
    )
  }
}))

# plot abundance for a couple fish
# get all fish
all_fish_sp <- unique(gaussian_fits$species_name)
# pick some random fish to plot
rand_fish <- sample(all_fish_sp, 5)

curr_rls <- rls_temp %>%
  filter(species_name %in% rand_fish)

curr_fits <- gaussian_fits %>%
  filter(species_name %in% rand_fish)

ggplot(data = curr_rls, aes(x = sst, y = total_biomass)) +
  geom_point(alpha = 0.5) +
  geom_line(data = curr_fits, aes(x = sst, y = predicted), col = "red") +
  labs(x = "Sea Surface Temperature (C)", y = "Biomass (g)") +
  facet_grid(~species_name) +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())


ggplot(data = rls_temp, aes(x = sst, y = total_biomass)) +
  geom_point() +
  geom_line(data = gaussian_fits, aes(x = sst, y = predicted)) +
  #stat_function(fun = dnorm, args = list(mean = mean(curr_rls$sst), sd = sd(curr_rls$sst))) +
  labs(x = "Sea Surface Temperature (C)", y = "Biomass (g)") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())

# plot a fish
ggplot(data = filter(rls_temp, species_name == "Lutjanus lutjanus"), aes(x = sst, y = total_biomass)) +
  geom_point() +
  #geom_line(data = fitted_curve, aes(x = sst, y = predicted)) +
  #stat_function(fun = dnorm, args = list(mean = mean(curr_rls$sst), sd = sd(curr_rls$sst))) +
  labs(x = "Sea Surface Temperature (C)", y = "Biomass (g)") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())



