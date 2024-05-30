# load libraries
library(tidyverse)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# define parameters
n_species <- 10
abundance_intercept_range <- c(1000,1000)
abundance_slope_range <- c(-10, 10)
function_intercept_range <- c(0, 0)
function_slope_range <- c(0, 10)
environment_vals <- seq(0, 100, 1 )
n_simulations <- 500

# generate random abundance functions
generate_random_functions <- function(n_species, abund_intercept_range, abund_slope_range, funct_intercept_range, funct_slope_range){
  # loop over each species in model and generate random linear abundance and "function" functions for each one
  # model:
    # abundance = slope * environment + intercept
    # function = slope * abundance + intercept
  do.call(rbind, lapply(1:n_species, function(i){
    data.frame(species_ID = i,
               abundance_intercept = runif(1, abund_intercept_range[1], abund_intercept_range[2]),
               abundance_slope = runif(1, abund_slope_range[1], abund_slope_range[2]),
               function_intercept = runif(1, funct_intercept_range[1], funct_intercept_range[2]),
               function_slope = runif(1, funct_slope_range[1], funct_slope_range[2]))
  }))
}

test <- generate_random_functions(n_species, abundance_intercept_range, abundance_slope_range, function_intercept_range, function_slope_range)

# function to run one simulation (one set of random models)
run_one_sim <- function(n_species, environment_vals, abundance_intercept_range, abundance_slope_range, function_intercept_range, function_slope_range){
  # generate a random set of models
  models <- generate_random_functions(n_species,
                                        abundance_intercept_range,
                                        abundance_slope_range,
                                        function_intercept_range,
                                        function_slope_range)
  
  # calculate abundances and functions for a range of environment values
  model_results <- do.call(rbind, lapply(environment_vals, function(E){
    # loop over every species in model
    do.call(rbind, lapply(1:nrow(models), function(i){
      curr_abundance <- ifelse(((models$abundance_slope[i] * E) + models$abundance_intercept[i]) > 0, ((models$abundance_slope[i] * E) + models$abundance_intercept[i]), 0)
      curr_function <- (models$function_slope[i] * curr_abundance) + models$function_intercept[i]
      data.frame(E = E,
                 species_ID = models$species_ID[i],
                 slope = models$abundance_slope[i],
                 abundance = curr_abundance,
                 funct = curr_function)
    }))
  }))
}

# function to run n simulations
run_n_sims <- function(n_sims, n_species, environment_vals, abundance_intercept_range, abundance_slope_range, function_intercept_range, function_slope_range){
  do.call(rbind, lapply(1:n_sims, function(n){
    curr_results <- run_one_sim(n_species, environment_vals, abundance_intercept_range, abundance_slope_range, function_intercept_range, function_slope_range) %>%
      mutate(sim_number = n)
  }))
}

# run n sets of models
all_results <- run_n_sims(n_simulations, n_species, environment_vals, abundance_intercept_range, abundance_slope_range, function_intercept_range, function_slope_range)

# calculate resilience for each set of models (as inverse variance in sum of functions across environmental gradient)
resilience <- all_results %>%
  group_by(sim_number, E) %>%
  summarize(total_function = sum(funct)) %>%
  group_by(sim_number) %>%
  summarize(resilience = 1/var(total_function),
            total_function = sum(total_function))

# calculate response diversity as variance in slopes for abundance as a function of environment
response_diversity <- all_results %>%
  group_by(sim_number, species_ID) %>%
  summarize(slope = median(slope),
            abs_slope = abs(slope)) %>%
  #ungroup() %>%
  group_by(sim_number) %>%
  summarize(response_diversity = var(slope),
            sum_sensitivity = sum(abs_slope),
            mean_sensitivity = mean(abs_slope))

# join together resilience and response diversity data
model_stats <- resilience %>%
  left_join(response_diversity, by = "sim_number")

# plot model results

# Fig 1, abundance vs environment
ggplot(data = filter(all_results, sim_number == 1), aes(x = E, y = abundance, col = as.factor(species_ID))) +
  geom_line(size = 3) +
  #geom_smooth(method = "lm") +
  labs(x = "E", y = "Abundance", color = "Species") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))

# Fig 2, function vs abundance
ggplot(data = filter(all_results, sim_number == 1), aes(x = abundance, y = funct, col = as.factor(species_ID))) +
  geom_line(size = 3) +
  #geom_smooth(method = "lm") +
  labs(x = "Abundance", y = "Function Provisioning", color = "Species") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))

# Fig 3, function vs environment
ggplot(data = filter(all_results, sim_number == 1), aes(x = E, y = funct, col = as.factor(species_ID))) +
  geom_line(size = 3) +
  #geom_point() +
  labs(x = "E", y = "Ecosystem Function", color = "Species") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))

# Fig 4, resilience vs response diversity
ggplot(data = model_stats, aes(x = response_diversity, y = log(resilience))) +
  geom_point(size=3) +
  geom_smooth(method = "lm") +
  labs(x = "Response Diversity", y = "log(Resilience)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))

ggplot(data = model_stats, aes(x = response_diversity, y = total_function)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

ggplot(data = model_stats, aes(x = sum_sensitivity, y = log(resilience))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

ggplot(data = model_stats, aes(x = mean_sensitivity, y = log(resilience))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

all_results %>%
  group_by(sim_number, E) %>%
  summarize(total_function = sum(funct)) %>%
  ungroup() %>%
  filter(sim_number %in% 1:50) %>%
  ggplot(aes(x = E, y = total_function, col = as.factor(sim_number))) +
  geom_line() +
  labs(x = "E", y = "Total Function", col = "Ecosystem ID") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))

# do some statistical analysis
resilience_response_model <- lm(log(resilience) ~ response_diversity, data = model_stats)
null_model <- lm(log(resilience) ~ 1, data = model_stats)

summary(resilience_response_model)
summary(null_model)

AIC(resilience_response_model)
AIC(null_model)

report::report(resilience_response_model)

