# set working directory
setwd("/projects/lude8513/response_diversity")
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd("../")

# load libraries
library(parallel)

# load helper functions
source("code/rd_utils.R")

# define default parameter ranges
# what proportion of species contribute to function
proportion_contribute <- seq(0.1, 1, 0.1)
# linear abundance function
abundance_slope_range <- c(-1, 1)
# gaussian abundance function
b_range <- c(0, 100)
# function function
function_intercept_range <- c(0, 0)
function_slope_range <- c(0, 10)
# standard deviation range to sample from
sd_range <- c(0.1, 30)
# environmental gradient vector
environment_vals <- seq(0, 100, 1)
# how many ecosystems to generate
n_simulations <- 1000


# generate data

# define ranges for parameters
# format: 
#  linear - list(function_intercept_range, function_slope_range, abundance_intercept_range, abundance_slope_range)
#  gaussian - list(function_intercept_range, function_slope_range, a_range, b_range, c_range)

# define parameter ranges for different scenarios
linear_small_int_ranges <- list(function_intercept_range, function_slope_range, c(5, 5), abundance_slope_range, sd_range)
linear_mid_int_ranges <- list(function_intercept_range, function_slope_range, c(50, 50), abundance_slope_range, sd_range)
linear_large_int_ranges <- list(function_intercept_range, function_slope_range, c(500, 500), abundance_slope_range, sd_range)
linear_rand_int_ranges <- list(function_intercept_range, function_slope_range, c(-100, 100), abundance_slope_range, sd_range)

gaussian_constant_ranges <- list(function_intercept_range, function_slope_range, c(30, 30), b_range, c(10, 10), sd_range)
gaussian_varied_ranges <- list(function_intercept_range, function_slope_range, c(0, 50), b_range, c(0, 20), sd_range)

p_contribute_ranges <- list(function_intercept_range, function_slope_range, c(0, 50), b_range, c(0, 20), sd_range)

crossing_large_ranges <- list(function_intercept_range, function_slope_range, c(0, 100), c(-1, -1), sd_range)
crossing_small_ranges <- list(function_intercept_range, function_slope_range, c(0, 100), c(-0.1, -0.1), sd_range)
crossing_rand_ranges <- list(function_intercept_range, function_slope_range, c(-100, 100), c(-1, 1), sd_range)

# generate data for proportion contribute sensitivity analysis
p_contribute_n10_results <- do.call(rbind, lapply(proportion_contribute, function(p){
  # run simulations for current value of p
  curr_p <- run_n_sims(n_simulations, n_species = 10, environment_vals, p_contribute_ranges, response_shape = 'gaussian', p_contribute = p)

  #print(paste0(typeof(curr_p), " curr type"))
  #print(head(curr_p))

  print("reached resilience calc")

  # calc resilience
  resilience <- curr_p %>%
    group_by(sim_number, E) %>%
    summarize(total_function = sum(funct)) %>%
    group_by(sim_number) %>%
    summarize(resilience = 1/var(total_function))

  print("post resilience calc")

  # calc weighted response diversity
  mean_weighted_response_diversity <- curr_p %>%
    group_by(sim_number, species_ID) %>%
    summarize(mean_abundance = mean(abundance), a = median(a), b = median(b), c = median(c)) %>%
    group_by(sim_number) %>%
    summarize(mean_weighted_response_diversity = Hmisc::wtd.var(b, mean_abundance))

  curr_stats <- resilience %>%
    left_join(mean_weighted_response_diversity, by = "sim_number") %>%
    mutate(p = p)

  print(paste0("finished sim: p = ", p))
  curr_stats
}))

p_contribute_n50_results <- do.call(rbind, lapply(proportion_contribute, function(p){
  # run simulations for current value of p
  curr_p <- run_n_sims(n_simulations, n_species = 50, environment_vals, p_contribute_ranges, response_shape = 'gaussian', p_contribute = p)

  # calc resilience
  resilience <- curr_p %>%
    group_by(sim_number, E) %>%
    summarize(total_function = sum(funct)) %>%
    group_by(sim_number) %>%
    summarize(resilience = 1/var(total_function))

  # calc weighted response diversity
  mean_weighted_response_diversity <- curr_p %>%
    group_by(sim_number, species_ID) %>%
    summarize(mean_abundance = mean(abundance), a = median(a), b = median(b), c = median(c)) %>%
    group_by(sim_number) %>%
    summarize(mean_weighted_response_diversity = Hmisc::wtd.var(b, mean_abundance))

  curr_stats <- resilience %>%
    left_join(mean_weighted_response_diversity, by = "sim_number") %>%
    mutate(p = p)

  print(paste0("finished sim: p = ", p))
  curr_stats
}))

# run simulations for main analyses
linear_small_int_results <- run_n_sims(n_simulations, n_species = 10, environment_vals, linear_small_int_ranges, 'linear', p_contribute = 1)
linear_mid_int_results <- run_n_sims(n_simulations, n_species = 10, environment_vals, linear_mid_int_ranges, 'linear', p_contribute = 1)
linear_large_int_results <- run_n_sims(n_simulations, n_species = 10, environment_vals, linear_large_int_ranges, 'linear', p_contribute = 1)
linear_rand_int_results <- run_n_sims(n_simulations, n_species = 10, environment_vals, linear_rand_int_ranges, 'linear', p_contribute = 1)
linear_rand_int_n50_results <- run_n_sims(n_simulations, n_species = 50, environment_vals, linear_rand_int_ranges, 'linear', p_contribute = 1)

gaussian_constant_results <- run_n_sims(n_simulations, n_species = 10, environment_vals, gaussian_constant_ranges, 'gaussian', p_contribute = 1)
gaussian_varied_n10_results <- run_n_sims(n_simulations, n_species = 10, environment_vals, gaussian_varied_ranges, 'gaussian', p_contribute = 1)
gaussian_varied_n50_results <- run_n_sims(n_simulations, n_species = 50, environment_vals, gaussian_varied_ranges, 'gaussian', p_contribute = 1)

crossing_large_slope <- run_n_sims(n_simulations, n_species = 10, environment_vals, crossing_large_ranges, 'linear', p_contribute = 1, perfectly_crossing = TRUE)
crossing_small_slope <- run_n_sims(n_simulations, n_species = 10, environment_vals, crossing_small_ranges, 'linear', p_contribute = 1, perfectly_crossing = TRUE)
crossing_rand_slope <- run_n_sims(n_simulations, n_species = 10, environment_vals, crossing_rand_ranges, 'linear', p_contribute = 1, perfectly_crossing = TRUE)

# write data to CSV's
write_csv(p_contribute_n10_results, "data/sim_data/p_contribute_n10.csv")
write_csv(p_contribute_n50_results, "data/sim_data/p_contribute_n50.csv")

write_csv(linear_small_int_results, "data/sim_data/linear_small_int_results.csv")
write_csv(linear_mid_int_results, "data/sim_data/linear_mid_int_results.csv")
write_csv(linear_large_int_results, "data/sim_data/linear_large_int_results.csv")
write_csv(linear_rand_int_results, "data/sim_data/linear_rand_int_results.csv")
write_csv(linear_rand_int_n50_results, "data/sim_data/linear_rand_int_n50_results.csv")

write_csv(gaussian_constant_results, "data/sim_data/gaussian_constant_results.csv")
write_csv(gaussian_varied_n10_results, "data/sim_data/gaussian_varied_n10_results.csv")
write_csv(gaussian_varied_n50_results, "data/sim_data/gaussian_varied_n50_results.csv")

write_csv(crossing_large_slope, "data/sim_data/crossing_large_slope_results.csv")
write_csv(crossing_small_slope, "data/sim_data/crossing_small_slope_results.csv")
write_csv(crossing_rand_slope, "data/sim_data/crossing_rand_slope_results.csv")

