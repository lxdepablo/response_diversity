# set working directory
setwd("/Users/lude8513/r_scripts/response_diversity")

# load libraries
library(tictoc)

# load helper functions
source("rd_utils.R")

tic()

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
linear_small_int_ranges <- list(function_intercept_range, function_slope_range, c(5, 5), abundance_slope_range)
linear_large_int_ranges <- list(function_intercept_range, function_slope_range, c(50, 50), abundance_slope_range)
linear_rand_int_ranges <- list(function_intercept_range, function_slope_range, c(0, 100), abundance_slope_range)

gaussian_constant_ranges <- list(function_intercept_range, function_slope_range, c(30, 30), b_range, c(10, 10))
gaussian_varied_ranges <- list(function_intercept_range, function_slope_range, c(0, 50), b_range, c(0, 20))

p_contribute_ranges <- list(function_intercept_range, function_slope_range, c(0, 50), b_range, c(0, 20))

# generate data for proportion contribute sensitivity analysis
p_contribute_n10_results <- do.call(rbind, lapply(proportion_contribute, function(p){
  # run simulations for current value of p
  curr_p <- run_n_sims(n_simulations, n_species = 10, environment_vals, p_contribute_ranges, response_shape = 'gaussian', p_contribute = p)
  
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
linear_large_int_results <- run_n_sims(n_simulations, n_species = 10, environment_vals, linear_large_int_ranges, 'linear', p_contribute = 1)
linear_rand_int_results <- run_n_sims(n_simulations, n_species = 10, environment_vals, linear_rand_int_ranges, 'linear', p_contribute = 1)

gaussian_constant_results <- run_n_sims(n_simulations, n_species = 10, environment_vals, gaussian_constant_ranges, 'gaussian', p_contribute = 1)
gaussian_varied_n10_results <- run_n_sims(n_simulations, n_species = 10, environment_vals, gaussian_varied_ranges, 'gaussian', p_contribute = 1)
gaussian_varied_n50_results <- run_n_sims(n_simulations, n_species = 50, environment_vals, gaussian_varied_ranges, 'gaussian', p_contribute = 1)


# write data to CSV's
write_csv(p_contribute_n10_results, "sim_data/p_contribute_n10.csv")
write_csv(p_contribute_n50_results, "sim_data/p_contribute_n50.csv")

write_csv(linear_small_int_results, "sim_data/linear_small_int_results.csv")
write_csv(linear_large_int_results, "sim_data/linear_large_int_results.csv")
write_csv(linear_rand_int_results, "sim_data/linear_rand_int_results.csv")

write_csv(gaussian_constant_results, "sim_data/gaussian_constant_results.csv")
write_csv(gaussian_varied_n10_results, "sim_data/gaussian_varied_n10_results.csv")
write_csv(gaussian_varied_n50_results, "sim_data/gaussian_varied_n50_results.csv")

toc()

