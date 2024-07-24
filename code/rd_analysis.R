# load libraries -------
library(tidyverse)
library(cowplot)

# set working directory --------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load helper functions ---------
source("rd_utils.R")

# load data ------------
p_contribute_n10_results <- read_csv("../data/big_sim_data/p_contribute_n10.csv")
p_contribute_n50_results <- read_csv("../data/big_sim_data/p_contribute_n50.csv")

linear_small_int_results <- read_csv("../data/big_sim_data/linear_small_int_results.csv")
linear_large_int_results <- read_csv("../data/big_sim_data/linear_large_int_results.csv")
linear_rand_int_results <- read_csv("../data/big_sim_data/linear_rand_int_results.csv")

gaussian_constant_results <- read_csv("../data/big_sim_data/gaussian_constant_results.csv")
gaussian_varied_n10_results <- read_csv("../data/big_sim_data/gaussian_varied_n10_results.csv")
gaussian_varied_n50_results <- read_csv("../data/big_sim_data/gaussian_varied_n50_results.csv")

# calculate stats ----------

## proportion contribute sensitivity analysis ---------
# get all values for p
p_vec <- unique(p_contribute_n10_results$p)
# generate models for each value and pull out model stats
p_contribute_n10_stats <- do.call(rbind, lapply(p_vec, function(curr_p){
  # subset data for current value of p
  curr_p_data <- filter(p_contribute_n10_results, p == curr_p) %>%
    filter(resilience != Inf)
  
  # build model
  curr_model <- lm(log(resilience) ~ mean_weighted_response_diversity, data = curr_p_data)
  
  conf_int <- confint(curr_model, 'mean_weighted_response_diversity', level=0.95)

  p_val <- summary(curr_model)$coefficients[,4][[2]]
  r_squared <- summary(curr_model)$r.squared
  estimate <- summary(curr_model)$coefficients[,1][[2]]
  
  data.frame(proportion_contribute = curr_p, p_val = p_val, r_squared = r_squared, estimate = estimate, conf_lower = conf_int[[1]], conf_upper = conf_int[[2]])
}))

p_contribute_n50_stats <- do.call(rbind, lapply(p_vec, function(curr_p){
  # subset data for current value of p
  curr_p_data <- filter(p_contribute_n50_results, p == curr_p) %>%
    filter(resilience != Inf)
  
  # build model
  curr_model <- lm(log(resilience) ~ mean_weighted_response_diversity, data = curr_p_data)
  
  conf_int <- confint(curr_model, 'mean_weighted_response_diversity', level=0.95)
  
  p_val <- summary(curr_model)$coefficients[,4][[2]]
  r_squared <- summary(curr_model)$r.squared
  estimate <- summary(curr_model)$coefficients[,1][[2]]
  
  data.frame(proportion_contribute = curr_p, p_val = p_val, r_squared = r_squared, estimate = estimate, conf_lower = conf_int[[1]], conf_upper = conf_int[[2]])
}))

# linear response shape ---------
# intercept constant and small
linear_small_int_stats <- calc_stats(linear_small_int_results, 'linear')
# intercept constant and large
linear_large_int_stats <- calc_stats(linear_large_int_results, 'linear')
# intercept random and varied
linear_rand_int_stats <- calc_stats(linear_rand_int_results, 'linear')

# gaussian response shape -----------
# SD and peak abundance constant
gaussian_constant_stats <- calc_stats(gaussian_constant_results, 'gaussian')
# SD and peak abundance varied
# n species = 10
gaussian_varied_n10_stats <- calc_stats(gaussian_varied_n10_results, 'gaussian')
# n species = 50
# weighted RD
gaussian_varied_n50_stats <- calc_stats(gaussian_varied_n50_results, 'gaussian')

# species removal experiment -----------
# pick a random value for environment
e_val <- sample(unique(gaussian_varied_n10_results$E), 1)
# remove the rarest species from each ecosystem
rarest_removed <- do.call(rbind, lapply(unique(gaussian_varied_n10_results$sim_number), function(n){
  # filter for current ecosystem
  curr_sim <- gaussian_varied_n10_results %>%
    filter(sim_number == n)
  
  # identify rarest species at given value of E
  curr_E <- curr_sim %>%
    filter(E == e_val)
  
  rarest_sp <- curr_E[which.min(curr_E$abundance), "species_ID"][[1]]
  
  # remove rarest species
  rarest_removed <- curr_sim %>%
    filter(species_ID != rarest_sp)
}))
# remove the most common species from each ecosystem
common_removed <- do.call(rbind, lapply(unique(gaussian_varied_n10_results$sim_number), function(n){
  # filter for current ecosystem
  curr_sim <- gaussian_varied_n10_results %>%
    filter(sim_number == n)
  
  # identify rarest species at given value of E
  curr_E <- curr_sim %>%
    filter(E == e_val)
  
  common_sp <- curr_E[which.max(curr_E$abundance), "species_ID"][[1]]
  
  # remove rarest species
  common_removed <- curr_sim %>%
    filter(species_ID != common_sp)
}))
# remove a random species from each ecosystem
rand_removed <- do.call(rbind, lapply(unique(gaussian_varied_n10_results$sim_number), function(n){
  # filter for current ecosystem
  curr_sim <- gaussian_varied_n10_results %>%
    filter(sim_number == n)
  
  # pick a random species
  rand_sp <- sample(unique(curr_sim$species_ID), 1)
  
  # remove rarest species
  rarest_removed <- curr_sim %>%
    filter(species_ID != rand_sp)
}))

# get pre-removal stats
pre_removal_stats <- gaussian_varied_n10_stats %>%
  rename_with(
    ~ paste0("pr_", .x))

# calculate stats
rarest_removed_stats <- calc_stats(rarest_removed, 'gaussian') %>%
  left_join(pre_removal_stats, by = c("sim_number" = "pr_sim_number")) %>%
  mutate(delta_wrd = w_response_diversity - pr_w_response_diversity,
         delta_urd = w_response_diversity - pr_w_response_diversity,
         delta_resilience = resilience - pr_resilience,
         delta_tf = total_function - pr_total_function,
         removal = "rarest")
common_removed_stats <- calc_stats(common_removed, 'gaussian') %>%
  left_join(pre_removal_stats, by = c("sim_number" = "pr_sim_number")) %>%
  mutate(delta_wrd = w_response_diversity - pr_w_response_diversity,
         delta_urd = w_response_diversity - pr_w_response_diversity,
         delta_resilience = resilience - pr_resilience,
         delta_tf = total_function - pr_total_function,
         removal = "most_common")
rand_removed_stats <- calc_stats(rand_removed, 'gaussian') %>%
  left_join(pre_removal_stats, by = c("sim_number" = "pr_sim_number")) %>%
  mutate(delta_wrd = w_response_diversity - pr_w_response_diversity,
         delta_urd = w_response_diversity - pr_w_response_diversity,
         delta_resilience = resilience - pr_resilience,
         delta_tf = total_function - pr_total_function,
         removal = "random")

species_removal_stats <- rarest_removed_stats %>%
  rbind(common_removed_stats) %>%
  rbind(rand_removed_stats)

# compare treatments
wrd_aov <- aov(delta_wrd ~ removal,
               data = species_removal_stats)
res_aov <- aov(delta_resilience ~ removal,
               data = species_removal_stats)
tf_aov <- aov(delta_tf ~ removal,
              data = species_removal_stats)

summary(wrd_aov)
summary(res_aov)
summary(tf_aov)

# make plots -----------

# species removal comparison
ggplot(data = species_removal_stats, aes(x = removal, y = delta_wrd)) +
  geom_boxplot() +
  labs(x = "Species Removed", y = "\u0394 Weighted Response Diversity") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())
ggplot(data = species_removal_stats, aes(x = removal, y = delta_resilience)) +
  geom_boxplot() +
  labs(x = "Species Removed", y = "\u0394 Resilience") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())
ggplot(data = species_removal_stats, aes(x = removal, y = delta_tf)) +
  geom_boxplot() +
  labs(x = "Species Removed", y = "\u0394 Total Function") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())


# Fig S1.
p_contribute_plot(p_contribute_n10_stats)
p_contribute_plot(p_contribute_n50_stats)

# Fig 1, abundance vs environment
abundance_environment_plot(linear_small_int_results)
abundance_environment_plot(linear_large_int_results)
ae1 <- abundance_environment_plot(linear_rand_int_results) +
  theme(legend.position = "none")
ae2 <- abundance_environment_plot(gaussian_varied_n10_results) +
  labs(y = "")
abundance_environment_plot(gaussian_varied_n50_results)

plot_grid(ae1, ae2, labels = c("Linear", "Gaussian"), hjust = -2, vjust = 2, label_size = 20)

# Fig 2, function vs abundance
function_abundance_plot(linear_small_int_results)
function_abundance_plot(linear_large_int_results)
function_abundance_plot(linear_rand_int_results)
function_abundance_plot(gaussian_varied_n10_results)
function_abundance_plot(gaussian_varied_n50_results)

# Fig 3, total function vs environment
total_function_environment_plot(linear_small_int_results)
total_function_environment_plot(linear_large_int_results)
tf1 <- total_function_environment_plot(linear_rand_int_results) +
  theme(legend.position = "none")
tf2 <- total_function_environment_plot(gaussian_varied_n10_results) +
  labs(y = "")
total_function_environment_plot(gaussian_varied_n50_results)

plot_grid(tf1, tf2, labels = c("Linear", "Gaussian"), hjust = -1.7, vjust = 2, label_size = 20, rel_widths = c(3, 4))

# Fig 4, resilience vs mean weighted response diversity
resilience_w_rd_plot(linear_small_int_stats)
resilience_w_rd_plot(linear_large_int_stats)
resilience_w_rd_plot(linear_rand_int_stats)
resilience_w_rd_plot(gaussian_varied_n10_stats)
resilience_w_rd_plot(gaussian_varied_n50_stats)

# Fig ?, resilience vs unweighted RD
resilience_u_rd_plot(linear_small_int_stats)
resilience_u_rd_plot(linear_large_int_stats)
resilience_u_rd_plot(linear_rand_int_stats)
resilience_u_rd_plot(gaussian_varied_n10_stats)
resilience_u_rd_plot(gaussian_varied_n50_stats)

# Fig ?, total function vs mean weighted RD
function_w_rd_plot(linear_small_int_stats)
function_w_rd_plot(linear_large_int_stats)
function_w_rd_plot(linear_rand_int_stats)
function_w_rd_plot(gaussian_varied_n10_stats)
function_w_rd_plot(gaussian_varied_n50_stats)

# Fig ?, total function vs unweighted RD
function_u_rd_plot(linear_small_int_stats)
function_u_rd_plot(linear_large_int_stats)
function_u_rd_plot(linear_rand_int_stats)
function_u_rd_plot(gaussian_varied_n10_stats)
function_u_rd_plot(gaussian_varied_n50_stats)

# Fig ?, resilience vs mean weighted RD vs total function
rfw1 <- resilience_function_w_rd_plot(linear_small_int_stats) +
  labs(x = "")
rfw2 <- resilience_function_w_rd_plot(linear_large_int_stats) +
  labs(x = "", y = "")
rfw3 <- resilience_function_w_rd_plot(linear_rand_int_stats) +
  labs(y = "")
rfw4 <- resilience_function_w_rd_plot(gaussian_varied_n10_stats)
rfw5 <- resilience_function_w_rd_plot(gaussian_varied_n50_stats) +
  labs(y = "")

plot_grid(rfw1, rfw2, rfw3, rfw4, rfw5, labels = c("Linear Small Intercept",
                                                   "Linear Large Intercept",
                                                   "Linear Random Intercept",
                                                   "Gaussian 10 Species",
                                                   "Gaussian 50 Species"))

# Fig ?, resilience vs unweighted RD vs total function
rfu1 <- resilience_function_u_rd_plot(linear_small_int_stats) +
  labs(x = "")
rfu2 <-resilience_function_u_rd_plot(linear_large_int_stats) +
  labs(x = "", y = "")
rfu3 <-resilience_function_u_rd_plot(linear_rand_int_stats) +
  labs(y = "")
rfu4 <-resilience_function_u_rd_plot(gaussian_varied_n10_stats)
rfu5 <-resilience_function_u_rd_plot(gaussian_varied_n50_stats) +
  labs(y = "")

plot_grid(rfu1, rfu2, rfu3, rfu4, rfu5, labels = c("Linear Small Intercept",
                                                   "Linear Large Intercept",
                                                   "Linear Random Intercept",
                                                   "Gaussian 10 Species",
                                                   "Gaussian 50 Species"))

# make models --------
## log(resilience) ~ WRD ----------
summary(lm(log(resilience) ~ w_response_diversity, data = linear_small_int_stats))
summary(lm(log(resilience) ~ w_response_diversity, data = linear_large_int_stats))
summary(lm(log(resilience) ~ w_response_diversity, data = linear_rand_int_stats))
summary(lm(log(resilience) ~ w_response_diversity, data = gaussian_varied_n10_stats))
summary(lm(log(resilience) ~ w_response_diversity, data = gaussian_varied_n50_stats))

## log(resilience) ~ URD ------------
summary(lm(log(resilience) ~ uw_response_diversity, data = linear_small_int_stats))
summary(lm(log(resilience) ~ uw_response_diversity, data = linear_large_int_stats))
summary(lm(log(resilience) ~ uw_response_diversity, data = linear_rand_int_stats))
summary(lm(log(resilience) ~ uw_response_diversity, data = gaussian_varied_n10_stats))
summary(lm(log(resilience) ~ uw_response_diversity, data = gaussian_varied_n50_stats))

## total_function ~ WRD ------------
summary(lm(total_function ~ uw_response_diversity, data = linear_small_int_stats))
summary(lm(total_function ~ uw_response_diversity, data = linear_large_int_stats))
summary(lm(total_function ~ uw_response_diversity, data = linear_rand_int_stats))
summary(lm(total_function ~ uw_response_diversity, data = gaussian_varied_n10_stats))
summary(lm(total_function ~ uw_response_diversity, data = gaussian_varied_n50_stats))

## total_function ~ URD ------------
summary(lm(total_function ~ w_response_diversity, data = linear_small_int_stats))
summary(lm(total_function ~ w_response_diversity, data = linear_large_int_stats))
summary(lm(total_function ~ w_response_diversity, data = linear_rand_int_stats))
summary(lm(total_function ~ w_response_diversity, data = gaussian_varied_n10_stats))
summary(lm(total_function ~ w_response_diversity, data = gaussian_varied_n50_stats))







