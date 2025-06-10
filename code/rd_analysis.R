# load libraries -------
library(tidyverse)
library(cowplot)
library(viridis)
library(scales)
library(matrixStats)
library(BBmisc)
library(e1071)

# set working directory --------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load helper functions ---------
source("rd_utils.R")

# load data ------------
p_contribute_n10_results <- read_csv("../data/sim_data/p_contribute_n10.csv")
p_contribute_n50_results <- read_csv("../data/sim_data/p_contribute_n50.csv")

linear_small_int_results <- read_csv("../data/sim_data/linear_small_int_results.csv")
linear_mid_int_results <- read_csv("../data/sim_data/linear_mid_int_results.csv")
linear_large_int_results <- read_csv("../data/sim_data/linear_large_int_results.csv")
linear_rand_int_results <- read_csv("../data/sim_data/linear_rand_int_results.csv")
linear_rand_int_n50_results <- read_csv("../data/sim_data/linear_rand_int_n50_results.csv")

gaussian_constant_results <- read_csv("../data/sim_data/gaussian_constant_results.csv")
gaussian_varied_n10_results <- read_csv("../data/sim_data/gaussian_varied_n10_results.csv")
gaussian_varied_n50_results <- read_csv("../data/sim_data/gaussian_varied_n50_results.csv")

crossing_small_slope_results <- read_csv("../data/sim_data/crossing_small_slope_results.csv")
crossing_large_slope_results <- read_csv("../data/sim_data/crossing_large_slope_results.csv")
crossing_rand_slope_results <- read_csv("../data/sim_data/crossing_rand_slope_results.csv")

# calculate stats ----------

## normalize function
linear_small_int_results <- linear_small_int_results %>%
  mutate(funct = rescale(funct, to = c(0, 1), from = c(0, max(linear_small_int_results$funct))))
linear_mid_int_results <- linear_mid_int_results %>%
  mutate(funct = rescale(funct, to = c(0, 1), from = c(0, max(linear_mid_int_results$funct))))
linear_large_int_results <- linear_large_int_results %>%
  mutate(funct = rescale(funct, to = c(0, 1), from = c(0, max(linear_large_int_results$funct))))
linear_rand_int_results <- linear_rand_int_results %>%
  mutate(funct = rescale(funct, to = c(0, 1), from = c(0, max(linear_rand_int_results$funct))))
linear_rand_int_n50_results <- linear_rand_int_n50_results %>%
  mutate(funct = rescale(funct, to = c(0, 1), from = c(0, max(linear_rand_int_n50_results$funct))))

gaussian_constant_results <- gaussian_constant_results %>%
  mutate(funct = rescale(funct, to = c(0, 1), from = c(0, max(gaussian_constant_results$funct))))
gaussian_varied_n10_results <- gaussian_varied_n10_results %>%
  mutate(funct = rescale(funct, to = c(0, 1), from = c(0, max(gaussian_varied_n10_results$funct))))
gaussian_varied_n50_results <- gaussian_varied_n50_results %>%
  mutate(funct = rescale(funct, to = c(0, 1), from = c(0, max(gaussian_varied_n50_results$funct))))

crossing_small_slope_results <- crossing_small_slope_results %>%
  mutate(funct = rescale(funct, to = c(0, 1), from = c(0, max(crossing_small_slope_results$funct))))
crossing_large_slope_results <- crossing_large_slope_results %>%
  mutate(funct = rescale(funct, to = c(0, 1), from = c(0, max(crossing_large_slope_results$funct))))
crossing_rand_slope_results <- crossing_rand_slope_results %>%
  mutate(funct = rescale(funct, to = c(0, 1), from = c(0, max(crossing_rand_slope_results$funct))))

## proportion contribute sensitivity analysis ---------
# get all values for p
p_vec <- unique(p_contribute_n10_results$p)
# generate models for each value and pull out model stats
p_contribute_n10_stats <- do.call(rbind, lapply(p_vec, function(curr_p){
  # subset data for current value of p
  curr_p_data <- dplyr::filter(p_contribute_n10_results, p == curr_p) %>%
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
linear_small_int_stats <- calc_stats(linear_small_int_results, 'linear') %>%
  filter(w_response_diversity != 0)
# intercept constant and intermediate
linear_mid_int_stats <- calc_stats(linear_mid_int_results, 'linear') %>%
  filter(w_response_diversity != 0)
# intercept constant and large
linear_large_int_stats <- calc_stats(linear_large_int_results, 'linear') %>%
  filter(w_response_diversity != 0)
# intercept random and varied
linear_rand_int_stats <- calc_stats(linear_rand_int_results, 'linear') %>%
  filter(w_response_diversity != 0)
# intercept random and varied (50 species)
linear_rand_int_n50_stats <- calc_stats(linear_rand_int_n50_results, 'linear') %>%
  filter(w_response_diversity != 0)

# gaussian response shape -----------
# SD and peak abundance constant
gaussian_constant_stats <- calc_stats(gaussian_constant_results, 'gaussian') %>%
  filter(w_response_diversity != 0)
# SD and peak abundance varied
# n species = 10
gaussian_varied_n10_stats <- calc_stats(gaussian_varied_n10_results, 'gaussian') %>%
  filter(w_response_diversity != 0)
# n species = 50
# weighted RD
gaussian_varied_n50_stats <- calc_stats(gaussian_varied_n50_results, 'gaussian') %>%
  filter(w_response_diversity != 0)

# perfectly crossing cases ------------
# slope is small and constant
crossing_small_slope_stats <- calc_stats(crossing_small_slope_results, 'linear') %>%
  filter(w_response_diversity != 0, resilience != Inf)

# slope is large and constant
crossing_large_slope_stats <- calc_stats(crossing_large_slope_results, 'linear') %>%
  filter(w_response_diversity != 0, resilience != Inf)

# slope is randomized
crossing_rand_slope_stats <- calc_stats(crossing_rand_slope_results, 'linear') %>%
  filter(w_response_diversity != 0, resilience != Inf)

# calculate summary stats for each case --------------
linear_summary <- data.frame(case = c("small_intercept", "mid_intercept", "large_intercept", "rand_intercept", "rand_intercept_n50"),
                             abund_int_mean = c(mean(linear_small_int_results$abundance_intercept),
                                                mean(linear_mid_int_results$abundance_intercept),
                                                mean(linear_large_int_results$abundance_intercept),
                                                mean(linear_rand_int_results$abundance_intercept),
                                                mean(linear_rand_int_n50_results$abundance_intercept)),
                             abund_int_sd = c(sd(linear_small_int_results$abundance_intercept),
                                              sd(linear_mid_int_results$abundance_intercept),
                                              sd(linear_large_int_results$abundance_intercept),
                                              sd(linear_rand_int_results$abundance_intercept),
                                              sd(linear_rand_int_n50_results$abundance_intercept)))

gaussian_summary <- data.frame(case = c("constant", "varied", "varied_n50"),
                             a_mean = c(mean(gaussian_constant_results$a),
                                                mean(gaussian_varied_n10_results$a),
                                                mean(gaussian_varied_n50_results$a)),
                             a_sd = c(sd(gaussian_constant_results$a),
                                              sd(gaussian_varied_n10_results$a),
                                              sd(gaussian_varied_n50_results$a)),
                             c_mean = c(mean(gaussian_constant_results$c),
                                        mean(gaussian_varied_n10_results$c),
                                        mean(gaussian_varied_n50_results$c)),
                             c_sd = c(sd(gaussian_constant_results$c),
                                      sd(gaussian_varied_n10_results$c),
                                      sd(gaussian_varied_n50_results$c)))

crossing_summary <- data.frame(case = c("small_slope", "large_slope", "rand_slope"),
                             slope_mean = c(mean(crossing_small_slope_results$abundance_slope),
                                                mean(crossing_large_slope_results$abundance_slope),
                                                mean(crossing_rand_slope_results$abundance_slope)),
                             slope_sd = c(sd(crossing_small_slope_results$abundance_slope),
                                              sd(crossing_large_slope_results$abundance_slope),
                                              sd(crossing_rand_slope_results$abundance_slope)))

## collate all cases for boxplots --------
# add case column
linear_cases <- list(linear_small_int_results = linear_small_int_results, 
                     linear_mid_int_results = linear_mid_int_results, 
                     linear_large_int_results = linear_large_int_results, 
                     linear_rand_int_results = linear_rand_int_results, 
                     linear_rand_int_n50_results = linear_rand_int_n50_results)

all_linears <- bind_rows(lapply(names(linear_cases), function(x) {
  df <- linear_cases[[x]]
  df$case <- x
  df
}))

gaussian_cases <- list(gaussian_constant_results = gaussian_constant_results, 
                       gaussian_varied_n10_results = gaussian_varied_n10_results, 
                       gaussian_varied_n50_results = gaussian_varied_n50_results)

all_gaussians <- bind_rows(lapply(names(gaussian_cases), function(x) {
  df <- gaussian_cases[[x]]
  df$case <- x
  df
}))

crossing_cases <- list(crossing_small_slope_results = crossing_small_slope_results, 
                       crossing_large_slope_results = crossing_large_slope_results, 
                       crossing_rand_slope_results = crossing_rand_slope_results)

all_crossings <- bind_rows(lapply(names(crossing_cases), function(x) {
  df <- crossing_cases[[x]]
  df$case <- x
  df
}))

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

# summary stats figs ----------
# linear cases
ggplot(data = all_linears, aes(x = case, y = abundance_intercept)) +
  geom_boxplot()
# gaussian cases
ggplot(data = all_gaussians, aes(x = case, y = a)) +
  geom_boxplot()
ggplot(data = all_gaussians, aes(x = case, y = c)) +
  geom_boxplot()
# crossing cases
ggplot(data = all_crossings, aes(x = case, y = abs(abundance_slope))) +
  geom_boxplot()
ggplot(data = all_crossings, aes(x = case, y = abundance_intercept)) +
  geom_boxplot()

# Fig S1. -----------
p_contribute_plot(p_contribute_n10_stats)
p_contribute_plot(p_contribute_n50_stats)

# Fig 1, abundance vs environment ------------
abundance_environment_plot(linear_small_int_results)
abundance_environment_plot(linear_large_int_results)
abundance_environment_plot(linear_rand_int_results)
abundance_environment_plot(crossing_small_slope_results)
abundance_environment_plot(crossing_large_slope_results)
abundance_environment_plot(crossing_rand_slope_results)

ae1 <- abundance_environment_plot(linear_rand_int_results) +
  theme(legend.position = "none")
ae2 <- abundance_environment_plot(gaussian_varied_n10_results) +
  labs(y = "") +
  theme(legend.position = "none")
ae3 <- abundance_environment_plot(crossing_large_slope_results) +
  labs(y = "")

plot_grid(ae1, ae2, ae3, labels = c("Linear", "Gaussian", "Crossing"),
          nrow = 1,
          hjust = -2,
          vjust = 2,
          rel_widths = c(3, 3, 4),
          label_size = 20)

# Fig 2, function vs abundance ----------
function_abundance_plot(linear_small_int_results)
function_abundance_plot(linear_large_int_results)
function_abundance_plot(linear_rand_int_results)
function_abundance_plot(gaussian_varied_n10_results)
function_abundance_plot(gaussian_varied_n50_results)
function_abundance_plot(crossing_small_slope_results)

# Fig 3, total function vs environment
total_function_environment_plot(linear_small_int_results)
total_function_environment_plot(linear_large_int_results)
total_function_environment_plot(crossing_large_slope_results)
tf1 <- total_function_environment_plot(linear_rand_int_results) +
  theme(legend.position = "none")
tf2 <- total_function_environment_plot(gaussian_varied_n10_results) +
  labs(y = "") +
  theme(legend.position = "none")
tf3 <- total_function_environment_plot(crossing_rand_slope_results) +
  labs(y = "")

plot_grid(tf1, tf2, tf3, labels = c("Linear", "Gaussian", "Crossing"), nrow = 1, hjust = -1.7, vjust = 2, label_size = 20, rel_widths = c(3, 3, 4.5))

# Fig 4, resilience vs mean weighted response diversity
resilience_w_rd_plot(linear_small_int_stats)
resilience_w_rd_plot(linear_large_int_stats)
resilience_w_rd_plot(linear_rand_int_stats)
resilience_w_rd_plot(gaussian_varied_n10_stats)
resilience_w_rd_plot(gaussian_varied_n50_stats)
resilience_w_rd_plot(crossing_small_slope_stats)

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
function_w_rd_plot(crossing_large_slope_stats)

# Fig ?, total function vs unweighted RD
function_u_rd_plot(linear_small_int_stats)
function_u_rd_plot(linear_large_int_stats)
function_u_rd_plot(linear_rand_int_stats)
function_u_rd_plot(gaussian_varied_n10_stats)
function_u_rd_plot(gaussian_varied_n50_stats)

# Fig ?, resilience vs mean weighted RD vs total function
rfw1 <- resilience_w_rd_plot(linear_small_int_stats) +
  labs(y = "", x = "")
rfw2 <- resilience_w_rd_plot(linear_mid_int_stats) +
  labs(y = "", x = "")
rfw3 <- resilience_w_rd_plot(linear_large_int_stats) +
  labs(y = "", x = "")
rfw4 <- resilience_w_rd_plot(linear_rand_int_stats)  +
  labs(y = "", x = "")
rfw5 <- resilience_w_rd_plot(linear_rand_int_n50_stats) +
  labs(y = "", x = "")
rfw6 <- resilience_w_rd_plot(gaussian_constant_stats) +
  labs(y = "", x = "")
rfw7 <- resilience_w_rd_plot(gaussian_varied_n10_stats) +
  labs(y = "", x = "")
rfw8 <- resilience_w_rd_plot(gaussian_varied_n50_stats) +
  labs(y = "", x = "")
rfw9 <- resilience_w_rd_plot(crossing_small_slope_stats) +
  labs(y = "", x = "")
rfw10 <- resilience_w_rd_plot(crossing_large_slope_stats) +
  labs(y = "", x = "")
rfw11 <- resilience_w_rd_plot(crossing_rand_slope_stats) +
  labs(y = "", x = "")

# build full figure for fig 2
linear_plots <- plot_grid(rfw1, rfw2, rfw3, rfw4, rfw5,
                          ncol = 2,
                          hjust = -3.8,
                          labels = c("Case 1",
                                     "Case 4",
                                     "Case 2",
                                     "Case 5",
                                     "Case 3"),
                          byrow = FALSE,
                          align = "hv",
                          axis = "tl")
linear_title <- get_plot_component(ggplot() +
                                     labs(title = "Linear Cases") +
                                     theme(plot.title = element_text(hjust = 0.5,
                                                                     size = 15)),
                                   "title",
                                   return_all=T)[[2]]
linear_plots_title <- plot_grid(linear_title, linear_plots,
                                nrow = 2,
                                rel_heights = c(0.05, 0.95))
gaussian_plots <- plot_grid(rfw6, rfw7, rfw8,
                          ncol = 1,
                          hjust = -3.3,
                          labels = c("Case 6",
                                     "Case 7",
                                     "Case 8"),
                          byrow = FALSE,
                          align = "hv",
                          axis = "tl")
gaussian_title <- get_plot_component(ggplot() +
                                       labs(title = "Gaussian Cases") +
                                       theme(plot.title = element_text(hjust = 0.5,
                                                                       size = 15)),
                                   "title",
                                   return_all=T)[[2]]
gaussian_plots_title <- plot_grid(gaussian_title, gaussian_plots,
                                nrow = 2,
                                rel_heights = c(0.05, 0.95))
crossing_plots <- plot_grid(rfw9,rfw10, rfw11,
                          ncol = 1,
                          hjust = c(-3.1, -3, -3),
                          labels = c("Case 9 ",
                                     "Case 10",
                                     "Case 11"),
                          byrow = FALSE,
                          align = "hv",
                          axis = "tl")
crossing_title <- get_plot_component(ggplot() +
                                       labs(title = "Perfectly Crossing Cases") +
                                       theme(plot.title = element_text(hjust = 0.5,
                                                                       size = 15)),
                                   "title",
                                   return_all=T)[[2]]
crossing_plots_title <- plot_grid(crossing_title, crossing_plots,
                                nrow = 2,
                                rel_heights = c(0.05, 0.95))

fig2_xlab <- get_plot_component(ggplot() +
                                  labs(x = "Response Diversity") +
                                  theme(axis.title = element_text(size = 25)),
                                "xlab-b")
fig2_ylab <- get_plot_component(ggplot() +
                                  labs(y = "log(Functional Stability)") +
                                  theme(axis.title = element_text(size = 25)),
                                "ylab-l")


fig2_plots_ylab <- plot_grid(fig2_ylab,
                             linear_plots_title,
                             gaussian_plots_title,
                             crossing_plots_title,
                             ncol = 4,
                             rel_widths = c(0.04, 0.48, 0.24, 0.24))
fig2_full <- plot_grid(fig2_plots_ylab, fig2_xlab,
                       nrow = 2,
                       rel_heights = c(0.96, 0.04))

# Fig 3
ggplot(data = gaussian_varied_n50_stats, aes(x = normalize(w_response_diversity, method = "range"), y = log(resilience), col = log(tolerance))) +
  geom_point(alpha=0.3, size = 5) +
  geom_smooth(method = "lm") +
  scale_color_viridis_c() +
  labs(x = "Normalized Response Diversity", y = "log(Functional Stability)", title = "Simulated") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))

ggplot(data = gaussian_varied_n10_stats, aes(x = tolerance)) +
  geom_histogram()


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

gaussian_risk_reward <- gaussian_varied_n50_stats %>%
  mutate(rd_group = ceiling(w_response_diversity/200),
         rd_group = ifelse(rd_group >= 4, 4, rd_group),
         rd_group = as.factor(rd_group),
         # convert resilience to risk (1/resilience)
         risk = 1/resilience)

# Get the returns and risk (variance) from the dataset
returns <- gaussian_risk_reward$total_function
risk <- gaussian_risk_reward$risk


# Generate a palette of colors from the Viridis palette for the levels of rd_group
rd_group_colors <- viridis(length(unique(gaussian_risk_reward$rd_group)), option = "C")

ggplot(gaussian_risk_reward, aes(x = (1/resilience), y = total_function)) +
  geom_point(aes(col = w_response_diversity), size = 2, alpha = 0.8) +
  geom_smooth(aes(group = rd_group, col = as.numeric(rd_group)), se = FALSE, span = 2, linewidth = 2, linetype = "dashed") +
  scale_color_viridis_c(name = "Response Diversity") +  # For the points
  labs(x = "Var(Function)", y = "Function") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))


ggplot(linear_rand_int_stats, aes(x = log(1/resilience), y = total_function, col = w_response_diversity)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis_c() +
  labs(x = "Var(Function)", y = "Total Function", col = "Response Diversity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))

# make models --------
## log(resilience) ~ WRD ----------
summary(lm(log(resilience) ~ w_response_diversity + a + c, data = linear_small_int_stats))
summary(lm(log(resilience) ~ w_response_diversity + tolerance, data = linear_large_int_stats))
summary(lm(log(resilience) ~ w_response_diversity + tolerance, data = linear_rand_int_stats))
summary(lm(log(resilience) ~ w_response_diversity + a*c, data = gaussian_varied_n10_stats))
summary(lm(log(resilience) ~ w_response_diversity + a*c, data = gaussian_varied_n50_stats))
summary(lm(log(resilience) ~ w_response_diversity + a*c + coverage_prop, data = gaussian_constant_stats))

total_model <- lm(log(resilience) ~ w_response_diversity
           + a * c
           + mean_f_slope
           + sd_f_slope
           + evenness,
           data = gaussian_varied_n10_stats)
summary(total_model)
car::vif(total_model, type = "predictor")

summary(lm(log(resilience) ~ w_response_diversity +
             a * c +
             tolerance +
             evenness +
             coverage_prop +
             min_spacing,
           data = gaussian_varied_n10_stats))

summary(lm(log(resilience) ~ w_response_diversity +
             a * c +
             redundancy_prop +
             redundancy_count +
             redundancy_score,
           data = gaussian_varied_n50_stats))

summary(lm(log(resilience) ~ w_response_diversity +
             a * c +
             evenness +
             coverage_prop +
             min_spacing +
             skew_spacing +
             redundancy_prop +
             redundancy_score,
           data = gaussian_varied_n50_stats))

summary(lm(log(resilience) ~ w_response_diversity + a * c +
     redundancy_score * c,
   data = gaussian_varied_n50_stats))

summary(lm(log(resilience) ~ w_response_diversity * evenness + a * c, data = gaussian_varied_n50_stats))


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

## total_function ~ stability -------------
summary(lm(total_function ~ log(resilience), data = linear_small_int_stats))
summary(lm(total_function ~ log(resilience), data = linear_large_int_stats))
summary(lm(total_function ~ log(resilience), data = linear_rand_int_stats))
summary(lm(total_function ~ log(resilience), data = gaussian_varied_n10_stats))
summary(lm(total_function ~ log(resilience), data = gaussian_varied_n50_stats))
summary(lm(total_function ~ log(resilience), data = gaussian_constant_stats))

## total_function ~ stability * response_diversity -------------
summary(lm(total_function ~ log(resilience) * w_response_diversity, data = linear_small_int_stats))
summary(lm(total_function ~ log(resilience) * w_response_diversity, data = linear_mid_int_stats))
summary(lm(total_function ~ log(resilience) * w_response_diversity, data = linear_large_int_stats))
summary(lm(total_function ~ log(resilience) * w_response_diversity, data = linear_rand_int_stats))
summary(lm(total_function ~ log(resilience) * w_response_diversity, data = gaussian_varied_n10_stats))
summary(lm(total_function ~ log(resilience) * w_response_diversity, data = gaussian_varied_n50_stats))
summary(lm(total_function ~ log(resilience) * w_response_diversity, data = gaussian_constant_stats))

## stability ~ total_function * response_diversity -------------
summary(lm(log(resilience) ~ total_function * w_response_diversity, data = linear_small_int_stats))
summary(lm(log(resilience) ~ total_function * w_response_diversity, data = linear_mid_int_stats))
summary(lm(log(resilience) ~ total_function * w_response_diversity, data = linear_large_int_stats))
summary(lm(log(resilience) ~ total_function * w_response_diversity, data = linear_rand_int_stats))
summary(lm(log(resilience) ~ total_function * w_response_diversity, data = linear_rand_int_n50_stats))
summary(lm(log(resilience) ~ total_function * w_response_diversity, data = crossing_small_slope_stats))
summary(lm(log(resilience) ~ total_function * w_response_diversity, data = crossing_large_slope_stats))
summary(lm(log(resilience) ~ total_function * w_response_diversity, data = crossing_rand_slope_stats))
summary(lm(log(resilience) ~ total_function * w_response_diversity, data = gaussian_varied_n10_stats))
summary(lm(log(resilience) ~ total_function * w_response_diversity, data = gaussian_varied_n50_stats))
summary(lm(log(resilience) ~ total_function * w_response_diversity, data = gaussian_constant_stats))


## rescale variables to study interaction ------------
gaussian_n10_rescaled <- gaussian_varied_n10_stats %>%
  mutate(total_function_scaled = scale(total_function)[,1],
         log_resilience_scaled = scale(log(resilience))[,1],
         rd_scaled = scale(w_response_diversity)[,1])

summary(lm(log_resilience_scaled ~ total_function_scaled * rd_scaled, data = gaussian_n10_rescaled))


