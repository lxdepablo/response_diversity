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

# Fig S1. -----------
p_contribute_plot(p_contribute_n10_stats)
p_contribute_plot(p_contribute_n50_stats)

# Fig S?, abundance vs environment ------------
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

# Fig S?, function vs abundance ----------
function_abundance_plot(linear_small_int_results)
function_abundance_plot(linear_large_int_results)
function_abundance_plot(linear_rand_int_results)
function_abundance_plot(gaussian_varied_n10_results)
function_abundance_plot(gaussian_varied_n50_results)
function_abundance_plot(crossing_small_slope_results)

# Fig S?, total function vs environment -----
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

plot_grid(tf1, tf2, tf3, labels = c("Linear", "Gaussian", "Crossing"),
          nrow = 1, hjust = -1.7, vjust = 2, label_size = 20, rel_widths = c(3, 3, 4.5))

# Fig 2, resilience vs mean weighted RD -------
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

# Fig ? -----
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

# Fig ?, risk-reward plot ----
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
summary(lm(log(resilience) ~ w_response_diversity + tolerance, data = linear_small_int_stats))
summary(lm(log(resilience) ~ w_response_diversity + tolerance, data = linear_mid_int_stats))
summary(lm(log(resilience) ~ w_response_diversity + tolerance, data = linear_large_int_stats))
summary(lm(log(resilience) ~ w_response_diversity + tolerance, data = linear_rand_int_stats))
summary(lm(log(resilience) ~ w_response_diversity + a*c, data = gaussian_varied_n10_stats))
summary(lm(log(resilience) ~ w_response_diversity + a*c, data = gaussian_varied_n50_stats))
summary(lm(log(resilience) ~ w_response_diversity + a*c + coverage_prop, data = gaussian_constant_stats))

# this model includes things related to function (slope and sd)
# including these covariates brings up the R2 quite a bit
total_model <- lm(log(resilience) ~ w_response_diversity
           + a * c
           + mean_f_slope
           + sd_f_slope
           + evenness,
           data = gaussian_varied_n10_stats)
summary(total_model)
# assess multicollinearity
car::vif(total_model, type = "predictor")

# these models incorporate metrics for evenness, coverage, and redundancy
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

## total_function ~ response diversity ------------
summary(lm(total_function ~ uw_response_diversity, data = linear_small_int_stats))
summary(lm(total_function ~ uw_response_diversity, data = linear_large_int_stats))
summary(lm(total_function ~ uw_response_diversity, data = linear_rand_int_stats))
summary(lm(total_function ~ uw_response_diversity, data = gaussian_varied_n10_stats))
summary(lm(total_function ~ uw_response_diversity, data = gaussian_varied_n50_stats))

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
