# load libraries -------
library(tidyverse)
library(vegan)
library(cluster)
library(cowplot)
library(missForest)

# set working directory --------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load helper functions ---------
source("rd_utils.R")

# load data -----------
gaussian_varied_n10_results <- read_csv("../data/big_sim_data/gaussian_varied_n10_results.csv")

# simulate missing data ---------
p_missing <- 0.8

# figure out which indices to remove
set.seed(123)  # Set a seed for reproducibility
n <- nrow(gaussian_varied_n10_results)
ind_to_remove <- sample(seq_len(n), size = floor(n * p_missing))

# remove data at specified indices and fill in NA's
missing_data <- gaussian_varied_n10_results %>%
  mutate(abundance = ifelse(row_number() %in% ind_to_remove, NA, abundance),
         funct = ifelse(row_number() %in% ind_to_remove, NA, funct))

# shannon RD ----------
# calculate relative abundances
relative_abundance <- gaussian_varied_n10_results %>%
  select(c(E, species_ID, sim_number, abundance)) %>%
  pivot_wider(names_from = species_ID, values_from = abundance) %>%
  mutate(across(
    -c(E, sim_number), function(i){i/rowSums(.[3:12])}
  ))

# calculate shannon diversity/entropy for each row
shannon_diversity <- diversity(relative_abundance[3:12], index = "shannon")

# calculate mean diversity/entropy across all environmental conditions
shannon_entropy <- relative_abundance %>%
  cbind(shannon_diversity) %>%
  group_by(sim_number) %>%
  summarise(mean_entropy = mean(shannon_diversity))

# euclidean RD -----------
# iterate over ecosystems and calculate euclidean RD for each one
euclidean_RD <- do.call(rbind, lapply(unique(relative_abundance$sim_number), function(i){
  # filter out data for current ecosystem
  curr_sim <- filter(relative_abundance, sim_number == i)
  
  # convert relative abundance data to matrix
  rel_abundance_mat <- as.matrix(curr_sim[,3:12])
  
  # transpose matrix
  rel_abundance_t <- t(rel_abundance_mat)
  
  # calculate distance matrix
  dist_mat <- dist(rel_abundance_t)
  
  # calculate mean distance
  mean_dist <- mean(dist_mat)
  
  data.frame(sim_number = i, e_response_diversity = mean_dist)
}))

# gower RD ----------------
# iterate over ecosystems and calculate euclidean RD for each one
gower_RD <- do.call(rbind, lapply(unique(relative_abundance$sim_number), function(i){
  # filter out data for current ecosystem
  curr_sim <- filter(relative_abundance, sim_number == i)
  
  # convert relative abundance data to matrix
  rel_abundance_mat <- as.matrix(curr_sim[,3:12])
  
  # transpose matrix
  rel_abundance_t <- t(rel_abundance_mat)
  
  # calculate distance matrix
  dist_mat <- daisy(rel_abundance_t, metric = "gower")
  
  # calculate mean distance
  mean_dist <- mean(dist_mat)
  
  data.frame(sim_number = i, g_response_diversity = mean_dist)
}))

# join in stats
# calculate stats
stats <- calc_stats(gaussian_varied_n10_results, "gaussian") %>%
  left_join(shannon_entropy, by = "sim_number") %>%
  left_join(euclidean_RD, by = "sim_number") %>%
  left_join(gower_RD, by = "sim_number")

# make some plots
ggplot(data = stats, aes(x = w_response_diversity, y = mean_entropy)) +
  geom_point()

ggplot(data = stats, aes(x = w_response_diversity, y = e_response_diversity)) +
  geom_point()

ggplot(data = stats, aes(x = mean_entropy, y = log(resilience))) +
  geom_point()

# compare 4 metrics for RD

p1 <- ggplot(data = stats, aes(x = w_response_diversity, y = log(resilience))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = "", y = "log(Resilience)") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())

p2 <- ggplot(data = stats, aes(x = uw_response_diversity, y = log(resilience))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())

p3 <- ggplot(data = stats, aes(x = mean_entropy, y = log(resilience))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = "Response Diversity", y = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())

p4 <- ggplot(data = stats, aes(x = e_response_diversity, y = log(resilience))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = "Response Diversity", y = "log(Resilience)") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())

p5 <- ggplot(data = stats, aes(x = g_response_diversity, y = log(resilience))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = "Response Diversity", y = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        panel.grid = element_blank())


plot_grid(p1, p2, p3, p4, p5,
          labels = c("Weighted Variance", "Unweighted Variance", "Shannon Entropy", "Euclidean Distance", "Gower Distance"),
          label_size = 20,
          hjust = -0.6,
          nrow = 2)







