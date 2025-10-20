# load libraries
library(tidyverse)
library(janitor)
library(lme4)
library(fixest)
library(BBmisc)
library(sjPlot)
library(cowplot)
library(parallel)
library(ggeffects)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# get some helper functions
source("cc_utils.R")

# read in regression-ready data
stab_rd <- read_csv("../data/stability_response_diversity.csv")
biomass_clim <- read_csv("../data/biomass_clim.csv")

# visualize data
ggplot(data = stab_rd, aes(x = response_diversity, y = log(stability), col = total_function, size = richness)) +
  geom_point(alpha=0.5) +
  #geom_smooth(method = "lm") +
  scale_color_viridis_c() +
  labs(x = "Response Diversity (Temperature)", y = "log(Functional Stability)", col = "Total Function (g/m^2)", size = "Richness") +
  theme_minimal()

# make plots for fig 3

# this should be a partial regression where richness is held constant
cc_fig3 <- ggplot(data = stab_rd, aes(x = response_diversity, y = log(stability))) +
  geom_point(alpha=0.3) +
  geom_smooth(method = "lm") +
  labs(x = "Response Diversity", y = "log(Functional Stability)", title = "Cedar Creek") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))

# build models
model <- lm(log_stab ~ response_diversity + richness + mean_a * mean_c,
                   data = stab_rd)
summary(model)

partials <- plot_model(model, type = "pred")
plot_grid(partials[[1]], partials[[2]], nrow = 2)
partials[[1]]

# make partial regression with ggplot for pretty figures
pred_rd <- ggpredict(model, terms = "response_diversity")
plot(pred_rd)

# plot marginal effects
ggplot() +
  # points from original data
  geom_point(data = stab_rd, aes(x = rd_norm, y = log_stab), alpha = 0.3) +
  # fitted line
  geom_line(data = pred_rd, aes(x = x, y = predicted), color = "blue", size = 1) +
  # confidence ribbon
  geom_ribbon(data = pred_rd, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(x = "Response Diversity", y = "log(Functional Stability)", title = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 25),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 25))

# test for heteroskedasticity
library(lmtest)

# run Breusch-Pagan test, p < 0.05 means there's heteroskedasticity
bptest(model)

# calculate more robust standard erros
library(sandwich)

# robust standard errors (Huber-White)
coeftest(model, vcov = vcovHC(model, type = "HC1"))

# check for unobserved confounders
library(sensemakr)

sense <- sensemakr(model = model, 
                   treatment = "response_diversity", 
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





  
