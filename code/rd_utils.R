# load libraries
if(require(tidyverse)){
  print("loaded tidyverse")
} else {
  print("trying to install tidyverse")
  install.packages("tidyverse")
  if(require(tidyverse)){
    print("tidyverse installed and loaded")
  } else {
    stop("could not install tidyverse")
  }
}

num_cores <- 64

# generate random abundance functions
# inputs:
# n_species: number of species in each ecosystem
# response_shape: shape of abundance function. either 'linear' or 'gaussian'
# params_ranges: vector containing ranges for abundance and function parameters.
#  c(function_intercept, function_slope, linear/gaussian parameters)
# response_shape: "linear" or "gaussian"
# e_max: maximum value for environment

generate_abundance_functions <- function(n_species, params_ranges, response_shape, e_max){
  if(response_shape == "linear") {
    # pull parameter ranges from vector
    funct_intercept_range <- params_ranges[[1]]
    funct_slope_range <- params_ranges[[2]]
    abund_intercept_range <- params_ranges[[3]]
    abund_slope_range <- params_ranges[[4]]
    sd_range <- params_ranges[[5]]
    
    # generate distributions to sample from for current community
    a_int_mean <- runif(1, abund_intercept_range[1], abund_intercept_range[2])
    a_int_sd <- runif(1, sd_range[1], sd_range[2])
    a_slope_mean <- runif(1, abund_slope_range[1], abund_slope_range[2])
    a_slope_sd <- runif(1, sd_range[1], sd_range[2])
    f_int_mean <- runif(1, funct_intercept_range[1], funct_intercept_range[2])
    f_int_sd <- runif(1, sd_range[1], sd_range[2])
    f_slope_mean <- runif(1, funct_slope_range[1], funct_slope_range[2])
    f_slope_sd <- runif(1, sd_range[1], sd_range[2])
    
    # generate abundance and function functions for each species by sampling from normal distributions
    all_species <- data.frame(species_ID = 1:n_species,
               abundance_intercept = rnorm(n_species, a_int_mean),
               abundance_slope = rnorm(n_species, a_slope_mean),
               function_intercept = truncnorm::rtruncnorm(n_species, a = 0, b = Inf, f_int_mean),
               function_slope = truncnorm::rtruncnorm(n_species, a = 0, b = Inf, f_slope_mean))
    
    # check if any species abundances are fully negative, and if so resample them
    all_species_resampled <- bind_rows(lapply(1:n_species, function(i){
      # filter to get current species
      curr_species <- filter(all_species, species_ID == i)
      
      while(TRUE){
        # if intercept > 0, species must have positive values
        if(curr_species$abundance_intercept > 0){
          return(curr_species)
        # check if species is ever positive
        # if intercept and slope both < 0, species will never be positive
        } else if(curr_species$abundance_slope < 0){
          # resample
          curr_species <- data.frame(species_ID = curr_species$species_ID,
                                     abundance_intercept = rnorm(1, a_int_mean),
                                     abundance_slope = rnorm(1, a_slope_mean),
                                     function_intercept = truncnorm::rtruncnorm(1, a = 0, b = Inf, f_int_mean),
                                     function_slope = truncnorm::rtruncnorm(1, a = 0, b = Inf, f_slope_mean))
        # check if highest value of abundance is positive
        } else if(linear_abundance(curr_species$abundance_slope, curr_species$abundance_intercept, e_max) > 0){
          return(curr_species)
        } else {
          # species must be bad
          # resample
          curr_species <- data.frame(species_ID = curr_species$species_ID,
                                     abundance_intercept = rnorm(1, a_int_mean),
                                     abundance_slope = rnorm(1, a_slope_mean),
                                     function_intercept = truncnorm::rtruncnorm(1, a = 0, b = Inf, f_int_mean),
                                     function_slope = truncnorm::rtruncnorm(1, a = 0, b = Inf, f_slope_mean))
        }
      }
    }))
    return(all_species_resampled)
    
  } else if(response_shape == "gaussian"){
    # pull parameter ranges from vector
    funct_intercept_range <- params_ranges[[1]]
    funct_slope_range <- params_ranges[[2]]
    a_range <- params_ranges[[3]]
    b_range <- params_ranges[[4]]
    c_range <- params_ranges[[5]]
    sd_range <- params_ranges[[6]]
    
    # generate distributions to sample from for current community
    a_mean <- runif(1, a_range[1], a_range[2])
    a_sd <- runif(1, sd_range[1], sd_range[2])
    b_mean <- runif(1, b_range[1], b_range[2])
    b_sd <- runif(1, sd_range[1], sd_range[2])
    c_mean <- runif(1, c_range[1], c_range[2])
    c_sd <- runif(1, sd_range[1], sd_range[2])
    f_int_mean <- runif(1, funct_intercept_range[1], funct_intercept_range[2])
    f_int_sd <- runif(1, sd_range[1], sd_range[2])
    f_slope_mean <- runif(1, funct_slope_range[1], funct_slope_range[2])
    f_slope_sd <- runif(1, sd_range[1], sd_range[2])
    
    # generate gaussian abundance functions and linear function functions
    data.frame(species_ID = 1:n_species,
               a = truncnorm::rtruncnorm(n_species, a = 0, b = Inf, a_mean, a_sd),
               b = truncnorm::rtruncnorm(n_species, a = 0, b = Inf, b_mean, b_sd),
               c = truncnorm::rtruncnorm(n_species, a = 0, b = Inf, c_mean, c_sd),
               function_intercept = truncnorm::rtruncnorm(n_species, a = 0, b = Inf, f_int_mean),
               function_slope = truncnorm::rtruncnorm(n_species, a = 0, b = Inf, f_slope_mean))
    
    # make sure a, b, and c are positive
    
  } else {
    print("Response shape must be either 'linear' or 'gaussian'.")
    return(NA)
  }
}

# Gaussian function to calculate abundance as a function of environment
# inputs:
# a: height of the curves peak
# b: position of the center of the peak
# c: standard deviation
# E: value for environment
gaussian_abundance <- function(a, b, c, E){
  abund <- (a*exp(1)^(-((E-b)^2)/(2*(c^2))))
  ifelse(abund > 0, abund, 0)
}

# linear function to calculate abundance as a function of environment
linear_abundance <- function(slope, intercept, E){
  abund <- ((slope*E) + intercept)
  ifelse(abund > 0, abund, 0)
}

# linear function to calculate function/services as a function of abundance (A)
linear_function <- function(slope, intercept, A){
  ((slope * A) + intercept)
}

# function to run one simulation (one set of random models)
# inputs:
# n_species: number of unique species in the ecosystem
# environment_vals: vector of values for environmental gradient/axis
# params_ranges: list of ranges for each species parameter (used in generating abundance functions)
# response_shape: "linear" or "gaussian"
# p_contribute: proportion of species contributing to ecosystem function
#     n_species*p_contribute must be a whole number
# perfectly_crossing: TRUE or FALSE, determines whether abundance functions should be perfectly crossing/mirrored
#     only works if n_species*p_contribute is even and response shape is linear
run_one_sim <- function(n_species, environment_vals, params_ranges, response_shape, p_contribute = 1, perfectly_crossing = FALSE){
  # get maximum environment values
  e_max <- max(environment_vals)
  
  # check if n_species*p_contribute is a whole number, if not throw error
  if(round((n_species*p_contribute)%%1, 3) > 0){
    #print(n_species*p_contribute)
    print((n_species*p_contribute)%%1)
    stop("n_species*p_contribute must be a whole number.")
  }
  
  # generate abundance functions for each species
  if(perfectly_crossing == FALSE){
    # generate a random set of models
    models <- generate_abundance_functions(n_species, params_ranges, response_shape, e_max) %>%
      # decide which species contribute to function
      mutate(function_intercept = ifelse(species_ID <= (p_contribute*n_species), function_intercept, 0),
             function_slope = ifelse(species_ID <= (p_contribute*n_species), function_slope, 0))
    
  } else if(perfectly_crossing == TRUE & round((n_species*p_contribute)%%2, 3) > 0) {
    stop("n_species*p_contribute must be even to generate perfectly crossing abundance functions.")
  } else if(perfectly_crossing == TRUE & response_shape == "gaussian"){
    stop("response shape must be 'linear' to generate perfectly crossing abundance functions.")
  } else {
    # generate a random set of models
    models <- generate_abundance_functions(n_species, params_ranges, response_shape, e_max) %>%
      # decide which species contribute to function
      # use half of desired species (n_species*p_contribute)/2
      mutate(function_intercept = ifelse(species_ID <= (p_contribute*n_species)/2, function_intercept, 0),
             function_slope = ifelse(species_ID <= (p_contribute*n_species)/2, function_slope, 0),
             abundance_intercept = ifelse(species_ID <= (p_contribute*n_species)/2, abundance_intercept, 0),
             abundance_slope = ifelse(species_ID <= (p_contribute*n_species)/2, abundance_slope, 0))
    
    # mirror every species about mean environment
    models <- bind_rows(lapply(1:nrow(models), function(n){
      # find max environment value
      mean_e <- mean(environment_vals)
      
      # if model exists, mirror it
      curr_model <- models[n,]
      if(curr_model$abundance_intercept != 0 & curr_model$abundance_slope != 0){
        new_model <- data.frame(species_ID = curr_model$species_ID + (p_contribute*n_species)/2,
                                abundance_intercept = (2 * curr_model$abundance_slope * mean_e) + curr_model$abundance_intercept,
                                abundance_slope = -1 * curr_model$abundance_slope,
                                function_intercept = curr_model$function_intercept,
                                function_slope = curr_model$function_slope)
        rbind(new_model, curr_model)
      } else {
        return(NULL)
      }
    }))
  }
  #print(paste0("models: ", models))
  
  # calculate abundances and functions for a range of environment values
  model_results <- do.call(rbind, lapply(environment_vals, function(E){
    # loop over every species in model
    do.call(rbind, lapply(1:nrow(models), function(i){
      # calculate abundance for current value of environment
      if(response_shape == 'linear'){
        curr_abundance <- linear_abundance(models$abundance_slope[i], models$abundance_intercept[i], E)
      } else if(response_shape == 'gaussian'){
        curr_abundance <- gaussian_abundance(models$a[i], models$b[i], models$c[i], E)
      } else {
        print("response shape must be either 'linear' or 'gaussian'.")
      }
      # calculate function
      curr_function <- linear_function(models$function_slope[i], models$function_intercept[i], curr_abundance)
      
      # generate output dataframe
      if(response_shape == 'linear') {
        data.frame(E = E,
                   species_ID = models$species_ID[i],
                   abundance_slope = models$abundance_slope[i],
                   abundance_intercept = models$abundance_intercept[i],
                   abundance = curr_abundance,
                   funct = curr_function)
      } else {
        data.frame(E = E,
                   species_ID = models$species_ID[i],
                   a = models$a[i],
                   b = models$b[i],
                   c = models$c[i],
                   abundance = curr_abundance,
                   funct = curr_function)
      }
    }))
  }))
  return(model_results)
}

# function to run n simulations
run_n_sims <- function(n_sims, n_species, environment_vals, params_ranges, response_shape, p_contribute, perfectly_crossing = FALSE){
  do.call(rbind, mclapply(1:n_sims, function(n){
    curr_results <- run_one_sim(n_species, environment_vals, params_ranges, response_shape, p_contribute, perfectly_crossing) %>%
      mutate(sim_number = n)
  }, mc.cores = num_cores))
}



# function to calculate weighted response diversity
weighted_rd <- function(df, response_shape){
  if(response_shape == 'linear'){
    rd <- df %>%
      group_by(sim_number, species_ID) %>%
      summarize(mean_abundance = mean(abundance), slope = median(abundance_slope)) %>%
      group_by(sim_number) %>%
      summarize(w_response_diversity = Hmisc::wtd.var(slope, mean_abundance)/mean(abs(slope)))
  } else if(response_shape == 'gaussian'){
    rd <- df %>%
      group_by(sim_number, species_ID) %>%
      summarize(mean_abundance = mean(abundance), a = median(a), b = median(b), c = median(c)) %>%
      group_by(sim_number) %>%
      summarize(w_response_diversity = Hmisc::wtd.var(b, mean_abundance)/mean(b))
  } else {
    print("response shape must be either 'linear' or 'gaussian'.")
    return(NULL)
  }
}

# function to calculate unweighted response diversity
unweighted_rd <- function(df, response_shape){
  if(response_shape == 'linear'){
    rd <- df %>%
      group_by(sim_number, species_ID) %>%
      summarize(slope = median(abundance_slope)) %>%
      group_by(sim_number) %>%
      summarize(uw_response_diversity = var(slope))
  } else if(response_shape == 'gaussian'){
    rd <- df %>%
      group_by(sim_number, species_ID) %>%
      summarize(a = median(a), b = median(b), c = median(c)) %>%
      group_by(sim_number) %>%
      summarize(uw_response_diversity = var(b))
  } else {
    print("response shape must be either 'linear' or 'gaussian'.")
    return(NULL)
  }
}

# function to calculate resilience and total function
resilience <- function(df){
  res <- df %>%
    group_by(sim_number, E) %>%
    summarize(total_function = sum(funct)) %>%
    group_by(sim_number) %>%
    summarize(resilience = 1/(var(total_function)/mean(total_function)),
              total_function = mean(total_function))
}

calc_stats <- function(df, response_shape){
  if(response_shape == 'linear'){
    # weighted RD
    wrd <- weighted_rd(df, 'linear')
    # unweighted RD
    urd <- unweighted_rd(df, 'linear')
    # resilience and total function
    res <- resilience(df)
    # join all together
    stats <- wrd %>%
      left_join(urd) %>%
      left_join(res)
  } else if(response_shape == 'gaussian'){
    # weighted RD
    wrd <- weighted_rd(df, 'gaussian')
    # unweighted RD
    urd <- unweighted_rd(df, 'gaussian')
    # resilience and total function
    res <- resilience(df)
    # join all together
    stats <- wrd %>%
      left_join(urd) %>%
      left_join(res)
  } else {
    print("response shape must be either 'linear' or 'gaussian'.")
    return(NULL)
  }
}

# functions to make plots
abundance_environment_plot <- function(df){
  ggplot(data = filter(df, sim_number == 1), aes(x = E, y = abundance, col = as.factor(species_ID))) +
    geom_line(size = 3) +
    #geom_smooth(method = "lm") +
    labs(x = "E", y = "Abundance", color = "Species") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25))
}

function_abundance_plot <- function(df){
  ggplot(data = filter(df, sim_number == 1), aes(x = abundance, y = funct, col = as.factor(species_ID))) +
    geom_line(size = 3) +
    #geom_smooth(method = "lm") +
    labs(x = "Abundance", y = "Function Provisioning", color = "Species") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25))
}

total_function_environment_plot <- function(df){
  df %>%
    group_by(sim_number, E) %>%
    summarize(total_function = sum(funct)) %>%
    ungroup() %>%
    filter(sim_number %in% 1:10) %>%
    ggplot(aes(x = E, y = total_function, col = as.factor(sim_number))) +
    geom_line(linewidth = 2) +
    labs(x = "E", y = "Total Function", col = "Ecosystem ID") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25))
}

resilience_w_rd_plot <- function(df){
  ggplot(data = df, aes(x = w_response_diversity, y = log(resilience))) +
    geom_point(size=3) +
    geom_smooth(method = "lm") +
    labs(x = "Weighted Response Diversity", y = "log(Resilience)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25))
}

resilience_u_rd_plot <- function(df){
  ggplot(data = df, aes(x = uw_response_diversity, y = log(resilience))) +
    geom_point(size=3) +
    geom_smooth(method = "lm") +
    labs(x = "Unweighted Response Diversity", y = "log(Resilience)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25))
}

function_w_rd_plot <- function(df){
  ggplot(data = df, aes(x = w_response_diversity, y = total_function)) +
    geom_point(size=3, alpha=0.3) +
    geom_smooth(method = "lm") +
    labs(x = "Weighted Response Diversity", y = "Total Function") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25))
}

function_u_rd_plot <- function(df){
  ggplot(data = df, aes(x = uw_response_diversity, y = total_function)) +
    geom_point(size=3, alpha=0.3) +
    geom_smooth(method = "lm") +
    labs(x = "Unweighted Response Diversity", y = "Total Function") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25))
}

resilience_function_w_rd_plot <- function(df){
  ggplot(data = df, aes(x = w_response_diversity, y = log(resilience), col = total_function)) +
    geom_point(size=3, alpha=0.8) +
    geom_smooth(method = "lm") +
    scale_color_viridis_c() +
    labs(x = "Response Diversity", y = "log(Stability)", col = "Ecosystem Function") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 15))
}

resilience_function_u_rd_plot <- function(df){
  ggplot(data = df, aes(x = uw_response_diversity, y = log(resilience), col = total_function)) +
    geom_point(size=3, alpha=0.8) +
    geom_smooth(method = "lm") +
    scale_color_viridis_c() +
    labs(x = "Unweighted Response Diversity", y = "log(Resilience)", col = "Total Function") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25))
}

p_contribute_plot <- function(df){
  ggplot(data = df, aes(x = proportion_contribute, y = estimate, col = p_val)) +
    geom_point(size = 5) +
    geom_errorbar(aes(ymin = conf_lower, ymax = conf_upper), linewidth = 2) +
    scale_color_viridis_c() +
    geom_smooth(method = "lm") +
    labs(x = "Proportion Species Contributing to Function", y = "Estimate", title = "Sensitivity Analysis") +
    theme_bw() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25),
          panel.grid = element_blank())
}







