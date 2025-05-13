library(dplyr)
library(purrr)
library(minpack.lm)
library(furrr)

# Define the Gaussian function
gaussian <- function(x, a, b, c) {
  a * exp(-((x - b)^2) / (2 * c^2))
}

# define function to fit response curve for a species
fit_response_curve <- function(df, sp){
  #print(x)
  curr_species <- df %>%
    filter(species ==  sp) %>%
    ungroup()
  
  # check if there's more than 5 observations
  if(nrow(curr_species) >= 10){
    # drop outliers
    # Calculate Q1 and Q3
    Q1 <- quantile(curr_species$biomass_g_m2, 0.25)
    Q3 <- quantile(curr_species$biomass_g_m2, 0.75)
    
    # Compute the IQR
    IQR <- Q3 - Q1
    
    # Define the bounds
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    
    # Filter the vector to remove outliers
    curr_species_filtered <- curr_species %>%
      filter(biomass_g_m2 >= lower_bound & biomass_g_m2 <= upper_bound)
    
    # define bounds for parameters
    med_n <- median(curr_species_filtered$mean_max_temp)
    sd_n <- sd(curr_species_filtered$mean_max_temp)
    
    # Initial parameter estimates
    start_params <- list(a = Q3, b = med_n, c = sd_n)
    
    fit_result <- tryCatch(
      {
        fit <- minpack.lm::nlsLM(
          biomass_g_m2 ~ gaussian(mean_max_temp, a, b, c),
          data = curr_species_filtered,
          start = start_params,
          lower = c(-Inf, med_n-sd_n, sd_n*0.5),
          upper = c(Inf, med_n+sd_n, sd_n*1.5),
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
      print(paste0("model did not converge for current species: ", sp))
      print(fit_result$error)
      fitted_curve <- data.frame(
        species = NULL,
        mean_max_temp = NULL,
        predicted = NULL,
        a = NULL,
        b = NULL,
        c = NULL,
        biomass_g_m2 = NULL
      )
    } else {
      # Proceed with the fit object
      fit <- fit_result$fit
      # print(paste0("species: ", sp))
      # print(fit)
      # Generate fitted values for the original data points
      curr_species_filtered$predicted <- predict(fit, newdata = curr_species_filtered)
      
      # Generate fitted values for a smooth curve
      E_seq <- seq(min(curr_species$mean_max_temp), max(curr_species$mean_max_temp), length.out = 100)
      fitted_curve <- data.frame(
        species = sp,
        mean_max_temp = E_seq,
        predicted = predict(fit, newdata = data.frame(mean_max_temp = E_seq)),
        a = coef(fit)[1],
        b = coef(fit)[2],
        c = coef(fit)[3],
        biomass_g_m2 = mean(curr_species$biomass_g_m2)
      )
    }
  } else {
    print(paste0("not enough observations for current species: ", sp))
    fitted_curve <- data.frame(
      species = NULL,
      mean_max_temp = NULL,
      predicted = NULL,
      a = NULL,
      b = NULL,
      c = NULL,
      biomass_g_m2 = NULL
    )
  }
}


# functions to fit response curves with DOY as covariate
## gaussian with doy effect
gaussian_doy <- function(temp, doy_c, a, b, c, d) {
  a * exp(-((temp - b)^2) / (2 * c^2)) + d * doy_c
}

## iqr outlier filter
iqr_filter <- function(x, k = 1.5) {
  q <- quantile(x, c(.25, .75), na.rm = TRUE)
  iqr <- diff(q)
  x >= (q[1] - k * iqr) & x <= (q[2] + k * iqr)
}

fit_response_curve_doy <- function(dat_sp,
                                   species,
                                   min_n = 10,
                                   center_doy = TRUE) {
  if (nrow(dat_sp) < min_n) return(NULL)
  
  # center day‑of‑year so d represents within‑year timing
  dat_sp <- dat_sp %>%
    mutate(doy_c = if (center_doy)
      day_of_year - mean(day_of_year, na.rm = TRUE)
      else
        day_of_year)
  
  # simple IQR filter (optional)
  q  <- quantile(dat_sp$biomass_g_m2, c(.25, .75), na.rm = TRUE)
  iqr <- diff(q)
  dat_sp <- dat_sp %>%
    filter(biomass_g_m2 >= q[1] - 1.5 * iqr,
           biomass_g_m2 <= q[2] + 1.5 * iqr)
  
  # starting values & bounds
  start_vals <- list(
    a = max(dat_sp$biomass_g_m2, na.rm = TRUE),
    b = weighted.mean(dat_sp$mean_max_temp, dat_sp$biomass_g_m2),
    c = sd(dat_sp$mean_max_temp, na.rm = TRUE),
    # consider c=span/6
    d = 0
  )
  
  span <- diff(range(dat_sp$mean_max_temp))
  
  lower <- c(a = 0,
             b = min(dat_sp$mean_max_temp),
             c = 0.1,
             d = -Inf)
  upper <- c(a = Inf,
             b = max(dat_sp$mean_max_temp),
             c = span/2,
             d = Inf)
  
  fit <- tryCatch(
    nlsLM(biomass_g_m2 ~ gaussian_doy(mean_max_temp, doy_c, a, b, c, d),
          data = dat_sp,
          start = start_vals,
          lower = lower,
          upper = upper,
          control = nls.control(maxiter = 1000, minFactor = 1/1024)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  ## smooth prediction (DOY fixed at its mean, i.e. doy_c = 0)
  temp_seq <- seq(min(dat_sp$mean_max_temp),
                  max(dat_sp$mean_max_temp),
                  length.out = 100)
  newdata  <- tibble(mean_max_temp = temp_seq,
                     doy_c = 0)
  tibble(
    species = species,
    mean_max_temp = temp_seq,
    predicted = predict(fit, newdata = newdata),
    a = coef(fit)["a"],
    b = coef(fit)["b"],
    c = coef(fit)["c"],
    d = coef(fit)["d"],
    biomass_g_m2 = mean(dat_sp$biomass_g_m2)
  )
}

