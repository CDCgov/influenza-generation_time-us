## assumed_parameters_inc_flu.m ################################################
assumed_parameters_inc_flu <- function() {
  # library(epiparameter)
  # eparams <- epiparameter::epiparam()                             # Load epidemic parameter data
  # influenza_incubation <- epiparameter::as_epidist(eparams[56, ]) # Extract incubation period information for influenza
  # influenza_incubation$summary_stats$quantiles                    # Get the quantiles of the incubation period summary statistics
  
  # Given percentiles and corresponding values from Table 3 of Lessler et al
  percentiles <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  observed_values <- c(0.7, 1.1, 1.4, 1.9, 2.8)  # Influenza A
  # observed_values <- c(0.3, 0.4, 0.6, 0.7, 1.1)  # Influenza B
  
  # Function to estimate gamma distribution parameters (inc_shape and inc_scale)
  estimate_parameters_gamma <- function(parameters) {
    # Two parameters of the gamma distribution
    inc_shape <- parameters[1]
    inc_scale <- parameters[2]
    
    # Express percentiles in terms of inc_shape and inc_scale
    estimated_values <- qgamma(p = percentiles, shape = inc_shape, scale = inc_scale)
    
    # Calculate the sum of squared differences between estimated and observed percentiles
    sum_squared_diff <- sum((estimated_values - observed_values)^2)
    
    return(sum_squared_diff)
  }
  
  # Function to estimate log-normal distribution parameters (inc_mu and inc_sigma)
  estimate_parameters_lnorm <- function(parameters) {
    # Two parameters of the log-normal distribution
    inc_mu    <- parameters[1]
    inc_sigma <- parameters[2]
    
    # Express percentiles in terms of inc_mu and inc_sigma
    estimated_values <- qlnorm(p = percentiles, meanlog = inc_mu, sdlog = inc_sigma)
    
    # Calculate the sum of squared differences between estimated and observed percentiles
    sum_squared_diff <- sum((estimated_values - observed_values)^2)
    
    return(sum_squared_diff)
  }
  
  # Initial parameter values for optimization
  initial_parameters <- c(0, 1)  # You can start with any initial values
  
  # Perform parameter estimation using optimization
  estimated_parameters_gamma <- optim(par = initial_parameters, fn = estimate_parameters_gamma)
  estimated_parameters_lnorm <- optim(par = initial_parameters, fn = estimate_parameters_lnorm)
  
  # Extract the estimated parameters (gamma distribution)
  inc_shape <- estimated_parameters_gamma$par[1]
  inc_scale <- estimated_parameters_gamma$par[2]
  # Print the estimates
  cat("(gamma distribution) Estimated shape:", inc_shape, "\n")
  cat("(gamma distribution) Estimated scale:", inc_scale, "\n")
  
  # Extract the estimated parameters (log-normal distribution)
  inc_mu    <- estimated_parameters_lnorm$par[1]
  inc_sigma <- estimated_parameters_lnorm$par[2]
  # Print the estimates
  cat("(log-normal distribution) Estimated inc_mu:",    inc_mu, "\n")
  cat("(log-normal distribution) Estimated inc_sigma:", inc_sigma, "\n")
  
  # The mean and variance of incubation period for gamma distribution
  inc_gamma_mean <- inc_shape * inc_scale
  inc_gamma_sd   <- sqrt(inc_shape * inc_scale^2)
  # Print the estimates
  cat("(gamma distribution) Estimated inc_mean:", inc_gamma_mean, "\n")
  cat("(gamma distribution) Estimated inc_sd:",   inc_gamma_sd, "\n")
  
  # The mean and variance of incubation period for log-normal distribution
  inc_lnorm_mean  <- exp(inc_mu + 0.5 * inc_sigma^2)
  inc_lnorm_sd    <- sqrt((exp(inc_sigma^2) - 1) * exp(2 * inc_mu + inc_sigma^2))
  # Print the estimates
  cat("(log-normal distribution) Estimated inc_mean:",    inc_lnorm_mean, "\n")
  cat("(log-normal distribution) Estimated inc_sd:",      inc_lnorm_sd, "\n")
  
  return(c(inc_shape, inc_scale,
           inc_mu, inc_sigma,
           inc_gamma_mean, inc_gamma_sd,
           inc_lnorm_mean, inc_lnorm_sd))
}
