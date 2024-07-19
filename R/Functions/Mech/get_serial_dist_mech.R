## get_serial_dist_mech.m ######################################################
get_serial_dist_mech <- function(params, int_max) {
  # Calculate the generation time distribution for our mechanistic
  # approach with parameters given by params.
  
  # This function requires the Chebfun package to execute (freely
  # available at https://www.chebfun.org/download/).
  
  gamma <- params[1]
  k_inc <- params[3]
  
  # Call the incubation distribution
  f_inc <- function(t_inc) dgamma(t_inc, shape = k_inc, scale = 1 / (k_inc * gamma))
  
  # Call the TOST distribution
  source("../R/Functions/Mech/f_tost_form_mech.R")
  f_tost <- function(t) f_tost_form_mech(params, t)
  
  # Define the serial interval distribution by convolution (Eq. 2)
  f_serial_mech <- function(x) {
    cubature::adaptIntegrate(function(t) f_tost(x - t) * f_inc(t), 0, Inf)$integral
    # cubature::adaptIntegrate(function(t) f_tost(t) * f_inc(x - t), 0, Inf)$integral
    # pracma::integral(function(t) f_tost(x - t) * f_inc(t), 0, int_max)
    # pracma::integral(function(t) f_tost(t) * f_inc(x - t), 0, int_max)
  }
  
  return(f_serial_mech)
}
