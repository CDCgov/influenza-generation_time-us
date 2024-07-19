## get_gen_dist_mech.m #########################################################
get_gen_dist_mech <- function(params, int_max) {
  # Calculate the generation time distribution for our mechanistic
  # approach with parameters given by params.
  
  # This function requires the Chebfun package to execute (freely
  # available at https://www.chebfun.org/download/).
  
  gamma <- params[1]
  mu    <- params[2]
  k_inc <- params[3]
  k_E   <- params[4]
  k_I   <- params[5]
  alpha <- params[6]
  k_P   <- k_inc - k_E
  
  # parameter (page 23)
  C <- k_inc * gamma * mu / (alpha * k_P * mu + k_inc * gamma)
  
  # pdf and cdf of duration of each state (E, P, I)
  f_E <- function(t) dgamma(t, shape = k_E, scale = 1 / (k_inc * gamma))
  f_P <- function(t) dgamma(t, shape = k_P, scale = 1 / (k_inc * gamma))
  F_P <- function(t) pgamma(t, shape = k_P, scale = 1 / (k_inc * gamma))
  F_I <- function(t) pgamma(t, shape = k_I, scale = 1 / (k_I * mu))
  
  # Define the distribution (page 24)
  f1 <- function(x) {
    cubature::adaptIntegrate(function(t) f_P(t) * F_I(x - t), 0, Inf)$integral
    # cubature::adaptIntegrate(function(t) f_P(x - t) * F_I(t), 0, Inf)$integral
    # pracma::integral(function(t) f_P(t) * F_I(x - t), 0, int_max)
    # pracma::integral(function(t) f_P(x - t) * F_I(t), 0, int_max)
  }
  
  # Define the distribution (page 24)
  f_star <- function(x) {
    (alpha * C * (1 - F_P(x)) + C * (F_P(x) - f1(x))) * (x > 0)
  }
  
  # Define the generation time distribution (page 24)
  f_gen_mech <- function(x) {
    cubature::adaptIntegrate(function(t) f_star(x - t) * f_E(t), 0, Inf)$integral
    # cubature::adaptIntegrate(function(t) f_star(t) * f_E(x - t), 0, Inf)$integral
    # pracma::quadv(function(t) f_star(x - t) * f_E(t), 0, int_max)$Q
    # pracma::quadv(function(t) f_star(t) * f_E(x - t), 0, int_max)$Q
  }
  
  return(f_gen_mech)
}
