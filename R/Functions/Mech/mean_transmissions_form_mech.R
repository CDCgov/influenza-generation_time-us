## mean_transmissions_form_mech.m ##############################################
mean_transmissions_form_mech <- function(t_inc, household_size, asymp, params) {
  # Expected infectiousness of a host, conditional on incubation period
  # t_inc, integrated over the entire course of infection, under the
  # mechanistic approach with parameters given by params.
  
  gamma <- params[1]
  mu    <- params[2]
  k_inc <- params[3]
  k_E   <- params[4]
  k_I   <- params[5]
  alpha <- params[6]
  beta0 <- params[7]
  rho   <- params[8]
  x_A   <- params[9]
  k_P   <- k_inc - k_E
  
  beta_indiv        <- beta0 / (household_size^rho)
  beta_indiv[asymp] <- x_A * beta_indiv[asymp]
  
  mean_transmissions <- gamma * beta_indiv * (alpha * k_P * mu * t_inc + k_inc) / (alpha * k_P * mu + k_inc * gamma)
  return(mean_transmissions)
}
