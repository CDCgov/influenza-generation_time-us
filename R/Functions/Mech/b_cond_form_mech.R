## b_cond_form_mech.m ##########################################################
b_cond_form_mech <- function(x, t_inc, household_size, asymp, params) {
  # Expected infectiousness of a host at time x since symptom onset,
  # conditional on incubation period t_inc, under the mechanistic
  # approach with parameters given by params.
  
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
  
  C                 <- k_inc * gamma * mu / (alpha * k_P * mu + k_inc * gamma)
  beta_indiv        <- beta0 / (household_size^rho)
  beta_indiv[asymp] <- x_A * beta_indiv[asymp]
  
  ind_m <- (x < 0)
  ind_p <- (!ind_m)
  
  x_m          <- x[ind_m]
  t_inc_m      <- t_inc[ind_m]
  beta_indiv_m <- beta_indiv[ind_m]
  
  x_p          <- x[ind_p]
  beta_indiv_p <- beta_indiv[ind_p]
  
  f_m <- alpha * C * beta_indiv_m * (1 - pbeta(-x_m / t_inc_m, k_P, k_E))
  f_p <- C * beta_indiv_p * (1 - pgamma(x_p, shape = k_I, scale = 1 / (k_I * mu)))
  
  b_cond        <- rep(0, length(x))
  b_cond[ind_m] <- f_m
  b_cond[ind_p] <- f_p
  
  return(b_cond)
}
