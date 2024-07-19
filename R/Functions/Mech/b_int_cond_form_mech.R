## b_int_cond_form_mech.m ######################################################
b_int_cond_form_mech <- function(x, t_inc, household_size, asymp, params) {
  # Expected infectiousness of a host, conditional on incubation period
  # t_inc, integrated from times (-infinity) to x since symptom onset,
  # under the mechanistic approach with parameters given by params.
  
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
  t_inc_p      <- t_inc[ind_p]
  beta_indiv_p <- beta_indiv[ind_p]
  
  f_m1 <- alpha * C * beta_indiv_m * x_m * (1 - pbeta(-x_m / t_inc_m, k_P, k_E))
  f_m2 <- (alpha * C * beta_indiv_m * k_P * t_inc_m / k_inc) * (1 - pbeta(-x_m / t_inc_m, k_P + 1, k_E))
  f_m  <- f_m1 + f_m2
  
  f_p1 <- C * beta_indiv_p * alpha * k_P * t_inc_p / k_inc
  f_p2 <- C * beta_indiv_p * x_p * (1 - pgamma(x_p, shape = k_I, scale = 1 / (k_I * mu)))
  f_p3 <- (C * beta_indiv_p / mu) * pgamma(x_p, shape = k_I + 1, scale = 1 / (k_I * mu))
  f_p  <- f_p1 + f_p2 + f_p3
  
  B_cond        <- rep(0, length(x))
  B_cond[ind_m] <- f_m
  B_cond[ind_p] <- f_p
  
  return(B_cond)
}
