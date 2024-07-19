## update_theta_fun.m ##########################################################
update_theta_fun <- function(theta_old, data_struct_augmented, ll_household_old, ll_household_form, theta_prop_cov_mat, prior_fun) {
  # Update the vector of fitted model parameters, theta, in the parameter
  # fitting procedure for either model
  
  # The log-likelihood using the proposed parameters
  no_households     <- length(ll_household_old)
  ll_household_prop <- rep(NaN, no_households)
  
  # Avoid NaN in log-likelihood due to theta[1] > 1 in the mechanistic model
  while (any(is.nan(ll_household_prop))) {
    # Propose new values of each entry of theta using a multivariate normal proposal distribution
    no_params_fitted <- length(theta_old)
    theta_prop       <- rep(-1, no_params_fitted)
    
    # Avoid negative theta
    while (any(theta_prop < 0)) {
      theta_prop       <- theta_old + MASS::mvrnorm(1, mu = rep(0, no_params_fitted), Sigma = theta_prop_cov_mat)
    }
    
    # Calculate the log-likelihood using the proposed parameters
    ll_household_prop <- ll_household_form(theta_prop, data_struct_augmented)
  }
  
  # Calculate the logarithm of the ratio between the proposed and previous likelihoods
  a1 <- exp(sum(ll_household_prop) - sum(ll_household_old))
  a1 <- ifelse(is.nan(a1), 0, a1)
  
  # Calculate the acceptance probability, a, for the new parameters, accounting for priors
  a <- a1 * prior_fun(theta_prop) / prior_fun(theta_old)
  
  # Accept the proposed parameters and onset times with probability a
  if (runif(1) < a) {
    # Accept the proposed parameters
    theta_new <- theta_prop
    
    # New likelihood
    ll_household_new <- ll_household_prop
    
    # Information about acceptance
    acceptance <- list(
      overall = 1
    )
    
  } else {
    # Reject the proposed parameters
    theta_new <- theta_old
    
    # New likelihood
    ll_household_new <- ll_household_old
    
    # Information about acceptance
    acceptance <- list(
      overall = 0
    )
  }
  
  return(list(theta        = theta_new,
              ll_household = ll_household_new,
              acceptance   = acceptance))
}
