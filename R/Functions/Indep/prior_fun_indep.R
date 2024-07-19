## prior_fun_indep.m ###########################################################
prior_fun_indep <- function(theta) {
  # Calculate the prior density for the independent transmission and
  # symptoms model with estimated parameters theta
  
  p_mean <- dlnorm(theta[1], meanlog = 1.6, sdlog = 0.35) # median 5, 95% CI [2.5,10]
  p_sd   <- dlnorm(theta[2], meanlog = 0.7, sdlog = 0.65) # median 2, 95% CI [0.6,7]
  p_beta <- dlnorm(theta[3], meanlog = 0.7, sdlog = 0.8)  # median 2, 95% CI [0.4,10]
  
  p <- p_mean * p_sd * p_beta
  return(p)
}
