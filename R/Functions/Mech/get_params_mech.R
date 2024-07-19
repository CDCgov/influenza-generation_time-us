## get_params_mech.m ###########################################################
get_params_mech <- function(theta, params_known) {
  # Recover the full vector of parameters for our mechanistic approach,
  # params = [gamma,mu,k_inc,k_E,k_I,alpha,rho,x_A], from the vector of
  # parameters updated in the MCMC model fitting procedure for the
  # variable infectiousness model, theta = [k_E/k_inc,1/mu,alpha,beta],
  # and the vector of known parameters, params_known =
  # [k_inc,gamma,k_I,rho,x_A]
  
  k_inc <- params_known[1]
  gamma <- params_known[2]
  k_I   <- params_known[3]
  rho   <- params_known[4]
  x_A   <- params_known[5]
  
  k_E   <- theta[1] * k_inc
  mu    <- 1 / theta[2]
  alpha <- theta[3]
  beta  <- theta[4]
  
  params <- c(gamma, mu, k_inc, k_E, k_I, alpha, beta, rho, x_A)
  return(params)
}
