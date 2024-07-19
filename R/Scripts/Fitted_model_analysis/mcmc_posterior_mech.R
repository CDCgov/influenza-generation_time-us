## mcmc_posterior_mech.m #######################################################
mcmc_posterior_mech <- function(scn) {
  # Import posterior parameter distributions for the mechanistic model
  # obtained using data augmentation MCMC, and calculate the posterior
  # distributions of the mean and standard deviation of generation times and
  # the proportion of presymptomatic transmissions
  
  # Load assumed parameters
  load(paste0("../R/Results-temp/", scn, "/RData/assumed_parameters.RData"))
  
  # Load output of MCMC fitting procedure
  load(paste0("../R/Results-temp/", scn, "/RData/param_fit_mech.RData"))
  theta_mat <- result$theta_mat
  ll_vec    <- result$ll_vec
  output    <- result$output
  
  # Calculate posterior distributions of individual model parameters
  p_E_post    <- theta_mat[, 1]
  k_E_post    <- k_inc * p_E_post
  k_P_post    <- k_inc - k_E_post
  mu_inv_post <- theta_mat[, 2]
  alpha_post  <- theta_mat[, 3]
  beta_post   <- theta_mat[, 4]
  
  # Point estimates of model parameters
  p_E_best    <- mean(p_E_post)
  k_E_best    <- k_inc * p_E_best
  k_P_best    <- k_inc - k_E_best
  mu_inv_best <- mean(mu_inv_post)
  alpha_best  <- mean(alpha_post)
  beta_best   <- mean(beta_post)
  
  theta_best <- c(p_E_best, mu_inv_best, alpha_best, beta_best) # fitted parameters
  
  source("../R/Functions/Mech/get_params_mech.R")
  params_best <- get_params_mech(theta_best, params_known) # all parameters
  params_best <- matrix(params_best, ncol = length(params_best))
  
  # Posterior and point estimates for the proportion of presymptomatic transmissions
  prob_presymp_post <- (alpha_post * k_P_post / (k_inc * gamma)) / ((alpha_post * k_P_post / (k_inc * gamma)) + mu_inv_post)
  prob_presymp_best <- (alpha_best * k_P_best / (k_inc * gamma)) / ((alpha_best * k_P_best / (k_inc * gamma)) + mu_inv_best)
  
  # Posterior and point estimates for the mean and standard deviation of generation times
  no_steps_kept <- length(k_E_post)
  k_inc_post    <- rep(params_known[1], no_steps_kept)
  gamma_post    <- rep(params_known[2], no_steps_kept)
  k_I_post      <- rep(params_known[3], no_steps_kept)
  rho_post      <- rep(params_known[4], no_steps_kept)
  x_A_post      <- rep(params_known[5], no_steps_kept)
  params_post   <- cbind(gamma_post, 1 / mu_inv_post, k_inc_post, k_E_post, k_I_post, alpha_post, beta_post, rho_post, x_A_post)
  
  source("../R/Functions/Mech/get_gen_mean_sd_mech.R")
  gen_post  <- get_gen_mean_sd_mech(params_post)
  gen_best  <- get_gen_mean_sd_mech(params_best)
  
  mean_post <- gen_post$m_gen
  mean_best <- gen_best$m_gen
  
  sd_post   <- gen_post$s_gen
  sd_best   <- gen_best$s_gen
  
  # Save results
  mcmc_posterior_mech <- list(
    p_E_post              = p_E_post,
    mu_inv_post           = mu_inv_post,
    alpha_post            = alpha_post,
    beta_post             = beta_post,
    prob_presymp_post     = prob_presymp_post,
    params_best           = params_best,
    mean_post             = mean_post,
    mean_best             = mean_best,
    sd_post               = sd_post,
    sd_best               = sd_best,
    ll_vec                = ll_vec,
    empirical_summary_mat = output$empirical_summary_mat
  )
  save(mcmc_posterior_mech, file = paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_mech.RData"))
}
