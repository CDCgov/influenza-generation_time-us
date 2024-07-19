## mcmc_posterior_indep.m ######################################################
mcmc_posterior_indep <- function(scn) {
  # Import posterior parameter distributions for the independent transmission
  # and symptoms model obtained using data augmentation MCMC, and calculate
  # the posterior distribution of the proportion of presymptomatic
  # transmissions
  
  # Load assumed parameters
  load(paste0("../R/Results-temp/", scn, "/RData/assumed_parameters.RData"))
  
  # Load output of MCMC fitting procedure
  load(paste0("../R/Results-temp/", scn, "/RData/param_fit_indep.RData"))
  theta_mat <- result$theta_mat
  ll_vec    <- result$ll_vec
  output    <- result$output
  
  # Calculate posterior distributions of individual model parameters
  mean_post <- theta_mat[, 1]
  sd_post   <- theta_mat[, 2]
  beta_post <- theta_mat[, 3]
  
  # Point estimates of model parameters
  mean_best <- mean(mean_post)
  sd_best   <- mean(sd_post)
  beta_best <- mean(beta_post)
  
  # Posterior and point estimates for the proportion of presymptomatic transmissions
  F_inc      <- function(t) plnorm(t, meanlog = inc_mu, sdlog = inc_sigma)
  F_inc      <- function(t) pgamma(t, shape = inc_shape, scale = inc_scale) # replaced by gamma distributed incubation period (!)
  
  logn_mu    <- function(m, s) log(m^2 / sqrt(s^2 + m^2))
  logn_sigma <- function(m, s) sqrt(log(1 + s^2 / m^2))
  
  mu_post    <- logn_mu(mean_post, sd_post)
  sigma_post <- logn_sigma(mean_post, sd_post)
  
  mu_best    <- logn_mu(mean_best, sd_best)
  sigma_best <- logn_sigma(mean_best, sd_best)
  
  source("../R/Functions/Indep/get_presymp_trans_probs_indep_logn.R")
  prob_presymp_post <- get_presymp_trans_probs_indep_logn(mu_post, sigma_post, F_inc)
  prob_presymp_best <- get_presymp_trans_probs_indep_logn(mu_best, sigma_best, F_inc)
  
  # Save results
  mcmc_posterior_indep <- list(
    mean_post             = mean_post,
    mean_best             = mean_best,
    sd_post               = sd_post,
    sd_best               = sd_best,
    beta_post             = beta_post,
    beta_best             = beta_best,
    prob_presymp_post     = prob_presymp_post,
    prob_presymp_best     = prob_presymp_best,
    mu_post               = mu_post,
    sigma_post            = sigma_post,
    ll_vec                = ll_vec,
    empirical_summary_mat = output$empirical_summary_mat
  )
  save(mcmc_posterior_indep, file = paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_indep.RData"))
}
