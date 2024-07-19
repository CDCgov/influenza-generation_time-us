## param_fit_mech.m ############################################################
param_fit_mech <- function(i, param_seq, scn) {
  # Fit the parameters of the mechanistic model using data augmentation MCMC
  
  # The vector of unknown model parameters, theta  =
  # [p_E,mu_inv,alpha,beta_0], is estimated in the model fitting procedure.
  # Here, p_E=k_E/k_inc is the ratio between the shape parameters of the
  # Gamma distributed latent (E) and incubation (combined E and P) periods,
  # mu_inv=1/mu is the mean symptomatic infectious (I) period, alpha=alpha_P
  # is the ratio of the transmission rates during the presymptomatic
  # infectious (P) and symptomatic infectious (I) periods, and beta_0 is the
  # overall infectiousness parameter
  
  # Assumed parameters
  no_steps <- param_seq$no_steps[i]
  inc_mean <- param_seq$inc_mean[i]
  inc_sd   <- param_seq$inc_sd[i]
  
  # Load data
  load(paste0("../R/Results-temp/", scn, "/RData/data.RData"))
  
  # Load incubation period distribution function and assumed parameters
  load(paste0("../R/Results-temp/", scn, "/RData/assumed_parameters.RData"))
  f_inc <- f_inc_gam
  
  # Number of steps in chain, keeping only some steps after burn-in and thinning
  # no_steps <- 1000 # 10000000 in the uk study
  steps_keep <- seq(((no_steps / 5) + 1), no_steps, by = 100)
  
  # Function handles giving (i) the expected infectiousness at each time
  # since symptom onset x, conditional on incubation period t_inc; (ii) the
  # integral of this infectiousness profile between times since onset
  # (-infinity) and x; and (iii) the integral over all times since infection;
  # for parameters theta
  source("../R/Functions/Mech/get_params_mech.R")
  source("../R/Functions/Mech/b_cond_form_mech.R")
  source("../R/Functions/Mech/b_int_cond_form_mech.R")
  source("../R/Functions/Mech/mean_transmissions_form_mech.R")
  b_cond_form             <- function(x, t_inc, household_size, asymp, theta) b_cond_form_mech(x, t_inc, household_size, asymp, get_params_mech(theta, params_known))
  B_cond_form             <- function(x, t_inc, household_size, asymp, theta) b_int_cond_form_mech(x, t_inc, household_size, asymp, get_params_mech(theta, params_known))
  mean_transmissions_form <- function(t_inc, household_size, asymp, theta) mean_transmissions_form_mech(t_inc, household_size, asymp, get_params_mech(theta, params_known))
  
  # Function handle giving the log-likelihood for parameters theta and
  # augmented data data_struct_augmented
  source("../R/Functions/Mech/log_likelihood_household_mech.R")
  ll_household_form <- function(theta, data_struct_augmented) {
    log_likelihood_household_mech(f_inc,
                                  function(x, t_inc, household_size, asymp) b_cond_form(x, t_inc, household_size, asymp, theta),
                                  function(x, t_inc, household_size, asymp) B_cond_form(x, t_inc, household_size, asymp, theta),
                                  function(t_inc, household_size, asymp) mean_transmissions_form(t_inc, household_size, asymp, theta),
                                  data_struct_augmented)
  }
  
  # Function handle giving the estimated mean and standard deviation of
  # realised household generation times, and the presymptomatic proportion of
  # transmissions, for parameters theta and augmented data
  # data_struct_augmented
  source("../R/Functions/Mech/empirical_summary_mech.R")
  empirical_summary_form <- function(theta, data_struct_augmented) empirical_summary_mech(function(x, t_inc, household_size, asymp) b_cond_form(x, t_inc, household_size, asymp, theta), data_struct_augmented)
  
  # Standard deviations of proposal distributions for individual model parameters
  sd_prop_p_E    <- 0.075
  sd_prop_mu_inv <- 10  * sd_prop_p_E
  sd_prop_alpha  <- 10  * sd_prop_p_E
  sd_prop_beta   <- 1.4 * sd_prop_p_E
  
  # Correlation between proposal distributions
  corr_p_E_mu_inv   <- 0
  corr_p_E_alpha    <- 0.5
  corr_p_E_beta     <- 0
  corr_mu_inv_alpha <- 0.5
  corr_mu_inv_beta  <- 0
  corr_alpha_beta   <- 0
  
  corr_prop_mat <- matrix(c(1,               corr_p_E_mu_inv,   corr_p_E_alpha,    corr_p_E_beta,
                            corr_p_E_mu_inv, 1,                 corr_mu_inv_alpha, corr_mu_inv_beta,
                            corr_p_E_alpha,  corr_mu_inv_alpha, 1,                 corr_alpha_beta,
                            corr_p_E_beta,   corr_mu_inv_beta,  corr_alpha_beta,   1),
                          nrow = 4, byrow = TRUE)
  
  # Covariance matrix of multivariate normal proposal distribution for the
  # vector of fitted model parameters, theta
  sd_prop_diag <- diag(c(sd_prop_p_E, sd_prop_mu_inv, sd_prop_alpha, sd_prop_beta))
  theta_prop_cov_mat <- sd_prop_diag %*% corr_prop_mat %*% sd_prop_diag
  
  # Initial values of model parameters
  theta_init <- c(0.5, 1 / 0.18, 3.5, 2)
  
  # Standard deviation of proposal distributions for infection times of
  # infected hosts, and for simultaneous updates of the times at which
  # asymptomatic hosts become infected and enter the I stage (holding the
  # difference between the two times constant; note the distinction between
  # the P and I stages has no epidemiological meaning for asymptomatic hosts,
  # but the observed data are augmented with the time of entering the I stage
  # when fitting the model for convenience)
  t_i_prop_sd     <- 9
  t_prop_sd_asymp <- 18
  
  # Function handle giving the prior density for parameters theta
  source("../R/Functions/Mech/prior_fun_mech.R")
  prior_fun <- function(theta) prior_fun_mech(theta)
  
  # Function handles used to update either the model parameters or augmented
  # data during the fitting procedure
  source("../R/Functions/MCMC/update_theta_fun.R")
  source("../R/Functions/MCMC/update_infection_fun_mech.R")
  source("../R/Functions/MCMC/update_onset_fun.R")
  source("../R/Functions/MCMC/update_asymp_fun_mech.R")
  update_theta     <- function(theta, data_struct_augmented, ll_household) update_theta_fun(theta, data_struct_augmented, ll_household, ll_household_form, theta_prop_cov_mat, prior_fun)
  update_infection <- function(theta, data_struct_augmented, ll_household) update_infection_fun_mech(theta, data_struct_augmented, ll_household, ll_household_form, t_i_prop_sd)
  update_onset     <- function(theta, data_struct_augmented, ll_household) update_onset_fun(theta, data_struct_augmented, ll_household, ll_household_form)
  update_asymp     <- function(theta, data_struct_augmented, ll_household) update_asymp_fun_mech(theta, data_struct_augmented, ll_household, ll_household_form, t_prop_sd_asymp)
  
  # Initialize augmented data and calculate the initial likelihood
  source("../R/Functions/MCMC/initialise_augmented_data_mech.R")
  data_struct_augmented_init <- initialise_augmented_data_mech(data_struct_observed)
  ll_household_init          <- ll_household_form(theta_init, data_struct_augmented_init)
  
  # Run main parameter fitting procedure
  tic <- Sys.time()
  plotting <- FALSE # set to TRUE to plot output of fitting procedure as it is carried out (updating the plot once for each 1% of steps completed)
  source("../R/Functions/MCMC/fit_params.R")
  result <- fit_params(scn, no_steps, steps_keep, update_theta, update_infection, update_onset, update_asymp, empirical_summary_form, theta_init, data_struct_augmented_init, ll_household_init, plotting)
  print(toc <- Sys.time() - tic)
  
  # Save results
  save(result, file = paste0("../R/Results-temp/", scn, "/RData/param_fit_mech.RData"))
}
