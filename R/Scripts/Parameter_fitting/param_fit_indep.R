## param_fit_indep.m ###########################################################
param_fit_indep <- function(i, param_seq, scn) {
  # Fit the parameters of the independent transmission and symptoms model
  # using data augmentation MCMC
  
  # The vector of unknown model parameters, theta  = [m_gen,s_gen,beta_0], is
  # estimated in the model fitting procedure. Here, m_gen and s_gen are the
  # mean and standard deviation of the generation time distribution,
  # respectively, and beta_0 is the overall infectiousness parameter
  
  # Assumed parameters
  no_steps <- param_seq$no_steps[i]
  inc_mean <- param_seq$inc_mean[i]
  inc_sd   <- param_seq$inc_sd[i]
  
  # Load data
  load(paste0("../R/Results-temp/", scn, "/RData/data.RData"))
  
  # Load incubation period distribution function and assumed parameters
  load(paste0("../R/Results-temp/", scn, "/RData/assumed_parameters.RData"))
  f_inc <- f_inc_logn
  f_inc <- f_inc_gam  # replaced by gamma distributed incubation period (!)
  
  # Number of steps in chain, keeping only some steps after burn-in and thinning
  # no_steps <- 1000 # 10000000 in the uk study
  steps_keep <- seq(((no_steps / 5) + 1), no_steps, by = 100)
  
  # Function handles giving the density and cumulative distribution functions
  # of the generation time for parameters theta
  logn_mu    <- function(m, s) log(m^2 / sqrt(s^2 + m^2))
  logn_sigma <- function(m, s) sqrt(log(1 + s^2 / m^2))
  
  f_gen_form <- function(t_gen, theta) dlnorm(t_gen, meanlog = logn_mu(theta[1], theta[2]), sdlog = logn_sigma(theta[1], theta[2]))
  F_gen_form <- function(t_gen, theta) plnorm(t_gen, meanlog = logn_mu(theta[1], theta[2]), sdlog = logn_sigma(theta[1], theta[2]))
  
  # Function handle giving the log-likelihood for parameters theta and augmented data data_struct_augmented
  source("../R/Functions/Indep/log_likelihood_household_indep.R")
  ll_household_form <- function(theta, data_struct_augmented) {
    log_likelihood_household_indep(f_inc, theta[3], rho, x_A,
                                   function(t_gen) f_gen_form(t_gen, theta),
                                   function(t_gen) F_gen_form(t_gen, theta),
                                   data_struct_augmented)
  }
  
  # Function handle giving the estimated mean and standard deviation of
  # realised household generation times, and the presymptomatic proportion of
  # transmissions, for parameters theta and augmented data data_struct_augmented
  source("../R/Functions/Indep/empirical_summary_indep.R")
  empirical_summary_form <- function(theta, data_struct_augmented) empirical_summary_indep(theta[3], rho, x_A, function(t_gen) f_gen_form(t_gen, theta), data_struct_augmented)
  
  # Standard deviations of proposal distributions for individual model parameters
  sd_prop_mean <- 0.3
  sd_prop_sd   <- 2.5  * sd_prop_mean
  sd_prop_beta <- 0.25 * sd_prop_mean
  
  # Correlation between proposal distributions
  corr_prop_mat <- diag(3)
  
  # Covariance matrix of multivariate normal proposal distribution for the
  # vector of fitted model parameters, theta
  sd_prop_diag       <- diag(c(sd_prop_mean, sd_prop_sd, sd_prop_beta))
  theta_prop_cov_mat <- sd_prop_diag %*% corr_prop_mat %*% sd_prop_diag
  
  # Initial values of model parameters
  theta_init <- c(5, 5, 2)
  
  # Standard deviation of proposal distributions for infection times of
  # symptomatic and asymptomatic infected hosts
  t_i_prop_sd_symp  <- 8
  t_i_prop_sd_asymp <- 13
  
  # Function handle giving the prior density for parameters theta
  source("../R/Functions/Indep/prior_fun_indep.R")
  prior_fun <- function(theta) prior_fun_indep(theta)
  
  # Function handles used to update either the model parameters or augmented
  # data during the fitting procedure
  source("../R/Functions/MCMC/update_theta_fun.R")
  source("../R/Functions/MCMC/update_infection_fun_indep.R")
  source("../R/Functions/MCMC/update_onset_fun.R")
  source("../R/Functions/MCMC/update_asymp_fun_indep.R")
  update_theta     <- function(theta, data_struct_augmented, ll_household) update_theta_fun(theta, data_struct_augmented, ll_household, ll_household_form, theta_prop_cov_mat, prior_fun)
  update_infection <- function(theta, data_struct_augmented, ll_household) update_infection_fun_indep(theta, data_struct_augmented, ll_household, ll_household_form, t_i_prop_sd_symp)
  update_onset     <- function(theta, data_struct_augmented, ll_household) update_onset_fun(theta, data_struct_augmented, ll_household, ll_household_form)
  update_asymp     <- function(theta, data_struct_augmented, ll_household) update_asymp_fun_indep(theta, data_struct_augmented, ll_household, ll_household_form, t_i_prop_sd_asymp)
  
  # Initialise augmented data and calculate the initial likelihood
  source("../R/Functions/MCMC/initialise_augmented_data_indep.R")
  data_struct_augmented_init <- initialise_augmented_data_indep(data_struct_observed)
  ll_household_init          <- ll_household_form(theta_init, data_struct_augmented_init)
  
  # Run main parameter fitting procedure
  tic <- Sys.time()
  plotting <- FALSE # set to TRUE to plot output of fitting procedure as it is carried out (updating the plot once for each 1% of steps completed)
  source("../R/Functions/MCMC/fit_params.R")
  result <- fit_params(scn, no_steps, steps_keep, update_theta, update_infection, update_onset, update_asymp, empirical_summary_form, theta_init, data_struct_augmented_init, ll_household_init, plotting)
  print(toc <- Sys.time() - tic)
  
  # Save results
  save(result, file = paste0("../R/Results-temp/", scn, "/RData/param_fit_indep.RData"))
}
