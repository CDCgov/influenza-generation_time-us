## assumed_parameters.m ########################################################
assumed_parameters <- function(i, param_seq, scn) {
        # Assumed values of model parameters (other than those parameter values
        # that were fitted to the data)
        
        # Assumed parameters
        no_steps <- param_seq$no_steps[i]
        inc_mean <- param_seq$inc_mean[i]
        inc_sd   <- param_seq$inc_sd[i]
        x_A      <- param_seq$x_A[i]
        
        # # Estimate log-normal incubation period from Table 3 of Lessler et al
        # source("../R/Parameters/assumed_parameters_inc_flu.R")
        # inc_param <- assumed_parameters_inc_flu()
        # inc_shape <- inc_param[1]
        # inc_scale <- inc_param[2]
        # inc_mu    <- inc_param[3]
        # inc_sigma <- inc_param[4]
        # inc_gamma_mean  <- inc_param[5]
        # inc_gamma_sd    <- inc_param[6]
        # inc_lnorm_mean  <- inc_param[7]
        # inc_lnorm_sd    <- inc_param[8]
        
        # Log-normal incubation period
        inc_mu    <- log(inc_mean^2 / sqrt(inc_mean^2 + inc_sd^2))
        inc_sigma <- sqrt(log(1 + inc_sd^2 / inc_mean^2))
        # inc_mean   <- exp(inc_mu + 0.5 * inc_sigma^2)
        # inc_var    <- (exp(inc_sigma^2) - 1) * exp(2 * inc_mu + inc_sigma^2)
        f_inc_logn <- function(t_inc) dlnorm(t_inc, meanlog = inc_mu, sdlog = inc_sigma)
        
        # Gamma incubation period with the same mean and standard deviation
        inc_var  <- inc_sd^2
        inc_shape <- inc_mean^2 / inc_var
        inc_scale <- inc_var / inc_mean
        f_inc_gam <- function(t_inc) dgamma(t_inc, shape = inc_shape, scale = inc_scale)
        
        # Assumed parameter values common to both models
        # x_A <- 0.35  # relative infectiousness of asymptomatic infected hosts
        rho <- 1  # transmission scales with (household size)^(-rho)
        
        # Assumed parameter values in the mechanistic model
        k_inc <- inc_shape
        gamma <- 1 / (k_inc * inc_scale)
        k_I   <- 1
        
        # vector of known parameters
        params_known <- c(k_inc, gamma, k_I, rho, x_A)
        
        # Save to .RData file
        save(scn,
             no_steps,
             inc_mean,
             inc_sd,
             inc_var,
             inc_shape,
             inc_scale,
             inc_mu,
             inc_sigma,
             f_inc_logn,
             f_inc_gam,
             params_known,
             k_inc,
             gamma,
             k_I,
             rho,
             x_A,
             file = paste0("../R/Results-temp/", scn, "/RData/assumed_parameters.RData"))
}
