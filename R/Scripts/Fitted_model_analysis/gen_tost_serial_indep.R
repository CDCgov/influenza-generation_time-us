## gen_tost_serial_indep.m #####################################################
gen_tost_serial_indep <- function(scn) {
        # Calculate distributions of the generation time, time from onset of
        # symptoms to transmission (TOST) and serial interval for the independent
        # transmission and symptoms model, using point estimate (posterior mean)
        # parameter values
        
        # This script requires the Chebfun package to run (freely available at
        # https://www.chebfun.org/download/)
        
        # Define constants for the distributions
        dt      <- 0.1
        ti      <- seq(-50, 50, by = dt)
        int_max <- max(ti)
        
        plotting <- FALSE
        
        # Load incubation period distribution
        load(paste0("../R/Results-temp/", scn, "/RData/assumed_parameters.RData"))
        f_inc <- f_inc_logn
        f_inc <- f_inc_gam  # replaced by gamma distributed incubation period (!)
        
        # Load point estimate parameter values
        load(paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_indep.RData"))
        mean_best <- mcmc_posterior_indep$mean_best
        sd_best   <- mcmc_posterior_indep$sd_best
        
        logn_mu    <- function(m, s) log(m^2 / sqrt(s^2 + m^2))
        logn_sigma <- function(m, s) sqrt(log(1 + s^2 / m^2))
        
        mu_best    <- logn_mu(mean_best, sd_best)
        sigma_best <- logn_sigma(mean_best, sd_best)
        
        ## Generation time #############################################################
        # Define the generation time distribution
        f_gen_indep <- function(t) dlnorm(t, meanlog = mu_best, sdlog = sigma_best)
        
        # Numerical approximation
        f_gen_indep_ti <- f_gen_indep(ti)
        
        # Checking
        if (sum(f_gen_indep_ti * dt) < 0.99) {
                warning("Difficulties with the integration of the generation time distribution, f_gen_indep (!)")
        } else {
                # Plot the generation time distribution
                if (plotting) {
                        plot(ti, f_gen_indep_ti,
                             xlim = c(0, 15),
                             main = paste0("Mean = ", sprintf("%.1f", sum(f_gen_indep_ti * ti) * dt),
                                           ", SD = ", sprintf("%.1f", sqrt(sum(f_gen_indep_ti * ti^2) * dt - (sum(f_gen_indep_ti * ti) * dt)^2))))
                }
        }
        
        ## TOST ########################################################################
        # Define the TOST distribution by convolution (Eq. 1)
        f_tost_indep <- function(x) {
                cubature::adaptIntegrate(function(t) f_gen_indep(x + t) * f_inc(t), 0, Inf)$integral
                # pracma::integral(function(t) f_gen_indep(x + t) * f_inc(t), 0, int_max)
        }
        
        # Numerical approximation
        tic <- Sys.time()
        f_tost_indep_ti <-
                tryCatch({
                        sapply(ti, f_tost_indep)
                }, error = function(err) {
                        warning("An error occurred: ", conditionMessage(err), "\n")
                        rep(0, length(ti))                               # Return a vector of zeros
                })
        print(toc <- Sys.time() - tic)
        
        # # Numerical approximation
        # tic <- Sys.time()
        # cl <- parallel::makeCluster(parallel::detectCores() - 1)         # Create a cluster with the detected number of CPU cores
        # parallel::clusterExport(cl, c(ls(), lsf.str()))                  # Export the all variables and functions to the cluster
        # f_tost_indep_ti <- parallel::parSapply(cl, ti, f_tost_indep)     # Use parSapply to apply the function in parallel
        # parallel::stopCluster(cl)                                        # Stop the cluster when done
        # print(toc <- Sys.time() - tic)
        
        # Checking
        if (sum(f_tost_indep_ti * dt) < 0.99) {
                warning("Difficulties with the integration of the TOST distribution, f_tost_indep (!)")
        } else {
                # Plot the TOST distribution
                if (plotting) {
                        plot(ti, f_tost_indep_ti,
                             xlim = c(-10, 10),
                             main = paste0("Mean = ", sprintf("%.1f", sum(f_tost_indep_ti * ti) * dt),
                                           ", SD = ", sprintf("%.1f", sqrt(sum(f_tost_indep_ti * ti^2) * dt - (sum(f_tost_indep_ti * ti) * dt)^2))))
                }
        }
        
        ## Serial interval #############################################################
        # Define the serial interval distribution by convolution (Eq. 2)
        f_serial_indep <- function(x) {
                cubature::adaptIntegrate(function(t) f_gen_indep(x + t[1] - t[2]) * f_inc(t[1]) * f_inc(t[2]), c(0, 0), c(Inf, Inf))$integral
                # cubature::adaptIntegrate(function(t) f_tost_indep(x - t) * f_inc(t), 0, Inf)$integral
                # pracma::integral2(function(t, t1) f_gen_indep(x + t1 - t) * f_inc(t1) * f_inc(t), 0, int_max, 0, int_max)$Q
                # pracma::quadv(function(t) f_tost_indep(x - t) * f_inc(t), 0, int_max)$Q
        }
        
        # Numerical approximation
        tic <- Sys.time()
        f_serial_indep_ti <-
                tryCatch({
                        sapply(ti, f_serial_indep)
                }, error = function(err) {
                        warning("An error occurred: ", conditionMessage(err), "\n")
                        rep(0, length(ti))                               # Return a vector of zeros
                })
        print(toc <- Sys.time() - tic)
        
        # # Numerical approximation
        # tic <- Sys.time()
        # cl <- parallel::makeCluster(parallel::detectCores() - 1)         # Create a cluster with the detected number of CPU cores
        # parallel::clusterExport(cl, c(ls(), lsf.str()))                  # Export the all variables and functions to the cluster
        # f_serial_indep_ti <-
        #         tryCatch({
        #                 parallel::parSapply(cl, ti, f_serial_indep)      # Use parSapply to apply the function in parallel
        #         }, error = function(err) {
        #                 warning("An error occurred: ", conditionMessage(err), "\n")
        #                 rep(0, length(ti))                               # Return a vector of zeros
        #         })
        # parallel::stopCluster(cl)                                        # Stop the cluster when done
        # print(toc <- Sys.time() - tic)
        
        # Checking
        if (sum(f_serial_indep_ti * dt) < 0.99) {
                warning("Difficulties with the integration of the serial interval distribution, f_serial_indep (!)")
        } else {
                # Plot the serial interval distribution
                if (plotting) {
                        plot(ti, f_serial_indep_ti,
                             xlim = c(-5, 20),
                             main = paste0("Mean = ", sprintf("%.1f", sum(f_serial_indep_ti * ti) * dt),
                                           ", SD = ", sprintf("%.1f", sqrt(sum(f_serial_indep_ti * ti^2) * dt - (sum(f_serial_indep_ti * ti) * dt)^2))))
                }
        }
        
        ## Save all three distributions ################################################
        distr_indep <- data.frame(time     = ti,
                                  f_gen    = f_gen_indep_ti,
                                  f_tost   = f_tost_indep_ti,
                                  f_serial = f_serial_indep_ti)
        # Save the numerical approximation to a file
        save(distr_indep,
             file = paste0("../R/Results-temp/", scn, "/RData/gen_tost_serial_indep.RData"))
}
