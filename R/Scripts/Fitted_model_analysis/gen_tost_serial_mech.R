## gen_tost_serial_mech.m ######################################################
gen_tost_serial_mech <- function(scn) {
  # Calculate distributions of the generation time, time from onset of
  # symptoms to transmission (TOST) and serial interval for the mechanistic
  # model, using point estimate (posterior mean) parameter values
  
  # This script requires the Chebfun package to run (freely available at
  # https://www.chebfun.org/download/)
  
  # Define constants for the distributions
  dt      <- 0.1
  ti      <- seq(-50, 50, by = dt)
  int_max <- max(ti)
  
  plotting <- FALSE
  
  # Load point estimate parameter values
  load(paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_mech.RData"))
  params_best <- mcmc_posterior_mech$params_best
  
  ## Generation time #############################################################
  source("../R/Functions/Mech/get_gen_dist_mech.R")
  f_gen_mech <- get_gen_dist_mech(params_best, int_max)
  
  # Numerical approximation
  tic <- Sys.time()
  f_gen_mech_ti <-
    tryCatch({
      sapply(ti, f_gen_mech)
    }, error = function(err) {
      warning("An error occurred: ", conditionMessage(err), "\n")
      rep(0, length(ti))                                           # Return a vector of zeros
    })
  print(toc <- Sys.time() - tic)
  
  # # Numerical approximation
  # tic <- Sys.time()
  # cl <- parallel::makeCluster(parallel::detectCores() - 1)         # Create a cluster with the detected number of CPU cores
  # parallel::clusterExport(cl, c(ls(), lsf.str()))                  # Export the all variables and functions to the cluster
  # f_gen_mech_ti <-
  #   tryCatch({
  #     parallel::parSapply(cl, ti, f_gen_mech)                      # Use parSapply to apply the function in parallel
  #   }, error = function(err) {
  #     warning("An error occurred: ", conditionMessage(err), "\n")
  #     rep(0, length(ti))                                           # Return a vector of zeros
  #   })
  # parallel::stopCluster(cl)                                        # Stop the cluster when done
  # print(toc <- Sys.time() - tic)
  
  # Checking
  if (sum(f_gen_mech_ti * dt) < 0.99) {
    warning("Difficulties with the integration of the generation time distribution, f_gen_mech (!)")
  } else {
    # Plot the generation time distribution
    if (plotting) {
      plot(ti, f_gen_mech_ti,
           xlim = c(0, 15),
           main = paste0("Mean = ", sprintf("%.1f", sum(f_gen_mech_ti * ti) * dt),
                         ", SD = ", sprintf("%.1f", sqrt(sum(f_gen_mech_ti * ti^2) * dt - (sum(f_gen_mech_ti * ti) * dt)^2))))
    }
  }
  
  ## TOST ########################################################################
  source("../R/Functions/Mech/f_tost_form_mech.R")
  f_tost_mech <- function(t) f_tost_form_mech(params_best, t)
  
  # Numerical approximation
  f_tost_mech_ti <- f_tost_mech(ti)
  
  # Checking
  if (sum(f_tost_mech_ti * dt) < 0.99) {
    warning("Difficulties with the integration of the TOST distribution, f_tost_mech (!)")
  } else {
    # Plot the TOST distribution
    if (plotting) {
      plot(ti, f_tost_mech_ti,
           xlim = c(-10, 10),
           main = paste0("Mean = ", sprintf("%.1f", sum(f_tost_mech_ti * ti) * dt),
                         ", SD = ", sprintf("%.1f", sqrt(sum(f_tost_mech_ti * ti^2) * dt - (sum(f_tost_mech_ti * ti) * dt)^2))))
    }
  }
  
  ## Serial interval #############################################################
  source("../R/Functions/Mech/get_serial_dist_mech.R")
  f_serial_mech <- get_serial_dist_mech(params_best, int_max)
  
  # Numerical approximation
  tic <- Sys.time()
  f_serial_mech_ti <-
    tryCatch({
      sapply(ti, f_serial_mech)
    }, error = function(err) {
      warning("An error occurred: ", conditionMessage(err), "\n")
      rep(0, length(ti))                                           # Return a vector of zeros
    })
  print(toc <- Sys.time() - tic)
  
  # # Numerical approximation
  # tic <- Sys.time()
  # cl <- parallel::makeCluster(parallel::detectCores() - 1)         # Create a cluster with the detected number of CPU cores
  # parallel::clusterExport(cl, c(ls(), lsf.str()))                  # Export the all variables and functions to the cluster
  # f_serial_mech_ti <- parallel::parSapply(cl, ti, f_serial_mech)   # Use parSapply to apply the function in parallel
  # parallel::stopCluster(cl)                                        # Stop the cluster when done
  # print(toc <- Sys.time() - tic)
  
  # Checking
  if (sum(f_serial_mech_ti * dt) < 0.99) {
    warning("Difficulties with the integration of the serial interval distribution, f_serial_mech (!)")
  } else {
    # Plot the serial interval distribution
    if (plotting) {
      plot(ti, f_serial_mech_ti,
           xlim = c(-5, 20),
           main = paste0("Mean = ", sprintf("%.1f", sum(f_serial_mech_ti * ti) * dt),
                         ", SD = ", sprintf("%.1f", sqrt(sum(f_serial_mech_ti * ti^2) * dt - (sum(f_serial_mech_ti * ti) * dt)^2))))
    }
  }
  
  ## Save all three distributions ################################################
  distr_mech <- data.frame(time     = ti,
                           f_gen    = f_gen_mech_ti,
                           f_tost   = f_tost_mech_ti,
                           f_serial = f_serial_mech_ti)
  # Save the numerical approximation to a file
  save(distr_mech,
       file = paste0("../R/Results-temp/", scn, "/RData/gen_tost_serial_mech.RData"))
}
