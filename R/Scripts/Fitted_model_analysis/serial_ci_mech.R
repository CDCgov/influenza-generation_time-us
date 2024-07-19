## new script ##################################################################
cat("\014")    # clear the console
rm(list=ls())  # remove all variables
graphics.off() # close all plots
set.seed(2023) # fix randomness
suppressMessages(library(ggplot2))
suppressMessages(library(ggdist))
suppressMessages(library(dplyr))

## input parameter #################################################################
args <- commandArgs(trailingOnly = TRUE)       # Get command-line arguments
if (length(args) == 0) { args <- c(1, 2) }     # Check if arguments are provided; if not, use a default value
i <- as.numeric(args[1])                       # Convert the first argument to numeric
j <- as.numeric(args[2])                       # Convert the second argument to numeric
print(paste0("i = ", i, " (index of job)"))    # Print the index
print(paste0("j = ", j, " (index of sample)")) # Print the index

## single sample ###############################################################
serial_ci_mech <- function(i, j) {
  set.seed(2023) # fix randomness
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggdist))
  suppressMessages(library(dplyr))
  
  # Create a data frame as baseline
  param_seq_flu <- data.frame(
    virus          = "flu",
    no_steps       = 1000000,      # MCMC iterations
    period         = c("Two"),     # two seasons or single season
    coprimary      = c(FALSE),     # data with or without co-primary cases
    hh_size        = c("1p"),      # household sizes
    type           = c("Both"),    # type A or B
    inc_mean       = c(1.55),      # mean of incubation period
    inc_sd         = c(0.66),      # sd   of incubation period
    x_A            = c(0.57)       # relative infectiousness of asymptomatic infected hosts
  )
  # sensitivity analyses to store the combinations of parameters
  param_seq_flu <-
    bind_rows(param_seq_flu,
              param_seq_flu %>% slice(rep(1, 3)) %>% mutate(period         = c("Two", "2021_22", "2022_23")),     # two seasons or single season
              param_seq_flu %>% slice(rep(1, 2)) %>% mutate(coprimary      = c(TRUE, FALSE)),                     # data with or without co-primary cases
              param_seq_flu %>% slice(rep(1, 3)) %>% mutate(hh_size        = c("1p", "2or3", "4p")),              # household sizes
              param_seq_flu %>% slice(rep(1, 3)) %>% mutate(type           = c("Both", "A", "B")),                # type A or B (the same incubation period)
              param_seq_flu %>% slice(rep(1, 12)) %>% mutate(type          = rep(c("Both", "A", "B"), each = 4),  # type AND incubation period
                                                             inc_mean      = rep(c(0.61, 1.55, 4.30, 5.80), times = 3),
                                                             inc_sd        = rep(c(0.25, 0.66, 1.25, 3.10), times = 3)),
              param_seq_flu %>% slice(rep(1, 4)) %>% mutate(inc_mean       = c(0.61, 1.55, 4.30, 5.80),           # mean of incubation period
                                                            inc_sd         = c(0.25, 0.66, 1.25, 3.10)),          # sd   of incubation period
              param_seq_flu %>% slice(rep(1, 3)) %>% mutate(x_A            = c(0.11, 0.57, 1.54)),                # relative infectiousness of asymptomatic infected hosts
    ) %>%
    distinct()
  
  # Create a data frame as baseline
  param_seq_cov <- data.frame(
    virus          = "cov",
    no_steps       = 1000000,          # MCMC iterations
    period         = c("All Omicron"), # full period or by variant
    coprimary      = c(FALSE),         # data with or without co-primary cases
    hh_size        = c("1p"),          # household sizes
    type           = c("All"),         # all types
    inc_mean       = c(2.60),          # mean of incubation period
    inc_sd         = c(1.00),          # sd   of incubation period
    x_A            = c(0.35)           # relative infectiousness of asymptomatic infected hosts
  )
  # sensitivity analyses to store the combinations of parameters
  param_seq_cov <-
    bind_rows(param_seq_cov,
              param_seq_cov %>% slice(rep(1, 4)) %>% mutate(period         = c("All Omicron",
                                                                               "BA1", "BA4", "XBB")),             # full period or by variant
              param_seq_cov %>% slice(rep(1, 12)) %>% mutate(period         = rep(c("BA1", "BA4", "XBB"), each = 4), # variant AND incubation period
                                                             inc_mean       = rep(c(2.60, 2.90, 3.70, 5.80), times = 3),
                                                             inc_sd         = rep(c(1.00, 1.30, 1.60, 3.10), times = 3)),
              param_seq_cov %>% slice(rep(1, 2)) %>% mutate(coprimary      = c(TRUE, FALSE)),                     # data with or without co-primary cases
              param_seq_cov %>% slice(rep(1, 3)) %>% mutate(hh_size        = c("1p", "2or3", "4p")),              # household sizes
              param_seq_cov %>% slice(rep(1, 1)) %>% mutate(type           = c("All")),                           # all types
              param_seq_cov %>% slice(rep(1, 4)) %>% mutate(inc_mean       = c(2.60, 2.90, 3.70, 5.80),           # mean of incubation period
                                                            inc_sd         = c(1.00, 1.30, 1.60, 3.10)),          # sd   of incubation period
              param_seq_cov %>% slice(rep(1, 3)) %>% mutate(x_A            = c(0.10, 0.35, 1.27)),                # relative infectiousness of asymptomatic infected hosts
    ) %>%
    distinct()
  
  # flu and COVID-19
  param_seq <- rbind(param_seq_flu, param_seq_cov)
  param_seq$index = seq_len(nrow(param_seq))       # add index
  
  # Parameters Scripts
  source("../R/Parameters/scenario_name.R")
  scn <- scenario_name(i, param_seq)
  
  # folder
  if (j == 1) {
    # unlink(    paste0("../R/Results-temp/", scn, "/RData_SI/"), recursive = TRUE) # Delete the folder and its contents
    dir.create(paste0("../R/Results-temp/", scn, "/RData_SI/"), recursive = TRUE) # Create the new folder
  }
  
  # Load assumed parameters
  load(paste0("../R/Results-temp/", scn, "/RData/assumed_parameters.RData"))
  
  # Load output of MCMC fitting procedure
  load(paste0("../R/Results-temp/", scn, "/RData/param_fit_mech.RData"))
  
  # posterior parameter
  params_post <- data.frame(gamma = params_known[2],
                            mu    = 1 / result$theta_mat[, 2],
                            k_inc = params_known[1],
                            k_E   = params_known[1] * result$theta_mat[, 1],
                            k_I   = params_known[3],
                            alpha = result$theta_mat[, 3],
                            beta  = result$theta_mat[, 4],
                            rho   = params_known[4],
                            x_A   = params_known[5])
  
  # Define constants for the distributions
  dt      <- 0.1
  ti      <- seq(-50, 50, by = dt)
  int_max <- max(ti)
  
  ## Serial interval #############################################################
  # One sample
  params_i = j
  params_post_i <- as.numeric(params_post[params_i, ])
  
  # Serial interval
  source("../R/Functions/Mech/get_serial_dist_mech.R")
  f_serial_mech <- get_serial_dist_mech(params_post_i, int_max)
  
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
  
  # mean and SD given one sample
  SI_mean <- sum(f_serial_mech_ti * ti) * dt
  SI_sd   <- sqrt(sum(f_serial_mech_ti * ti^2) * dt - (sum(f_serial_mech_ti * ti) * dt)^2)
  
  # Print mean and SD for the current sample
  cat(sprintf("Sample %d: Mean = %.2f, SD = %.2f\n", params_i, SI_mean, SI_sd))
  
  # Save results
  df_SI <- data.frame(params_i = params_i,
                      SI_mean  = SI_mean,
                      SI_sd    = SI_sd)
  save(df_SI, file = paste0("../R/Results-temp/", scn, "/RData_SI/serial_mech_", params_i, ".RData"))
}

# run the function
serial_ci_mech(i, j)
