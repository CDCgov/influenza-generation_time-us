## new script ##################################################################
cat("\014")    # clear the console
rm(list=ls())  # remove all variables
graphics.off() # close all plots
set.seed(2023) # fix randomness
suppressMessages(library(ggplot2))
suppressMessages(library(ggdist))
suppressMessages(library(dplyr))

## main script #################################################################

# Running one scenario
source("../R/main_scenario.R")

# Running scenarios in parallel
tic <- Sys.time()
cl <- parallel::makeCluster(parallel::detectCores(), outfile = "")                  # Create a cluster with the detected number of CPU cores
parallel::clusterExport(cl, c(ls(), lsf.str()))                                     # Export the all variables and functions to the cluster
parallel::parSapply(cl, 1:nrow(param_seq), main_scenario)                           # Use parSapply to apply the function in parallel
parallel::stopCluster(cl)                                                           # Stop the cluster when done
print(toc <- Sys.time() - tic)

# Plotting estimates of all scenarios
source("../R/Plotting_code/figure_scenario.R")
figure_scenario(param_seq)

# Plotting parameters of all scenarios (mech model only)
source("../R/Plotting_code/figure_scenario_mech.R")
figure_scenario_mech(param_seq)
