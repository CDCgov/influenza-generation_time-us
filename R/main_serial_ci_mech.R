## new script ##################################################################
cat("\014")    # clear the console
rm(list=ls())  # remove all variables
graphics.off() # close all plots
set.seed(2023) # fix randomness
suppressMessages(library(ggplot2))
suppressMessages(library(ggdist))
suppressMessages(library(dplyr))

## main script #################################################################

# Running one sample
source("../R/Scripts/Fitted_model_analysis/serial_ci_mech.R")

# Running samples in parallel
tic <- Sys.time()
cl <- parallel::makeCluster(parallel::detectCores(), outfile = "")                  # Create a cluster with the detected number of CPU cores
parallel::clusterExport(cl, c(ls(), lsf.str()))                                     # Export the all variables and functions to the cluster
parallel::parSapply(cl, 1:nrow(params_post), function(j) serial_ci_mech(1, j))      # Use parSapply to apply the function in parallel
parallel::stopCluster(cl)                                                           # Stop the cluster when done
print(toc <- Sys.time() - tic)

# Load each scenario in each folder
dir_folders <- list.dirs(path = "../R/Results-temp", recursive = FALSE, full.names = TRUE)
dir_folders <- dir_folders[grepl("/\\d+-", dir_folders)]
for(dir_folder in dir_folders) {
  
  # Load file from each sample
  df_SI_list <- list()
  dir_files <- list.files(path = file.path(dir_folder, "RData_SI"), pattern = "^serial_mech_.*\\.RData$", full.names = TRUE)
  for(dir_file in dir_files) {
    # load data frame
    load(dir_file)
    # check for zeros and replace with NA
    df_SI[df_SI == 0] <- NA
    # save in the list
    df_SI_list[[dir_file]] <- df_SI
  }
  
  # Combine all data frames into one
  df_SI <- data.table::rbindlist(df_SI_list)
  
  # Calculate mean and 95% CI, and format the results
  df_SI_mean <- df_SI %>%
    tidyr::gather(variable, value, -params_i) %>%
    group_by(variable) %>%
    summarize(
      Mean = mean(value),
      CI_2p5 = quantile(value, probs = 0.025),
      CI_97p5 = quantile(value, probs = 0.975)
    ) %>%
    mutate(
      Mean = sprintf("%.1f", Mean),
      CI_2p5 = sprintf("%.1f", CI_2p5),
      CI_97p5 = sprintf("%.1f", CI_97p5),
    )
  
  # Save results as a CSV file
  write.csv(df_SI_mean, file = paste0(dir_folder, "/RData_SI/df_serial_mech_", nrow(df_SI), ".csv"))
}
