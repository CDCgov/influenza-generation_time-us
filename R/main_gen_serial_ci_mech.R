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
source("../R/Scripts/Fitted_model_analysis/gen_serial_ci_mech.R")

# Running samples in parallel
tic <- Sys.time()
cl <- parallel::makeCluster(parallel::detectCores(), outfile = "")                  # Create a cluster with the detected number of CPU cores
parallel::clusterExport(cl, c(ls(), lsf.str()))                                     # Export the all variables and functions to the cluster
parallel::parSapply(cl, 1:nrow(params_post), function(j) gen_serial_ci_mech(1, j))  # Use parSapply to apply the function in parallel
parallel::stopCluster(cl)                                                           # Stop the cluster when done
print(toc <- Sys.time() - tic)

# Load each scenario in each folder
dir_folders <- list.dirs(path = "../R/Results-temp", recursive = FALSE, full.names = TRUE)
dir_folders <- dir_folders[grepl("/\\d+-", dir_folders)]
for(dir_folder in dir_folders) {
  
  # Define constants for the distributions
  dt      <- 0.1
  ti      <- seq(-50, 50, by = dt)
  int_max <- max(ti)
  
  # Load file from each sample of generation time
  df_gen_list <- list()
  dir_files <- list.files(path = file.path(dir_folder, "RData_gen_serial"), pattern = "^gen_mech_.*\\.RData$", full.names = TRUE)
  for(dir_file in dir_files) {
    # load data frame
    load(dir_file)
    # check for zeros and replace with NA
    if (all(f_gen_mech_ti == 0)) f_gen_mech_ti[] <- NA
    # save in the list
    df_gen_list[[dir_file]] <- data.frame(Days         = ti,
                                          Density      = f_gen_mech_ti,
                                          Samples      = dir_file,
                                          Distribution = "Generation time")
  }
  
  # Load file from each sample of serial interval
  df_serial_list <- list()
  dir_files <- list.files(path = file.path(dir_folder, "RData_gen_serial"), pattern = "^serial_mech_.*\\.RData$", full.names = TRUE)
  for(dir_file in dir_files) {
    # load data frame
    load(dir_file)
    # check for zeros and replace with NA
    if (all(f_serial_mech_ti == 0)) f_serial_mech_ti[] <- NA
    # save in the list
    df_serial_list[[dir_file]] <- data.frame(Days         = ti,
                                             Density      = f_serial_mech_ti,
                                             Samples      = dir_file,
                                             Distribution = "Serial interval")
  }
  
  # Combine all data frames into one
  df_gen_mech_ti    <- data.table::rbindlist(df_gen_list)
  df_serial_mech_ti <- data.table::rbindlist(df_serial_list)
  # Both distributions
  df_dist_mech_ti   <- rbind(df_gen_mech_ti, df_serial_mech_ti)
  # Save the distributions to a file
  save(df_dist_mech_ti,
       file = paste0(dir_folder, "/RData/gen_serial_mech_ci.RData"))
  
  # Calculate mean and 95% CI of mean and SD over all samples
  df_dist_mech_csv <- df_dist_mech_ti %>%
    group_by(Distribution, Samples) %>%
    summarise(Mean = sum(Density * Days) * dt,
              SD   = sqrt(sum(Density * Days^2) * dt - (sum(Density * Days) * dt)^2)) %>%
    group_by(Distribution) %>%
    summarise(Mean_mean = mean(Mean),
              Mean_q025 = quantile(Mean, probs = 0.025),
              Mean_q975 = quantile(Mean, probs = 0.975),
              SD_mean = mean(SD),
              SD_q025 = quantile(SD, probs = 0.025),
              SD_q975 = quantile(SD, probs = 0.975)) %>%
    mutate(Mean_mean = sprintf("%.1f", Mean_mean),
           Mean_q025 = sprintf("%.1f", Mean_q025),
           Mean_q975 = sprintf("%.1f", Mean_q975),
           SD_mean = sprintf("%.1f", SD_mean),
           SD_q025 = sprintf("%.1f", SD_q025),
           SD_q975 = sprintf("%.1f", SD_q975)) %>%
    ungroup()
  # Save summary as a CSV file
  write.csv(df_dist_mech_csv,
            file = paste0(dir_folder, "/RData_gen_serial/gen_serial_mech_ci", ".csv"))
  
  # Calculate mean of density over all samples
  df_dist_mech_ti_csv = df_dist_mech_ti %>%
    group_by(Distribution, Days) %>%
    summarise(Density = mean(Density)) %>%
    ungroup() %>%
    mutate(Days = round(Days, 1)) # Round Days to 1 decimal place
  # Save summary as a CSV file
  write.csv(df_dist_mech_ti_csv,
            file = paste0(dir_folder, "/RData_gen_serial/gen_serial_mech_density_mean", ".csv"))
  
  # plot
  p <- df_dist_mech_ti %>%
    filter(
      (Distribution == "Generation time" & Days >= 0 & Days <= 10) |
        (Distribution == "Serial interval" & Days >= -1 & Days <= 10)
    ) %>%
    ggplot() +
    tidybayes::stat_lineribbon(
      aes(x        = Days,
          y        = Density,
          col      = Distribution,
          fill     = Distribution,
          linetype = Distribution),
      size = 1,
      alpha = 0.25,
      .width = c(.95)) +
    stat_summary(
      aes(x        = Days,
          y        = Density,
          col      = Distribution,
          linetype = Distribution),
      fun = median,
      geom = "line",
      size = 1
    ) +
    scale_color_manual(name = "Distribution",
                       values = c("Generation time" = "steelblue3",
                                  "Serial interval"  = "orange3")) +
    scale_fill_manual(name = "Distribution",
                      values = c("Generation time" = "steelblue3",
                                 "Serial interval"  = "orange3")) +
    scale_linetype_manual(name = "Distribution",
                          values = c("Generation time" = "solid",
                                     "Serial interval" = "dashed")) +
    scale_x_continuous(limits = c(-1, 10)) +
    labs(x = "Days",
         y = "Density") +
    theme_bw(base_size = 24) +
    theme(legend.position = c(0.8, 0.8),
          legend.key.width = unit(10, "line"))
  # save figure
  print( name_plot <- paste0(dir_folder, "/Figures/", "figure_gen_serial_ci-mech", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16, height = 9, type = "cairo")
}
