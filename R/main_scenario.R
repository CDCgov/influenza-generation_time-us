## new script ##################################################################
cat("\014")    # clear the console
rm(list=ls())  # remove all variables
graphics.off() # close all plots
set.seed(2023) # fix randomness
suppressMessages(library(ggplot2))
suppressMessages(library(ggdist))
suppressMessages(library(dplyr))

## input parameter #################################################################
args <- commandArgs(trailingOnly = TRUE)    # Get command-line arguments
if (length(args) == 0) { args <- c(1) }     # Check if arguments are provided; if not, use a default value
i <- as.numeric(args[1])                    # Convert the first argument to numeric
print(paste0("i = ", i, " (index of job)")) # Print the index

## single scenario ###############################################################
main_scenario <- function(i) {
  set.seed(2023) # fix randomness
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggdist))
  suppressMessages(library(dplyr))
  
  # Create a data frame as baseline
  param_seq_flu <- data.frame(
    virus          = "flu",
    no_steps       = 1000,         # MCMC iterations (here, 1000 is for demo; 1000000 was used in the paper)
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
    no_steps       = 1000,             # MCMC iterations (here, 1000 is for demo; 1000000 was used in the paper)
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
  source("../R/Parameters/scenario_folder.R")
  scenario_folder(scn)
  source("../R/Parameters/assumed_parameters.R")
  assumed_parameters(i, param_seq, scn)
  # Plotting Scripts
  source("../R/Plotting_code/figure_inc_pdf.R")
  figure_inc_pdf(scn)
  source("../R/Plotting_code/figure_inc_cdf.R")
  figure_inc_cdf(scn)
  
  # Data Script to load data (the household data is not available on github)
  if(param_seq$virus[i] == "flu") source("../R/Data/import_data.R") # replaced by the demo data from https://doi.org/10.7554/elife.70767
  if(param_seq$virus[i] == "cov") source("../R/Data/import_data.R") # replaced by the demo data from https://doi.org/10.7554/elife.70767
  import_data(i, param_seq, scn)
  # Data Script to format data
  source("../R/Data/format_data.R")
  format_data(i, param_seq, scn)
  
  # MCMC Scripts
  source("../R/Scripts/Parameter_fitting/param_fit_indep.R")
  param_fit_indep(i, param_seq, scn)
  source("../R/Scripts/Parameter_fitting/param_fit_mech.R")
  param_fit_mech(i, param_seq, scn)
  
  # Posterior Scripts
  source("../R/Scripts/Fitted_model_analysis/mcmc_posterior_indep.R")
  mcmc_posterior_indep(scn)
  source("../R/Scripts/Fitted_model_analysis/mcmc_posterior_mech.R")
  mcmc_posterior_mech(scn)
  # Plotting Scripts
  source("../R/Plotting_code/figure_post_intrinsic_gt.R")
  figure_post_intrinsic_gt(scn)
  source("../R/Plotting_code/figure_post_intrinsic_indep.R")
  figure_post_intrinsic_indep(scn)
  source("../R/Plotting_code/figure_post_intrinsic_mech.R")
  figure_post_intrinsic_mech(scn)
  source("../R/Plotting_code/figure_post_realized_gt.R")
  figure_post_realized_gt(scn)
  source("../R/Plotting_code/figure_post_gt.R")
  figure_post_gt(scn)
  
  # Distribution Scripts
  source("../R/Scripts/Fitted_model_analysis/gen_tost_serial_indep.R")
  gen_tost_serial_indep(scn)
  source("../R/Scripts/Fitted_model_analysis/gen_tost_serial_mech.R")
  gen_tost_serial_mech(scn)
  # Plotting Scripts
  source("../R/Plotting_code/figure_gen_tost_serial.R")
  figure_gen_tost_serial(scn)
}

# Run the scenario
main_scenario(i)
