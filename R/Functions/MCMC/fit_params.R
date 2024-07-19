## fit_params.m ################################################################
fit_params <- function(scn, no_steps, steps_keep, update_theta, update_infection, update_onset, update_asymp, empirical_summary_form, theta_init, data_struct_augmented_init, ll_household_init, plotting) {
  # Use data augmentation MCMC to fit the parameters, theta, of the model
  # of infectiousness under consideration, to the household transmission
  # data.
  
  # library
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggdist))
  suppressMessages(library(dplyr))
  
  no_params_fitted <- length(theta_init)
  
  # Initialize the vector of fitted parameters, theta, and the structure
  # array containing the augmented data, data_struct_augmented
  theta                 <- theta_init
  data_struct_augmented <- data_struct_augmented_init
  
  # Calculate initial likelihood contributions from each household
  ll_household <- ll_household_init
  
  # Vectors to hold output of fitting procedure (at steps at
  # which this is kept after burn-in and thinning)
  no_steps_kept         <- length(steps_keep)
  theta_mat             <- matrix(NA, nrow = no_steps_kept, ncol = no_params_fitted) # parameters at each step
  ll_vec                <- rep(NA, no_steps_kept)                                    # likelihood at each step
  empirical_summary_mat <- matrix(NA, nrow = no_steps_kept, ncol = 3)                # estimated mean/sd of realized household generation times, and the presymptomatic proportion of transmissions
  
  # Vectors to hold output of fitting procedure at every step
  theta_mat_all <- matrix(NA, nrow = no_steps, ncol = no_params_fitted)
  ll_vec_all    <- rep(NA, no_steps)
  
  # Vectors to hold information about which steps are accepted
  acceptance_vec                 <- rep(NA, no_steps)
  acceptance_vec_theta           <- rep(NA, no_steps)
  acceptance_vec_infection_symp  <- rep(NA, no_steps)
  acceptance_vec_onset           <- rep(NA, no_steps)
  acceptance_vec_asymp           <- rep(NA, no_steps) # updating just infection times of asymptomatic hosts in the independent transmission and symptoms model, but both infection times of asymptomatic hosts and times of entry into the I stage in the mechanistic model
  acceptance_vec_infection_asymp <- rep(NA, no_steps) # updating just infection times of asymptomatic hosts
  
  # Break steps into 100 groups to record progress
  step_no_mat <- matrix(1:no_steps, nrow = no_steps / 100, ncol = 100)
  
  # If plotting enabled, set up a figure to plot output
  if (plotting) {
    # Create an empty plot with burn-in period in shaded area
    p0 <- ggplot() +
      geom_rect(data = data.frame(xmin = 0, xmax = steps_keep[1] - 1, ymin = -Inf, ymax = Inf),
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "darkgray", alpha = 0.2) +
      facet_grid(variable ~ ., scales = "free_y") +
      labs(x = "Sample index",
           y = "Posterior") +
      theme_bw(base_size = 18)
  }
  
  # Run chain
  step_no_kept <- 0
  
  for (j in 1:ncol(step_no_mat)) {
    for (i in 1:nrow(step_no_mat)) {
      
      step_no <- step_no_mat[i, j]
      
      # Carry out one of 4 different steps
      if (step_no %% 4 < 0.5) {
        # Update the model parameters, theta
        result                        <- update_theta(theta, data_struct_augmented, ll_household)
        theta                         <- result$theta
        ll_household                  <- result$ll_household
        acceptance_vec_theta[step_no] <- result$acceptance$overall
        
      } else if (step_no %% 4 < 1.5) {
        # Update times of infection (for symptomatic infected hosts
        # in the independent transmission and symptoms model, and
        # for all infected hosts in the mechanistic model)
        result                                  <- update_infection(theta, data_struct_augmented, ll_household)
        data_struct_augmented                   <- result$data_struct_augmented
        ll_household                            <- result$ll_household
        acceptance_vec_infection_symp[step_no]  <- result$acceptance$symp
        acceptance_vec_infection_asymp[step_no] <- result$acceptance$asymp
        
      } else if (step_no %% 4 < 2.5) {
        # Update times of symptom onset for hosts who developed symptoms
        result                        <- update_onset(theta, data_struct_augmented, ll_household)
        data_struct_augmented         <- result$data_struct_augmented
        ll_household                  <- result$ll_household
        acceptance_vec_onset[step_no] <- result$acceptance$overall
        
      } else {
        # Update augmented data for asymptomatic infected hosts
        # (just infection times for the independent transmission
        # and symptoms model, times of both infection and entry
        # into the I stage for the mechanistic model)
        result                                  <- update_asymp(theta, data_struct_augmented, ll_household)
        data_struct_augmented                   <- result$data_struct_augmented
        ll_household                            <- result$ll_household
        acceptance_vec_asymp[step_no]           <- result$acceptance$overall
        acceptance_vec_infection_asymp[step_no] <- result$acceptance$infection
      }
      
      # Record parameters/likelihood/acceptance from the current step
      theta_mat_all[step_no,] <- theta
      ll_vec_all[step_no]     <- sum(ll_household)
      acceptance_vec[step_no] <- result$acceptance$overall
      
      # If step kept after burn-in and thinning, record further output
      if (step_no %in% steps_keep) {
        step_no_kept                         <- step_no_kept + 1
        theta_mat[step_no_kept,]             <- theta
        ll_vec[step_no_kept]                 <- sum(ll_household)
        empirical_summary                    <- empirical_summary_form(theta, data_struct_augmented)
        empirical_summary_mat[step_no_kept,] <- empirical_summary
      }
    }
    
    # Display progress when an integer percentage of steps has been completed
    print(sprintf("%3d%% completed (%s)", 100 * step_no / no_steps, Sys.time()))
    
    # If plotting enabled, plot output up to current step
    if (plotting) {
      # Plotting all theta
      p <- p0 +
        geom_line(data = reshape2::melt(data.frame(id            = 1:nrow(theta_mat_all),
                                                   theta_mat_all = theta_mat_all,
                                                   ll_vec_all    = ll_vec_all), id.vars = "id"),
                  aes(x = id, y = value),
                  alpha = 0.5, col = "black")
      
      if (step_no >= steps_keep[1]) {
        # Plotting thinned theta
        p <- p +
          geom_point(data = reshape2::melt(data.frame(id            = steps_keep,
                                                      theta_mat_all = theta_mat,
                                                      ll_vec_all    = ll_vec), id.vars = "id"),
                     aes(x = id, y = value),
                     alpha = 0.2, col = "black")
      }
      
      # save figure
      print( name_plot <- paste0("../R/Results-temp/flu-", scn, "/Figures/", "figure1_supp", ifelse(ncol(theta_mat_all) == 3, "1", "2"), "-mcmc-all", ".png") )
      ggplot2::ggsave(file = name_plot, plot = p, width = 11, height = 8.5, type = "cairo")
    }
  }
  
  # Calculate acceptance rates: overall, and for the different update steps
  acceptance_vec_theta           <- acceptance_vec_theta[!is.na(acceptance_vec_theta)]
  acceptance_vec_infection_symp  <- acceptance_vec_infection_symp[!is.na(acceptance_vec_infection_symp)]
  acceptance_vec_asymp           <- acceptance_vec_asymp[!is.na(acceptance_vec_asymp)]
  acceptance_vec_onset           <- acceptance_vec_onset[!is.na(acceptance_vec_onset)]
  acceptance_vec_infection_asymp <- acceptance_vec_infection_asymp[!is.na(acceptance_vec_infection_asymp)]
  
  acceptance_rate_overall         <- mean(acceptance_vec)
  acceptance_rate_theta           <- mean(acceptance_vec_theta)
  acceptance_rate_infection_symp  <- mean(acceptance_vec_infection_symp)
  acceptance_rate_onset           <- mean(acceptance_vec_onset)
  acceptance_rate_asymp           <- mean(acceptance_vec_asymp)
  acceptance_rate_infection_asymp <- mean(acceptance_vec_infection_asymp)
  
  acceptance_rate <- list(
    overall         = acceptance_rate_overall,
    theta           = acceptance_rate_theta,
    infection_symp  = acceptance_rate_infection_symp,
    onset           = acceptance_rate_onset,
    asymp           = acceptance_rate_asymp,
    infection_asymp = acceptance_rate_infection_asymp
  )
  
  # Structure array containing additional output from the fitting procedure
  output <- list(
    acceptance_rate             = acceptance_rate,
    data_struct_augmented_final = data_struct_augmented,
    empirical_summary_mat       = empirical_summary_mat
  )
  
  return(list(theta_mat = theta_mat,
              ll_vec    = ll_vec,
              output    = output))
}
