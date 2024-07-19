## figure_gen_tost_serial.m ####################################################
figure_gen_tost_serial <- function(scn) {
  # library
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggdist))
  suppressMessages(library(dplyr))
  
  # Load results
  load(paste0("../R/Results-temp/", scn, "/RData/gen_tost_serial_indep.RData"))
  load(paste0("../R/Results-temp/", scn, "/RData/gen_tost_serial_mech.RData"))
  
  distr_indep$model <- "indep"
  distr_mech$model  <- "mech"
  
  # Combine the data frames
  df <- rbind(distr_indep, distr_mech)
  
  # Reshape the data into long format
  df <- reshape2::melt(df, id.vars = c("model", "time"))
  
  # Define a custom labeller function to include the mean value in the facet title
  custom_labeller <- function(variable) {
    # Mean data frame
    df_mean <- df %>%
      group_by(variable, model) %>%
      summarize(mean_value  = sprintf("%.1f", sum(value * time) * diff(df$time)[1]),
                sd_value    = sprintf("%.1f", sqrt(sum(value * time^2) * diff(df$time)[1] - (sum(value * time) * diff(df$time)[1])^2)),
                # sd_value    = sprintf("%.1f", sqrt(sum((time - sum(value * time) * diff(df$time)[1])^2 * value) * diff(df$time)[1])),
                .groups = "drop") %>%
      group_by(variable) %>%
      summarise(text = paste(model, paste0("(Mean = ", mean_value, ", SD = ", sd_value, ")"), sep = ": ", collapse = "\n"),
                .groups = "drop")
    
    # labeller text
    label <- paste0(variable, "\n",
                    df_mean$text[df_mean$variable == variable])
    return(label)
  }
  
  # Range of plotting
  df_1 <- df %>%
    filter((variable == "f_gen" & time >= 0 & time <= 15) |
             (variable == "f_tost" & time >= -10 & time <= 10) |
             (variable == "f_serial" & time >= -5 & time <= 20))
  
  # Discontinuous line for plotting TOST distribution
  df_0 <- df_1
  df_0$value[df_0$model == "mech" & df_0$variable == "f_tost" & df_0$time == 0] <- NA
  
  # Plotting
  p <-
    ggplot() +
    geom_line(data = df_1, aes(x = time, y = value, col = model),
              size = 1, linetype = "dashed") +
    geom_line(data = df_0, aes(x = time, y = value, col = model),
              size = 2) +
    facet_wrap(~ variable, scales = "free",
               labeller = labeller(variable = custom_labeller)) +
    labs(x = "Days",
         y = "Density") +
    theme_bw(base_size = 24)
  
  ## Currently used in Rt estimation #############################################
  # Load data
  load(paste0("../R/Results-temp/", scn, "/RData/data.RData"))
  
  # Check if the data is from US flu data with 838 households
  if (length(data_struct_observed$household_no) == 838) {
    
    # Define the parameters for the current gamma distribution
    gen_mean  <- 3.6
    gen_var   <- 1.6^2
    gen_shape <- gen_mean^2 / gen_var
    gen_scale <- gen_var / gen_mean
    
    # Create a new data frame with the gamma distribution
    df_gen_gam <- df_1 %>%
      filter(variable == "f_gen" & model == "indep") %>%
      mutate(model = "gamma",
             value = dgamma(time, shape = gen_shape, scale = gen_scale))
    
    # Add the gamma distribution to the plot
    p <- p +
      geom_line(data = df_gen_gam, aes(x = time, y = value),
                size = 1, linetype = "dashed", col = "darkgray")
  }
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_gen_tost_serial", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16*2, height = 9, type = "cairo")
  
  ## without TOST ##############################################################
  # Plotting
  p <-
    ggplot() +
    geom_line(data = df_1 %>%
                filter(variable != "f_tost"),
              aes(x = time, y = value, col = model, linetype = variable),
              size = 2) +
    scale_color_manual(name = "Model",
                       values = c("indep" = "orange3",
                                  "mech"  = "steelblue3"),
                       labels = c("indep" = "Indep",
                                  "mech"  = "Mech")) +
    scale_linetype_manual(name = "Distribution",
                          values = c("f_gen"    = "solid",
                                     "f_serial" = "dashed"),
                          labels = c("f_gen"    = "Generation time",
                                     "f_serial" = "Serial interval")) +
    scale_x_continuous(limits = c(-1, 10)) +
    labs(x = "Days",
         y = "Density") +
    theme_bw(base_size = 24) +
    theme(legend.position = c(0.8, 0.8),
          legend.key.width = unit(10, "line"))
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_gen_serial", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16, height = 9, type = "cairo")
  
  ## without TOST (mech only) ##################################################
  # Plotting
  p <-
    ggplot() +
    geom_line(data = df_1 %>%
                filter(model == "mech") %>%
                filter(variable != "f_tost"),
              aes(x = time, y = value, linetype = variable),
              col = "steelblue3",
              size = 2) +
    scale_linetype_manual(name = "Distribution",
                          values = c("f_gen"    = "solid",
                                     "f_serial" = "dashed"),
                          labels = c("f_gen"    = "Generation time",
                                     "f_serial" = "Serial interval")) +
    scale_x_continuous(limits = c(-1, 10)) +
    labs(x = "Days",
         y = "Density") +
    theme_bw(base_size = 24) +
    theme(legend.position = c(0.8, 0.8),
          legend.key.width = unit(10, "line"))
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_gen_serial-mech", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16, height = 9, type = "cairo")
}
