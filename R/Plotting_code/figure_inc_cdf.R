## figure_inc_cdf.m ############################################################
figure_inc_cdf <- function(scn) {
  # library
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggdist))
  suppressMessages(library(dplyr))
  
  # Load incubation period distribution
  load(paste0("../R/Results-temp/", scn, "/RData/assumed_parameters.RData"))
  
  # Define constants for the distributions
  dt      <- 0.01
  ti      <- seq(0, 5, by = dt)
  int_max <- max(ti)
  
  # Log-normal incubation period
  F_inc_logn <- function(t_inc) plnorm(t_inc, meanlog = inc_mu, sdlog = inc_sigma)
  
  # Gamma incubation period
  F_inc_gam <- function(t_inc) pgamma(t_inc, shape = inc_shape, scale = inc_scale)
  
  # Numerical approximation
  F_inc_logn_ti <- F_inc_logn(ti)
  F_inc_gam_ti  <- F_inc_gam(ti)
  
  # As data frame
  distr_inc <- data.frame(time     = ti,
                          F_inc_logn_ti = F_inc_logn_ti,
                          F_inc_gam_ti  = F_inc_gam_ti
  )
  
  # Reshape the data into long format
  df <- reshape2::melt(distr_inc, id.vars = c("time"))
  df$model <- factor(df$variable,
                     levels = c("F_inc_logn_ti", "F_inc_gam_ti"),
                     labels = c("indep"        , "mech"))
  df$variable <- "F_inc"
  
  # Define a custom labeller function to include the mean value in the facet title
  custom_labeller <- function(variable) {
    # Mean data frame
    df_mean <- df %>%
      group_by(variable, model) %>%
      summarize(mean_value  = sprintf("%.1f", inc_mean),
                sd_value    = sprintf("%.1f", inc_sd),
                .groups = "drop") %>%
      mutate(range = ifelse(model == "indep",
                            sprintf("%.1f-%.1f",
                                    qlnorm(0.05, meanlog = inc_mu, sdlog = inc_sigma),
                                    qlnorm(0.95, meanlog = inc_mu, sdlog = inc_sigma)),
                            sprintf("%.1f-%.1f",
                                    qgamma(0.05, shape = inc_shape, scale = inc_scale),
                                    qgamma(0.95, shape = inc_shape, scale = inc_scale)))) %>%
      group_by(variable) %>%
      summarise(text = paste(model, paste0("(Mean = ", mean_value, ", SD = ", sd_value, ", Range = ", range, ")"), sep = ": ", collapse = "\n"),
                .groups = "drop")
    
    # labeller text
    label <- paste0(variable, "\n",
                    df_mean$text[df_mean$variable == variable])
    return(label)
  }
  
  # Given percentiles and corresponding values from Table 3 of Lessler et al
  df_flu <- data.frame(
    percentiles = c(0.05, 0.25, 0.50, 0.75, 0.95),
    observed_values_A = c(0.7, 1.1, 1.4, 1.9, 2.8),  # Influenza A
    observed_values_B = c(0.3, 0.4, 0.6, 0.7, 1.1)  # Influenza B
  )
  
  # Plotting
  p <-
    ggplot() +
    geom_point(data = df_flu, aes(x = observed_values_A, y = percentiles),
               col = "blue", size = 5, shape = 1) +
    geom_point(data = df_flu, aes(x = observed_values_B, y = percentiles),
               col = "red", size = 5, shape = 2) +
    geom_line(data = df, aes(x = time, y = value, col = model),
              size = 1, linetype = "dashed") +
    facet_wrap(~ variable, scales = "free",
               labeller = labeller(variable = custom_labeller)) +
    labs(x = "Days",
         y = "Cumulative Probability") +
    theme_bw(base_size = 24)
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_inc_cdf", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 11, height = 8.5, type = "cairo")
  
  # Plotting (gamma only)
  p <-
    ggplot() +
    geom_line(data = df %>%
                filter(model == "mech"),
              aes(x = time, y = value, col = model),
              size = 2) +
    scale_color_manual(name = "Model",
                       values = c("indep" = "orange3",
                                  "mech"  = "steelblue3"),
                       labels = c("indep" = "Indep",
                                  "mech"  = "Mech")) +
    geom_point(data = df_flu, aes(x = observed_values_A, y = percentiles),
               col = "black", size = 5, shape = 1, stroke = 2) +
    scale_x_continuous(limits = c(0, 4)) +
    labs(x = "Days",
         y = "Cumulative Probability") +
    theme_bw(base_size = 24) +
    theme(legend.position = "none")
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_inc_cdf-mech", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 11, height = 8.5, type = "cairo")
  
  ## combine png #################################################################
  # load files
  print(file_names <- c(paste0("../R/Results-temp/", scn, "/Figures/", "figure_inc_",
                               c("pdf", "cdf"),
                               "-mech", ".png")))
  
  # all figures
  rl = lapply(file_names, png::readPNG)
  gl = lapply(rl, grid::rasterGrob)
  
  # plot
  p = cowplot::plot_grid(plotlist = gl,
                         nrow = 1,
                         labels = "AUTO",
                         hjust = 0,
                         label_size = 24, label_fontface = "plain",
                         align = "hv", axis = "b")
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_inc_2-mech", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 11*2, height = 8.5, type = "cairo")
}
