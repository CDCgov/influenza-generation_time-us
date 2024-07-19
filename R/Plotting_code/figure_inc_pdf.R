## figure_inc_pdf.m ############################################################
figure_inc_pdf <- function(scn) {
  # library
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggdist))
  suppressMessages(library(dplyr))
  
  # Load incubation period distribution
  load(paste0("../R/Results-temp/", scn, "/RData/assumed_parameters.RData"))
  
  # # The mean and variance of incubation period
  # inc_mean   <- exp(inc_mu + 0.5 * inc_sigma^2)
  # inc_var    <- (exp(inc_sigma^2) - 1) * exp(2 * inc_mu + inc_sigma^2)
  
  # Define constants for the distributions
  dt      <- 0.01
  ti      <- seq(0, 15, by = dt)
  int_max <- max(ti)
  
  # Numerical approximation
  f_inc_logn_ti <- f_inc_logn(ti)
  f_inc_gam_ti  <- f_inc_gam(ti)
  
  # As data frame
  distr_inc <- data.frame(time     = ti,
                          f_inc_logn_ti = f_inc_logn_ti,
                          f_inc_gam_ti  = f_inc_gam_ti
  )
  
  # Reshape the data into long format
  df <- reshape2::melt(distr_inc, id.vars = c("time"))
  df$model <- factor(df$variable,
                     levels = c("f_inc_logn_ti", "f_inc_gam_ti"),
                     labels = c("indep"        , "mech"))
  df$variable <- "f_inc"
  
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
  
  # Plotting
  p <-
    ggplot() +
    geom_line(data = df, aes(x = time, y = value, col = model),
              size = 1, linetype = "dashed") +
    facet_wrap(~ variable, scales = "free",
               labeller = labeller(variable = custom_labeller)) +
    labs(x = "Days",
         y = "Probability Density") +
    theme_bw(base_size = 24)
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_inc_pdf", ".png") )
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
    scale_x_continuous(limits = c(0, 4)) +
    labs(x = "Days",
         y = "Probability Density") +
    theme_bw(base_size = 24) +
    theme(legend.position = "none")
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_inc_pdf-mech", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 11, height = 8.5, type = "cairo")
}
