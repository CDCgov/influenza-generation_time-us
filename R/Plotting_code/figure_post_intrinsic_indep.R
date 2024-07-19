## figure1_supp1.m #############################################################
figure_post_intrinsic_indep <- function(scn) {
  # library
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggdist))
  suppressMessages(library(dplyr))
  
  # Load results
  load(paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_indep.RData"))
  
  df_post      <- data.frame(
    mean    = mcmc_posterior_indep$mean_post,
    sd      = mcmc_posterior_indep$sd_post,
    beta    = mcmc_posterior_indep$beta_post,
    distr   = "post"
  )
  df_post$id <- 1:nrow(df_post)
  
  df_prior      <- data.frame(
    theta   = seq(0, 10, by = 0.01),
    mean    = dlnorm(seq(0, 10, by = 0.01), meanlog = 1.6, sdlog = 0.35),
    sd      = dlnorm(seq(0, 10, by = 0.01), meanlog = 0.7, sdlog = 0.65),
    beta    = dlnorm(seq(0, 10, by = 0.01), meanlog = 0.7, sdlog = 0.8),
    distr   = "prior"
  )
  
  # Reshape the data into long format
  df_post  <- reshape2::melt(df_post,  id.vars = c("distr", "id"))
  df_prior <- reshape2::melt(df_prior, id.vars = c("distr", "theta")) %>%
    filter(variable != "beta" | theta < 2.5)
  
  # Define a custom labeller function to include the mean value in the facet title
  custom_labeller <- function(variable) {
    # Mean data frame
    df_mean <- df_post %>%
      group_by(variable) %>%
      summarize(mean_value  = sprintf("%.2f", mean(value)),
                upper_value = sprintf("%.2f", quantile(value, probs = 0.025)),
                lower_value = sprintf("%.2f", quantile(value, probs = 0.975)),
                .groups = "drop") %>%
      group_by(variable) %>%
      summarise(text = paste("indep", paste0(mean_value, " (", upper_value, ", ", lower_value, ")"), sep = ": ", collapse = "\n"),
                .groups = "drop")
    
    # labeller text
    label <- paste0(variable, "\n",
                    df_mean$text[df_mean$variable == variable])
    return(label)
  }
  
  # Plotting
  p <-
    ggplot() +
    geom_histogram(data = df_post, aes(x = value, y = ..density..),
                   fill = "blue", color = "black", alpha = 0.5) +
    geom_vline(data = df_post %>%
                 group_by(variable, distr) %>%
                 summarize(mean_value = mean(value)),
               aes(xintercept = mean_value),
               linetype = "dashed", color = "red", size = 1.5) +
    geom_line(data = df_prior, aes(x = theta, y = value),
              linetype = "dashed", color = "black", size = 1) +
    facet_wrap(~ variable, scales = "free", nrow = 2,
               labeller = labeller(variable = custom_labeller)) +
    labs(x = "Theta",
         y = "Density") +
    theme_bw(base_size = 24)
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_post_intrinsic_indep-prior", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 11, height = 8.5, type = "cairo")
  
  # Plotting
  p <-
    ggplot() +
    geom_line(data = df_prior, aes(x = theta, y = value),
              linetype = "dashed", color = "black") +
    geom_ribbon(data = df_prior, aes(x = theta, ymin = 0, ymax = value),
                alpha = 0.5, fill = "gray") +
    geom_density(data = df_post, aes(x = value), adjust = 2,
                 linetype = "solid", color = "black", size = 1) +
    geom_area(data = df_post, aes(x = value), adjust = 2,
              stat = "density", position = "identity",
              alpha = 0.5, fill = "darkgray") +
    facet_wrap(variable ~ ., scales = "free", nrow = 2,
               labeller = labeller(variable = c("mean"    = "Mean generation time (days)",
                                                "sd"      = "SD of generation time (days)",
                                                "presymp" = "Proportion of transmission \n before symptom onset",
                                                "alpha"    = "Ratio of pre-symptomatic and \n symptomatic transmission rates",
                                                "beta"     = "Overall infectiousness",
                                                "presymp"  = "Proportion of transmission \n before symptom onset",
                                                "p_E"      = "Ratio of mean latent and incubation periods",
                                                "period_E" = "Mean latent period (days)",
                                                "period_P" = "Mean pre-symptomatic infectious period (days)",
                                                "mu_inv"   = "Mean symptomatic infectious period (days)"))) +
    labs(x = "Parameters",
         y = "Probability Density") +
    theme_bw(base_size = 24)
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_post_intrinsic_indep-violin", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 11, height = 8.5, type = "cairo")
  
  # Plotting
  p <-
    ggplot() +
    geom_line(data = df_post, aes(x = id, y = value),
              alpha = 0.5) +
    geom_hline(data = df_post %>%
                 group_by(variable) %>%
                 summarize(mean_value = mean(value)),
               aes(yintercept = mean_value),
               linetype = "dashed", size = 1.5) +
    facet_wrap(~ variable, scales = "free_y", nrow = 2,
               labeller = labeller(variable = custom_labeller)) +
    labs(x = "Sample index",
         y = "Posterior") +
    theme_bw(base_size = 24)
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_post_intrinsic_indep-mcmc", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 11, height = 8.5, type = "cairo")
}
