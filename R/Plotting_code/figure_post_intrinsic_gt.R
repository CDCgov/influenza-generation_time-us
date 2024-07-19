## figure1.m ###################################################################
figure_post_intrinsic_gt <- function(scn) {
  # library
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggdist))
  suppressMessages(library(dplyr))
  
  # Load results
  load(paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_indep.RData"))
  load(paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_mech.RData"))
  
  df_indep      <- data.frame(
    mean    = mcmc_posterior_indep$mean_post,
    sd      = mcmc_posterior_indep$sd_post,
    presymp = mcmc_posterior_indep$prob_presymp_post,
    beta    = mcmc_posterior_indep$beta_post,
    model   = "indep"
  )
  df_indep$id <- 1:nrow(df_indep)
  
  df_mech      <- data.frame(
    mean    = mcmc_posterior_mech$mean_post,
    sd      = mcmc_posterior_mech$sd_post,
    presymp = mcmc_posterior_mech$prob_presymp_post,
    beta    = mcmc_posterior_mech$beta_post,
    model   = "mech"
  )
  df_mech$id <- 1:nrow(df_mech)
  
  # Combine the data frames
  df <- rbind(df_indep, df_mech)
  
  # Reshape the data into long format
  df <- reshape2::melt(df, id.vars = c("model", "id"))
  
  # ensemble model
  df_ensemble <- df %>%
    bind_rows(mutate(., model = "ensemble")) %>%
    mutate(model = factor(model, levels = c("indep", "mech", "ensemble")))
  
  # Define a custom labeller function to include the mean value in the facet title
  custom_labeller <- function(variable) {
    # Mean data frame
    df_mean <- df_ensemble %>%
      group_by(variable, model) %>%
      summarize(mean_value  = sprintf("%.2f", mean(value)),
                upper_value = sprintf("%.2f", quantile(value, probs = 0.025)),
                lower_value = sprintf("%.2f", quantile(value, probs = 0.975)),
                .groups = "drop") %>%
      group_by(variable) %>%
      summarise(text = paste(model, paste0(mean_value, " (", upper_value, ", ", lower_value, ")"), sep = ": ", collapse = "\n"),
                .groups = "drop")
    
    # labeller text
    label <- paste0(variable, "\n",
                    df_mean$text[df_mean$variable == variable])
    return(label)
  }
  
  # Plotting
  p <-
    ggplot(df_ensemble, aes(x = model, y = value, fill = model)) +
    geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9)) +
    facet_wrap(~ variable, scales = "free_y",
               labeller = labeller(variable = custom_labeller)) +
    labs(x = "Model",
         y = "Posterior") +
    theme_bw(base_size = 24)
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_post_intrinsic_gt-ensemble", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 11, height = 8.5, type = "cairo")
  
  # Plotting
  p <-
    ggplot() +
    geom_line(data = df, aes(x = id, y = value, col = model),
              alpha = 0.5) +
    geom_hline(data = df %>%
                 group_by(variable, model) %>%
                 summarize(mean_value = mean(value)),
               aes(yintercept = mean_value, col = model),
               linetype = "dashed", size = 1.5) +
    facet_wrap(~ variable, scales = "free_y",
               labeller = labeller(variable = custom_labeller)) +
    labs(x = "Sample index",
         y = "Posterior") +
    theme_bw(base_size = 24)
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_post_intrinsic_gt-mcmc", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 11, height = 8.5, type = "cairo")
}
