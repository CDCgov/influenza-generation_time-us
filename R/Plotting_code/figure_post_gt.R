## figure1_supp5.m #############################################################
figure_post_gt <- function(scn) {
  # library
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggdist))
  suppressMessages(library(dplyr))
  
  # Load results
  load(paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_indep.RData"))
  load(paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_mech.RData"))
  
  # intrinsic
  df_indep      <- data.frame(
    mean    = mcmc_posterior_indep$mean_post,
    sd      = mcmc_posterior_indep$sd_post,
    presymp = mcmc_posterior_indep$prob_presymp_post,
    beta    = mcmc_posterior_indep$beta_post,
    model   = "Indep"
  )
  df_indep$id <- 1:nrow(df_indep)
  
  df_mech      <- data.frame(
    mean    = mcmc_posterior_mech$mean_post,
    sd      = mcmc_posterior_mech$sd_post,
    presymp = mcmc_posterior_mech$prob_presymp_post,
    beta    = mcmc_posterior_mech$beta_post,
    model   = "Mech"
  )
  df_mech$id <- 1:nrow(df_mech)
  
  # Combine the data frames (intrinsic)
  df_intrinsic <- rbind(df_indep, df_mech)
  df_intrinsic$gen <- "Intrinsic"
  
  # realized
  df_indep      <- data.frame(
    mean    = mcmc_posterior_indep$empirical_summary_mat[, 1],
    sd      = mcmc_posterior_indep$empirical_summary_mat[, 2],
    presymp = mcmc_posterior_indep$empirical_summary_mat[, 3],
    beta    = NA,
    model   = "Indep"
  )
  df_indep$id <- 1:nrow(df_indep)
  
  df_mech      <- data.frame(
    mean    = mcmc_posterior_mech$empirical_summary_mat[, 1],
    sd      = mcmc_posterior_mech$empirical_summary_mat[, 2],
    presymp = mcmc_posterior_mech$empirical_summary_mat[, 3],
    beta    = NA,
    model   = "Mech"
  )
  df_mech$id <- 1:nrow(df_mech)
  
  # Combine the data frames (realized)
  df_realized <- rbind(df_indep, df_mech)
  df_realized$gen <- "Realized"
  
  # Combine the data frames (intrinsic and realized)
  df <- rbind(df_intrinsic, df_realized)
  
  # Reshape the data into long format
  df <- reshape2::melt(df, id.vars = c("model", "id", "gen"))
  
  # ensemble model
  df_ensemble <- df %>%
    bind_rows(mutate(., model = "Ensemble")) %>%
    mutate(model = factor(model, levels = c("Indep", "Mech", "Ensemble")))
  
  # Plotting
  p <-
    ggplot(df %>%
             filter(variable == "mean" | variable == "sd"),
           aes(x = gen, y = value, fill = model)) +
    geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9), adjust = 2) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
    facet_wrap(variable ~ ., scales = "free_y",
               labeller = labeller(variable = c("mean"    = "Mean generation time (days)",
                                                "sd"      = "SD of generation time (days)",
                                                "presymp" = "Proportion of transmission before symptom onset",
                                                "beta"    = "Overall infectiousness"))
    ) +
    scale_fill_manual(name = "Model",
                      values = c("Indep" = "orange3",
                                 "Mech"  = "steelblue3"),
                      labels = c("Indep" = "Indep",
                                 "Mech"  = "Mech")) +
    labs(x = "Generation time",
         y = "Posterior",
         fill = "Model") +
    theme_bw(base_size = 15) +
    theme(# axis.text.x = element_blank(),
      # axis.title.x = element_blank(),
      legend.position = "bottom")
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_gt", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16, height = 9, type = "cairo")
  
  ## (mech only) ###############################################################
  # Plotting
  p <-
    ggplot(df_ensemble %>%
             filter(variable == "mean" | variable == "sd") %>%
             filter(model == "Mech"),
           aes(x = gen, y = value)) +
    geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9), adjust = 2,
                fill = "steelblue3") +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA,
                 fill = "steelblue3") +
    facet_wrap(variable ~ ., scales = "free_y",
               labeller = labeller(variable = c("mean"    = "Mean generation time (days)",
                                                "sd"      = "SD of generation time (days)",
                                                "presymp" = "Proportion of transmission before symptom onset",
                                                "beta"    = "Overall infectiousness"))
    ) +
    labs(x = "Generation time",
         y = "Posterior",
         fill = "Model") +
    theme_bw(base_size = 15)
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_gt-mech", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16, height = 9, type = "cairo")
  
  # Filter the data for the "Mech" model and "mean" variable
  df_mech <- df_ensemble %>%
    filter(model == "Mech" & variable %in% c("mean"))
  # Separate the data into "Intrinsic" and "Realized" groups
  df_mech_intrinsic <- df_mech %>%
    filter(gen == "Intrinsic")
  df_mech_realized <- df_mech %>%
    filter(gen == "Realized")
  
  # Perform paired t-test
  t_test_result <- t.test(df_mech_intrinsic$value, df_mech_realized$value, paired = TRUE)
  # Format and print the results
  p_value     <- sprintf("%.2f", t_test_result$p.value)
  mean_value  <- sprintf("%.2f", t_test_result$estimate)
  ci_lower    <- sprintf("%.2f", t_test_result$conf.int[1])
  ci_upper    <- sprintf("%.2f", t_test_result$conf.int[2])
  result_string <- paste0("diff = ", mean_value, " (", ci_lower, ", ", ci_upper, ") with ", "p-value = ", p_value)
  print(result_string)
  
  # Plotting (mech-mean only)
  p <-
    ggplot(df_mech,
           aes(x = gen, y = value)) +
    geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9), adjust = 2,
                fill = "steelblue3") +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA,
                 fill = "steelblue3") +
    labs(x = "Generation time",
         y = "Posterior of mean generation time (days)") +
    theme_bw(base_size = 24)
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", scn, "/Figures/", "figure_gt-mech-mean", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16, height = 9, type = "cairo")
}
