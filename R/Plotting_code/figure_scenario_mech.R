## figure_scenario_mech.m ######################################################
figure_scenario_mech <- function(param_seq) {
  # library
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggdist))
  suppressMessages(library(dplyr))
  
  # Delete and recreate folders for Figures and Results.
  no_steps <- param_seq$no_steps[1]
  no_steps_text <- format(no_steps, scientific = FALSE)
  unlink(    paste0("../R/Results-temp/", no_steps_text, "/figure_scenario_mech/"), recursive = TRUE) # Delete the folder and its contents
  dir.create(paste0("../R/Results-temp/", no_steps_text, "/figure_scenario_mech/"), recursive = TRUE) # Create the new folder
  
  # scenarios
  df_scenario <- list()
  for(i in 1:nrow(param_seq)) {
    # Assumed parameters
    index     <- param_seq$index[i]
    virus     <- param_seq$virus[i]
    period    <- param_seq$period[i]
    coprimary <- param_seq$coprimary[i]
    hh_size   <- param_seq$hh_size[i]
    type      <- param_seq$type[i]
    no_steps <- param_seq$no_steps[i]
    inc_mean <- param_seq$inc_mean[i]
    inc_sd   <- param_seq$inc_sd[i]
    x_A      <- param_seq$x_A[i]
    
    # scenario name
    source("../R/Parameters/scenario_name.R")
    scn <- scenario_name(i, param_seq)
    
    tryCatch({
      # Load results
      load(paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_mech.RData"))
      
      df_post      <- data.frame(
        alpha   = mcmc_posterior_mech$alpha_post,
        # beta    = mcmc_posterior_mech$beta_post,
        presymp = mcmc_posterior_mech$prob_presymp_post,
        p_E     = mcmc_posterior_mech$p_E_post,
        period_E = inc_mean *      mcmc_posterior_mech$p_E_post,
        period_P = inc_mean * (1 - mcmc_posterior_mech$p_E_post),
        mu_inv  = mcmc_posterior_mech$mu_inv_post,
        model   = "mech"
      )
      df_post$id <- 1:nrow(df_post)
      
      # Reshape the data into long format
      df <- reshape2::melt(df_post, id.vars = c("model", "id"))
      
      # add columns
      df_scenario[[i]] <- df %>%
        mutate(no_steps = no_steps,
               index     = index,
               virus     = virus,
               period    = period,
               coprimary = coprimary,
               hh_size   = hh_size,
               type      = type,
               inc_mean = inc_mean,
               inc_sd   = inc_sd,
               x_A      = x_A)
      
    }, error = function(err) {
      warning("An error occurred: ", conditionMessage(err), "\n")
    })
  }
  
  # Combine the data frames of all scenarios
  df_scenarios <- data.table::rbindlist(df_scenario)
  
  ## (sensitivity analysis) plotting distribution ##############################
  # each parameter
  targets = c("period", "type", "hh_size", "coprimary",
              "inc_mean", "x_A",
              "type_inc_mean",
              "inc_mean_typeA", "inc_mean_typeB")
  for(target in targets){
    for(virus in c("flu", "cov")){
      # reset data frame
      df_scenarios_subset <- df_scenarios %>%
        mutate(type_inc_mean = paste0(type, "\n(", inc_mean, ", ", inc_sd, ")")) %>%
        mutate(inc_mean_typeA = inc_mean) %>%
        mutate(inc_mean_typeB = inc_mean)
      
      # flu subset
      if(virus == "flu"){                df_scenarios_subset = subset(df_scenarios_subset, virus     == "flu")
      if(!grepl("period",       target)) df_scenarios_subset = subset(df_scenarios_subset, period    == "Two")
      if(!grepl("coprimary",    target)) df_scenarios_subset = subset(df_scenarios_subset, coprimary == FALSE)
      if(!grepl("hh_size",      target)) df_scenarios_subset = subset(df_scenarios_subset, hh_size == "1p")
      if(!grepl("type",         target)) df_scenarios_subset = subset(df_scenarios_subset, type == "Both")
      if(!grepl("inc_mean",     target)) df_scenarios_subset = subset(df_scenarios_subset, abs(inc_mean - 1.55) < 0.001)
      if(!grepl("x_A",          target)) df_scenarios_subset = subset(df_scenarios_subset, abs(x_A      - 0.57) < 0.001)
      
      if(grepl("typeA",         target)) df_scenarios_subset = subset(df_scenarios_subset, type == "A")
      if(grepl("typeB",         target)) df_scenarios_subset = subset(df_scenarios_subset, type == "B")
      if(grepl("type_inc_mean", target)) df_scenarios_subset = subset(df_scenarios_subset, type_inc_mean %in% c("Both\n(1.55, 0.66)", "A\n(1.55, 0.66)", "B\n(0.61, 0.25)"))
      }
      
      # cov subset
      if(virus == "cov"){         df_scenarios_subset = subset(df_scenarios_subset, virus     == "cov")
      if(target != "period")      df_scenarios_subset = subset(df_scenarios_subset, period    == "All Omicron")
      if(target != "coprimary")   df_scenarios_subset = subset(df_scenarios_subset, coprimary == FALSE)
      if(target != "hh_size")     df_scenarios_subset = subset(df_scenarios_subset, hh_size == "1p")
      if(target != "type")        df_scenarios_subset = subset(df_scenarios_subset, type == "All")
      if(target != "inc_mean")    df_scenarios_subset = subset(df_scenarios_subset, abs(inc_mean - 2.60) < 0.001)
      if(target != "x_A")         df_scenarios_subset = subset(df_scenarios_subset, abs(x_A      - 0.35) < 0.001)
      }
      
      # print table
      df_scenarios_subset_table <-
        df_scenarios_subset %>%
        group_by(variable, index, virus, period, coprimary, hh_size, type, inc_mean, inc_sd, x_A) %>%
        summarise(mean_value  = sprintf("%.2f", mean(value)),
                  upper_value = sprintf("%.2f", quantile(value, probs = 0.025)),
                  lower_value = sprintf("%.2f", quantile(value, probs = 0.975)))
      
      # save csv
      print( name_csv <- paste0("../R/Results-temp/", no_steps_text, "/figure_scenario_mech/", "df_scenarios_subset_table_", virus, "_", target, ".csv") )
      write.csv(df_scenarios_subset_table, name_csv, row.names = FALSE)
      
      # plot geom_violin
      p <-
        ggplot(df_scenarios_subset %>%
                 filter(inc_mean < 5.8),
               aes(x = as.factor(get(target)), y = value, fill = model, group = get(target))) +
        geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9), adjust = 2) +
        geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
        scale_fill_manual(name = "Model",
                          values = c("indep" = "orange3",
                                     "mech"  = "steelblue3"),
                          labels = c("indep" = "Indep",
                                     "mech"  = "Mech")) +
        facet_wrap(variable ~ ., scales = "free", ncol = 2,
                   labeller = labeller(variable = c("alpha"    = "Ratio of pre-symptomatic and \n symptomatic transmission rates",
                                                    "beta"     = "Overall infectiousness",
                                                    "presymp"  = "Proportion of transmission \n before symptom onset",
                                                    "p_E"      = "Ratio of mean latent and \n incubation periods",
                                                    "period_E" = "Mean latent period (days)",
                                                    "period_P" = "Mean pre-symptomatic \n infectious period (days)",
                                                    "mu_inv"   = "Mean symptomatic \n infectious period (days)"))) +
        scale_x_discrete(labels = c("TRUE" = "With",
                                    "FALSE" = "Without",
                                    "1p" = " All",
                                    "2or3" = " 2/3",
                                    "4p" = " 4+",
                                    "5p" = "5+",
                                    "2.6" = "(2.6, 1.0)",
                                    "2.9" = "(2.9, 1.3)",
                                    "3.7" = "(3.7, 1.6)",
                                    "5.8" = "(5.8, 3.1)",
                                    "0.61" = "(0.61, 0.25)",
                                    "1.55" = "(1.55, 0.66)",
                                    "4.3"  = "(4.30, 1.25)",
                                    "Two"     = "Both",
                                    "2021_22" = "2021/2022",
                                    "2022_23" = "2022/2023",
                                    "BA1" = "BA1/2",
                                    "BA4" = "BA4/5")) +
        labs(x = c("period" = "Season",
                   "coprimary" = "Co-primary cases",
                   "hh_size" = "Household sizes",
                   "type" = "Type",
                   "type_inc_mean" = "Type (mean and standard deviation of incubation period)",
                   "inc_mean" = "Incubation period (mean and standard deviation)",
                   "inc_mean_typeA" = "Incubation period (mean and standard deviation)",
                   "inc_mean_typeB" = "Incubation period (mean and standard deviation)",
                   "x_A" = "Relative infectiousness of asymptomatic individuals")[target],
             y = "Posterior") +
        theme_bw(base_size = 28) +
        theme(#axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
      
      # save figure
      print( name_plot <- paste0("../R/Results-temp/", no_steps_text, "/figure_scenario_mech/", "figure_scenario_mech_violin_", virus, "_", target, ".png") )
      ggplot2::ggsave(file = name_plot, plot = p, width = 13, height = 16, type = "cairo")
    }
  }
  
  ## combine png #################################################################
  for(tg in c(1, 2)){
    for(virus in c("flu", "cov")){
      # parameters
      if(tg == 1) targets = c("period", "type_inc_mean", "hh_size", "coprimary")
      if(tg == 2) targets = c("inc_mean", "inc_mean_typeA", "inc_mean_typeB", "x_A")
      
      # load files
      print(file_names <- c(paste0("../R/Results-temp/", no_steps_text,
                                   "/figure_scenario_mech/", "figure_scenario_mech_violin_", virus,
                                   "_", targets, ".png")))
      
      # all figures
      rl = lapply(file_names, png::readPNG)
      gl = lapply(rl, grid::rasterGrob)
      
      # plot
      p = cowplot::plot_grid(plotlist = gl,
                             nrow = 2,
                             labels = "AUTO",
                             hjust = 0,
                             label_size = 28, label_fontface = "plain",
                             align = "hv", axis = "b")
      
      # save figure
      print( name_plot <- paste0("../R/Results-temp/", no_steps_text, "/figure_scenario_mech/", "figure_scenario_mech_violin_", virus, "_", tg, ".png") )
      ggplot2::ggsave(file = name_plot, plot = p, width = 13*2, height = 16*2, type = "cairo")
    }
  }
}
