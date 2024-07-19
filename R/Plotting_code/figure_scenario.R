## figure_scenario.m ###########################################################
figure_scenario <- function(param_seq) {
  # library
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggdist))
  suppressMessages(library(dplyr))
  
  # Delete and recreate folders for Figures and Results.
  no_steps <- param_seq$no_steps[1]
  no_steps_text <- format(no_steps, scientific = FALSE)
  unlink(    paste0("../R/Results-temp/", no_steps_text, "/figure_scenario/"), recursive = TRUE) # Delete the folder and its contents
  dir.create(paste0("../R/Results-temp/", no_steps_text, "/figure_scenario/"), recursive = TRUE) # Create the new folder
  
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
      load(paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_indep.RData"))
      load(paste0("../R/Results-temp/", scn, "/RData/mcmc_posterior_mech.RData"))
      
      df_indep      <- data.frame(
        mean    = mcmc_posterior_indep$mean_post,
        sd      = mcmc_posterior_indep$sd_post,
        presymp = mcmc_posterior_indep$prob_presymp_post,
        beta    = mcmc_posterior_indep$beta_post,
        # mean_hh    = mcmc_posterior_indep$empirical_summary_mat[, 1],
        # sd_hh      = mcmc_posterior_indep$empirical_summary_mat[, 2],
        # presymp_hh = mcmc_posterior_indep$empirical_summary_mat[, 3],
        ll_vec  = mcmc_posterior_indep$ll_vec,
        model   = "indep"
      )
      df_indep$id <- 1:nrow(df_indep)
      
      df_mech      <- data.frame(
        mean    = mcmc_posterior_mech$mean_post,
        sd      = mcmc_posterior_mech$sd_post,
        presymp = mcmc_posterior_mech$prob_presymp_post,
        beta    = mcmc_posterior_mech$beta_post,
        # mean_hh    = mcmc_posterior_mech$empirical_summary_mat[, 1],
        # sd_hh      = mcmc_posterior_mech$empirical_summary_mat[, 2],
        # presymp_hh = mcmc_posterior_mech$empirical_summary_mat[, 3],
        ll_vec  = mcmc_posterior_mech$ll_vec,
        model   = "mech"
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
        ll_vec  = NA,
        model   = "indep"
      )
      df_indep$id <- 1:nrow(df_indep)
      
      df_mech      <- data.frame(
        mean    = mcmc_posterior_mech$empirical_summary_mat[, 1],
        sd      = mcmc_posterior_mech$empirical_summary_mat[, 2],
        presymp = mcmc_posterior_mech$empirical_summary_mat[, 3],
        beta    = NA,
        ll_vec  = NA,
        model   = "mech"
      )
      df_mech$id <- 1:nrow(df_mech)
      
      # Combine the data frames (realized)
      df_realized <- rbind(df_indep, df_mech)
      df_realized$gen <- "Realized"
      
      # Combine the data frames (intrinsic and realized)
      df <- rbind(df_intrinsic, df_realized)
      
      # Reshape the data into long format
      df <- reshape2::melt(df, id.vars = c("model", "id", "gen"))
      
      # ensemble
      df_scenario[[i]] <- df %>%
        # bind_rows(mutate(., model = "ensemble")) %>% # add ensemble model
        mutate(model = factor(model, levels = c("indep", "mech", "ensemble")),
               index     = index,
               virus     = virus,
               period    = period,
               coprimary = coprimary,
               hh_size   = hh_size,
               type      = type,
               no_steps = no_steps,
               inc_mean = inc_mean,
               inc_sd   = inc_sd,
               x_A      = x_A)
      
    }, error = function(err) {
      warning("An error occurred: ", conditionMessage(err), "\n")
    })
  }
  
  # Combine the data frames of all scenarios
  df_scenarios <- data.table::rbindlist(df_scenario)
  
  ## (sensitivity analysis) plotting GT only ###################################
  # plot geom_violin (flu; data stratifications)
  p <-
    ggplot(df_scenarios %>%
             mutate(type_inc_mean = paste0(type, "\n(", inc_mean, ", ", inc_sd, ")")) %>%
             mutate(scenario = factor(case_when(
               period != "Two" ~ "period",
               type != "Both" ~ "type",
               hh_size != "1p" ~ "hh_size",
               coprimary != FALSE ~ "coprimary",
               TRUE ~ "baseline"  # Default case
             ), levels = c("baseline", "period", "type", "hh_size", "coprimary")),
             scenario_label = case_when(
               scenario == "period" ~ as.character(period),
               scenario == "type" ~ as.character(type_inc_mean),
               scenario == "hh_size" ~ as.character(hh_size),
               scenario == "coprimary" ~ as.character(coprimary),
               scenario == "baseline" ~ as.character(period)
             )) %>%
             filter(virus == "flu",
                    model == "mech",
                    variable == "mean",
                    gen == "Intrinsic",
                    (abs(inc_mean - 1.55) < 0.001 & type != "B") |
                      (abs(inc_mean - 0.61) < 0.001 & type == "B"),
                    abs(x_A - 0.57) < 0.001),
           aes(x = scenario_label, y = value, fill = model)) +
    geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9), adjust = 2) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
    scale_fill_manual(name = "Model",
                      values = c("indep" = "orange3",
                                 "mech"  = "steelblue3"),
                      labels = c("indep" = "Indep",
                                 "mech"  = "Mech")) +
    facet_grid(. ~ scenario,
               scales = "free",
               labeller = labeller(scenario = c("baseline" = "Baseline",
                                                "period" = "Season",
                                                "coprimary" = "Co-primary cases",
                                                "hh_size" = "Household sizes",
                                                "type" = "Type",
                                                "type_inc_mean" = "Type (mean and standard deviation of incubation period)",
                                                "inc_mean" = "Incubation period (mean and standard deviation)",
                                                "inc_mean_typeA" = "Incubation period (mean and standard deviation)",
                                                "inc_mean_typeB" = "Incubation period (mean and standard deviation)",
                                                "x_A" = "Relative infectiousness of asymptomatic individuals"))) +
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
                                "BA1 (2.9, 1.3)" = "BA1/2 (2.9, 1.3)",
                                "BA4" = "BA4/5")) +
    labs(x = "Scenarios",
         y = "Posterior of mean intrinsic generation time (days)") +
    theme_bw(base_size = 22) +
    theme(#axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none")
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", no_steps_text, "/figure_scenario/", "figure_scenario_violin_", "flu", "_1_GT", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16, height = 9, type = "cairo")
  
  # data frame copied from the above ggplot (to calculate overlap)
  df_scenarios_gt <- df_scenarios %>%
    mutate(type_inc_mean = paste0(type, "\n(", inc_mean, ", ", inc_sd, ")")) %>%
    mutate(scenario = factor(case_when(
      period != "Two" ~ "period",
      type != "Both" ~ "type",
      hh_size != "1p" ~ "hh_size",
      coprimary != FALSE ~ "coprimary",
      TRUE ~ "baseline"  # Default case
    ), levels = c("baseline", "period", "type", "hh_size", "coprimary")),
    scenario_label = case_when(
      scenario == "period" ~ as.character(period),
      scenario == "type" ~ as.character(type_inc_mean),
      scenario == "hh_size" ~ as.character(hh_size),
      scenario == "coprimary" ~ as.character(coprimary),
      scenario == "baseline" ~ as.character(period)
    )) %>%
    filter(virus == "flu",
           model == "mech",
           variable == "mean",
           gen == "Intrinsic",
           (abs(inc_mean - 1.55) < 0.001 & type != "B") |
             (abs(inc_mean - 0.61) < 0.001 & type == "B"),
           abs(x_A - 0.57) < 0.001)
  
  # as list
  df_scenarios_gt <- df_scenarios_gt %>%
    group_by(scenario_label) %>%
    summarise(values = list(value))
  # matrix of overlapping index
  overlap_matrix <- matrix(NA,
                           nrow = nrow(df_scenarios_gt),
                           ncol = nrow(df_scenarios_gt),
                           dimnames = list(df_scenarios_gt$scenario_label, df_scenarios_gt$scenario_label))
  for (i in 1:nrow(df_scenarios_gt)) {
    for (j in 1:nrow(df_scenarios_gt)) {
      # overlapping index
      overlap_matrix[i, j] <- overlapping::overlap(data.frame(df_scenarios_gt$values[[i]],
                                                              df_scenarios_gt$values[[j]]))$OV
    }
  }
  print(overlap_matrix)
  
  # plot geom_violin (flu; incubation period)
  p <-
    ggplot(df_scenarios %>%
             filter(virus == "flu",
                    period    == "Two",
                    coprimary == FALSE,
                    hh_size == "1p",
                    type == "Both",
                    model == "mech",
                    variable == "mean",
                    gen == "Intrinsic",
                    inc_mean < 5.8,
                    abs(x_A - 0.57) < 0.001),
           aes(x = as.factor(inc_mean), y = value, fill = model)) +
    geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9), adjust = 2) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
    scale_fill_manual(name = "Model",
                      values = c("indep" = "orange3",
                                 "mech"  = "steelblue3"),
                      labels = c("indep" = "Indep",
                                 "mech"  = "Mech")) +
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
                                "BA1 (2.9, 1.3)" = "BA1/2 (2.9, 1.3)",
                                "BA4" = "BA4/5")) +
    labs(x = "Incubation period (mean and standard deviation)",
         y = "Posterior of mean intrinsic generation time (days)") +
    theme_bw(base_size = 22) +
    theme(#axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none")
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", no_steps_text, "/figure_scenario/", "figure_scenario_violin_", "flu", "_2_GT", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16, height = 9, type = "cairo")
  
  # plot geom_violin (cov)
  p <-
    ggplot(df_scenarios %>%
             filter(virus == "cov",
                    model == "mech",
                    variable == "mean",
                    gen == "Intrinsic",
                    coprimary == FALSE,
                    hh_size == "1p",
                    abs(inc_mean - 2.60) < 0.001,
                    abs(x_A - 0.35) < 0.001),
           aes(x = period, y = value, fill = model)) +
    geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9), adjust = 2) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
    scale_fill_manual(name = "Model",
                      values = c("indep" = "orange3",
                                 "mech"  = "steelblue3"),
                      labels = c("indep" = "Indep",
                                 "mech"  = "Mech")) +
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
                                "BA1 (2.9, 1.3)" = "BA1/2 (2.9, 1.3)",
                                "BA4" = "BA4/5")) +
    labs(x = "Period",
         y = "Posterior of mean intrinsic generation time (days)") +
    theme_bw(base_size = 26) +
    theme(#axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none")
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", no_steps_text, "/figure_scenario/", "figure_scenario_violin_", "cov", "_1_GT", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16, height = 9, type = "cairo")
  
  # plot geom_violin (cov)
  p <-
    ggplot(df_scenarios %>%
             filter(virus == "cov",
                    model == "mech",
                    variable == "beta",
                    gen == "Intrinsic",
                    coprimary == FALSE,
                    hh_size == "1p",
                    abs(inc_mean - 2.60) < 0.001,
                    abs(x_A - 0.35) < 0.001),
           aes(x = period, y = value, fill = model)) +
    geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9), adjust = 2) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
    scale_fill_manual(name = "Model",
                      values = c("indep" = "orange3",
                                 "mech"  = "steelblue3"),
                      labels = c("indep" = "Indep",
                                 "mech"  = "Mech")) +
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
                                "BA1 (2.9, 1.3)" = "BA1/2 (2.9, 1.3)",
                                "BA4" = "BA4/5")) +
    labs(x = "Period",
         y = "Overall infectiousness") +
    theme_bw(base_size = 26) +
    theme(#axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none")
  
  # save figure
  print( name_plot <- paste0("../R/Results-temp/", no_steps_text, "/figure_scenario/", "figure_scenario_violin_", "cov", "_2_beta", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16, height = 9, type = "cairo")
  
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
        group_by(variable, model, gen, index, virus, period, coprimary, hh_size, type, inc_mean, inc_sd, x_A) %>%
        summarise(mean_value  = sprintf("%.1f", mean(value, na.rm = T)),
                  upper_value = sprintf("%.1f", quantile(value, probs = 0.025, na.rm = T)),
                  lower_value = sprintf("%.1f", quantile(value, probs = 0.975, na.rm = T)))
      # save csv
      print( name_csv <- paste0("../R/Results-temp/", no_steps_text, "/figure_scenario/", "df_scenarios_subset_table_", virus, "_", target, ".csv") )
      write.csv(df_scenarios_subset_table, name_csv, row.names = FALSE)
      
      # plot geom_violin
      p <-
        ggplot(df_scenarios_subset %>%
                 filter(variable != "ll_vec") %>%
                 filter(variable != "beta" | gen == "Intrinsic") %>%
                 filter(inc_mean < 5.8),
               aes(x = as.factor(get(target)), y = value, fill = model)) +
        geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9), adjust = 2) +
        geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
        scale_fill_manual(name = "Model",
                          values = c("indep" = "orange3",
                                     "mech"  = "steelblue3"),
                          labels = c("indep" = "Indep",
                                     "mech"  = "Mech")) +
        facet_wrap(gen ~ variable, scales = "free", nrow = 2,
                   labeller = labeller(variable = c("mean"    = "Mean generation time (days)",
                                                    "sd"      = "SD of generation time (days)",
                                                    "presymp" = "Proportion of transmission \n before symptom onset",
                                                    "beta"    = "Overall infectiousness"))) +
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
                                    "BA1 (2.9, 1.3)" = "BA1/2 (2.9, 1.3)",
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
        theme_bw(base_size = 18) +
        theme(#axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
      
      # save figure
      print( name_plot <- paste0("../R/Results-temp/", no_steps_text, "/figure_scenario/", "figure_scenario_violin_", virus, "_", target, ".png") )
      ggplot2::ggsave(file = name_plot, plot = p, width = 16, height = 9, type = "cairo")
      
      # plot geom_violin (mech only)
      p <-
        ggplot(df_scenarios_subset %>%
                 filter(model == "mech") %>%
                 filter(variable != "ll_vec") %>%
                 filter(variable != "beta" | gen == "Intrinsic") %>%
                 filter(inc_mean < 5.8),
               aes(x = as.factor(get(target)), y = value, fill = model)) +
        geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9), adjust = 2) +
        geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
        scale_fill_manual(name = "Model",
                          values = c("indep" = "orange3",
                                     "mech"  = "steelblue3"),
                          labels = c("indep" = "Indep",
                                     "mech"  = "Mech")) +
        facet_wrap(variable ~ gen, scales = "free", ncol = 2,
                   labeller = labeller(variable = c("mean"    = "Mean generation time (days)",
                                                    "sd"      = "SD of generation time (days)",
                                                    "presymp" = "Proportion of transmission \n before symptom onset",
                                                    "beta"    = "Overall infectiousness"))) +
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
      print( name_plot <- paste0("../R/Results-temp/", no_steps_text, "/figure_scenario/", "figure_scenario_violin_", virus, "_", target, "-mech.png") )
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
                                   "/figure_scenario/", "figure_scenario_violin_", virus,
                                   "_", targets, "-mech.png")))
      
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
      print( name_plot <- paste0("../R/Results-temp/", no_steps_text, "/figure_scenario/", "figure_scenario_violin_", virus, "_", tg, ".png") )
      ggplot2::ggsave(file = name_plot, plot = p, width = 13*2, height = 16*2, type = "cairo")
    }
  }
  
  ## combine png #################################################################
  unlink(    paste0("../R/Results-temp/", no_steps_text, "/figure_1ABCD/"), recursive = TRUE) # Delete the folder and its contents
  dir.create(paste0("../R/Results-temp/", no_steps_text, "/figure_1ABCD/"), recursive = TRUE) # Create the new folder
  # load files
  print(file_names <- c(paste0("../R/Results-temp/", "1-flu-Two-FALSE-hh_1p-Both-1000000-1p55-0p66-0p57",
                               "/Figures/",
                               c("figure_gen_serial-mech", "figure_gt-mech-mean"), ".png"),
                        paste0("../R/Results-temp/", no_steps_text,
                               "/figure_scenario/", "figure_scenario_violin_", "flu",
                               "_", c(1:2), "_GT", ".png")))
  
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
  print( name_plot <- paste0("../R/Results-temp/", no_steps_text, "/figure_1ABCD/", "figure1", ".png") )
  ggplot2::ggsave(file = name_plot, plot = p, width = 16*2, height = 9*2, type = "cairo")
}
