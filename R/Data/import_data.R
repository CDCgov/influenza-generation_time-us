## import_data.m ###############################################################
import_data <- function(i, param_seq, scn) {
  
  # Assumed data set
  virus     <- param_seq$virus[i]
  period    <- param_seq$period[i]
  coprimary <- param_seq$coprimary[i]
  hh_size   <- param_seq$hh_size[i]
  type      <- param_seq$type[i]
  
  ## import_data.m ###############################################################
  # Import household transmission data from the file
  # "Supplementary_Data.xlsx", and apply preliminary formatting steps.
  
  # Import data
  df <- readxl::read_excel('../R/Data/Supplementary_Data.xlsx')
  
  # Extract the household number, household size, and symptom onset date for
  # each host (day zero taken to be date on which index case swabbed)
  df$household_no_all   <- df$hoconumber      # household id
  df$household_size_all <- df$householdsize   # number in household
  df$d_s_all            <- df$case_swab_ill   # days from index case being swabbed to the date if illness (cases and contacts)
  df$month_all          <- df$month_case_swab # month the case was swabbed and found positive to be recruited into the study
  
  # Logical vectors indicating whether each host
  # infected/uninfected/inconclusive, asymptomatic/symptomatic
  # (will also show a 1 if symptoms even if uninfected
  # but for all intervals below the date of first symptoms was set to missing
  # if someone had symptoms but was uninfected)
  df$inf_all          <- (df$infected == 1)                # 1 = PCR or antibody testing positive
  df$uninf_all        <- (df$infected == 0)                # 0 = no evidence but antibody testing done
  df$inconclusive_all <- !(df$inf_all | df$uninf_all)      # 9 = PCR neg or missing and antibody testing not done
  df$symp_all         <- (df$inf_all & (df$symptoms == 1)) # symptomatic
  df$asymp_all        <- (df$inf_all & (df$symptoms == 0)) # asymptomatic
  
  # Uninfected or asymptomatic individuals considered to develop symptoms at time infinity
  df$d_s_all[df$asymp_all] <- Inf
  df$d_s_all[df$uninf_all] <- Inf
  
  # Find index cases
  df$status_all <- df$status            # CASE = Index, Contact = Contact
  df$index_all <- (df$status == "CASE") # CASE = Index
  
  # Obtain right bounds for infection date of each host
  df$d_iR_all            <- rep(Inf, length(df$d_s_all))
  df$d_iR_all[df$symp_all]  <- df$d_s_all[df$symp_all]             # cannot be infected after day of onset
  df$d_iR_all[df$index_all] <- pmin(df$d_iR_all[df$index_all], 0)  # since date where index swabbed positive taken as time zero
  
  # swab results
  # 0 = neg, 1 = pos, NA = not done
  df$swab1_positive <- (df$swab1 == 1 & !is.na(df$swab1))
  df$swab2_positive <- (df$swab2 == 1 & !is.na(df$swab2))
  
  # Obtain right bounds for infection date of each host
  # using the days from index case being swabbed to the first/second study swab
  df$d_iR_all[df$swab1_positive] <- pmin(df$d_iR_all[df$swab1_positive], df$case_swab_swab1[df$swab1_positive])  # must have been infected before positive test
  df$d_iR_all[df$swab2_positive] <- pmin(df$d_iR_all[df$swab2_positive], df$case_swab_swab2[df$swab2_positive])  # must have been infected before positive test
  
  # Discard data from hosts not known to have been infected or not
  df$keep <- (df$symp_all | df$asymp_all | df$uninf_all) # exclude those (df$infected == 9)
  
  df_keep <- data.frame(
    household_no_all   = df$household_no_all[df$keep],    # household id
    household_size_all = df$household_size_all[df$keep],  # number in household (includes thrown out hosts)
    month_all          = df$month_all[df$keep],           # month the case was swabbed and found positive to be recruited into the study
    d_s_all            = df$d_s_all[df$keep],             # days from index case being swabbed to the date if illness (cases and contacts)
    d_iR_all           = df$d_iR_all[df$keep],            # right bounds for infection date of each host
    symp_all           = df$symp_all[df$keep],            # symptomatic
    asymp_all          = df$asymp_all[df$keep],           # asymptomatic
    uninf_all          = df$uninf_all[df$keep],           # uninfected
    inf_all            = !df$uninf_all[df$keep],          # infected
    index_all          = df$index_all[df$keep]            # index case
  )
  
  # Order data in each household by symptom onset date, also putting asymptomatic hosts before uninfected hosts
  df_keep$order <- order(df_keep$household_no_all, df_keep$uninf_all, df_keep$d_s_all)
  
  df_order <- data.frame(
    household_no_all   = df_keep$household_no_all[df_keep$order],
    household_size_all = df_keep$household_size_all[df_keep$order],
    month_all          = df_keep$month_all[df_keep$order],
    d_s_all            = df_keep$d_s_all[df_keep$order],
    d_iR_all           = df_keep$d_iR_all[df_keep$order],
    symp_all           = df_keep$symp_all[df_keep$order],
    asymp_all          = df_keep$asymp_all[df_keep$order],
    uninf_all          = df_keep$uninf_all[df_keep$order],
    inf_all            = df_keep$inf_all[df_keep$order],
    index_all          = df_keep$index_all[df_keep$order]
  )
  
  # Households kept once data from asymptomatic/inconclusive hosts thrown out
  households_vec_all <- unique(df_order$household_no_all)   # household id
  no_households_all  <- length(households_vec_all) # number in household
  
  # Assumed maximum possible gap between successive symptom onset dates within a household
  max_onsetdiff <- 28
  
  # Separate out clusters when the maximum onset difference is exceeded
  no_households_incl            <- 0
  households_incl_old           <- c()
  households_incl_new           <- c()
  household_no_new_incl         <- c()
  household_sizes_old_incl      <- c()  # counting thrown out hosts of unknown status
  household_sizes_new_incl      <- c()  # not counting thrown out hosts
  household_size_indiv_old_incl <- c()
  household_size_indiv_new_incl <- c()
  household_months_incl         <- c()
  d_s_incl                      <- c()
  d_iR_incl                     <- c()
  symp_incl                     <- c()
  asymp_incl                    <- c()
  uninf_incl                    <- c()
  onset_diffs_incl              <- c()
  
  for (household in households_vec_all) {
    # Extract indices, symptom onset dates, and status of all individuals in household
    household_inds  <- which(df_order$household_no_all == household) # individual id
    household_d_s   <- df_order$d_s_all[household_inds]              # days from index case being swabbed to the date if illness (cases and contacts)
    household_d_iR  <- df_order$d_iR_all[household_inds]             # right bounds for infection date of each host
    household_symp  <- df_order$symp_all[household_inds]             # symptomatic
    household_asymp <- df_order$asymp_all[household_inds]            # asymptomatic
    household_uninf <- df_order$uninf_all[household_inds]            # uninfected
    
    # Household size, number who developed symptoms, month in which index case was recruited
    no_in_household_old  <- df_order$household_size_all[household_inds[1]]  # counting thrown out hosts of unknown status
    no_in_household      <- length(household_inds)                          # not counting thrown out hosts
    no_symp_in_household <- sum(household_symp)
    month_household      <- df_order$month_all[household_inds[1]]
    
    # Only keep the data from the household if
    # (i) it contains more than one household member of known infection status, and
    # (ii) each gap between successive symptom onset dates is no larger than the permitted maximum
    household_onset_diffs    <- diff(household_d_s[household_symp])
    household_onset_diff_max <- ifelse(length(household_onset_diffs) == 0, NA, max(household_onset_diffs))
    
    if ((no_in_household > 1) && (is.na(household_onset_diff_max) || household_onset_diff_max <= max_onsetdiff)) {
      no_households_incl            <- no_households_incl + 1                                                        # number of households
      households_incl_old           <- c(households_incl_old, household)                                             # household id (old)
      households_incl_new           <- c(households_incl_new, no_households_incl)                                    # household id (new)
      household_no_new_incl         <- c(household_no_new_incl, rep(no_households_incl, no_in_household))            # household id (new)
      household_sizes_old_incl      <- c(household_sizes_old_incl, no_in_household_old)                              # number in household (old)
      household_sizes_new_incl      <- c(household_sizes_new_incl, no_in_household)                                  # number in household (new)
      household_size_indiv_old_incl <- c(household_size_indiv_old_incl, rep(no_in_household_old, no_in_household))   # number in household (old)
      household_size_indiv_new_incl <- c(household_size_indiv_new_incl, rep(no_in_household, no_in_household))       # number in household (new)
      household_months_incl         <- c(household_months_incl, rep(month_household, no_in_household))               # month the case was swabbed and found positive to be recruited into the study
      d_s_incl                      <- c(d_s_incl, household_d_s)                                                    # days from index case being swabbed to the date if illness (cases and contacts)
      d_iR_incl                     <- c(d_iR_incl, household_d_iR)                                                  # right bounds for infection date of each host
      symp_incl                     <- c(symp_incl, household_symp)                                                  # symptomatic
      asymp_incl                    <- c(asymp_incl, household_asymp)                                                # asymptomatic
      uninf_incl                    <- c(uninf_incl, household_uninf)                                                # uninfected
      onset_diffs_incl              <- c(onset_diffs_incl, household_onset_diffs)                                    # ???
    } else {
      print(household_inds) # discarded individual id
    }
  }
  
  # Convert vectors indicating infection/symptom status to logical vectors
  symp_incl  <- as.logical(symp_incl)
  asymp_incl <- as.logical(asymp_incl)
  uninf_incl <- as.logical(uninf_incl)
  
  # Data of each individual
  df_data_individual <- data.frame(
    # numeric vector
    household_no_new_incl,
    household_size_indiv_old_incl,
    household_size_indiv_new_incl,
    # logical vectors
    d_s_incl,
    d_iR_incl,
    symp_incl,
    asymp_incl,
    uninf_incl
  ) %>%
    # specific household sizes
    filter(case_when(
      hh_size == "2" ~ household_size_indiv_new_incl == "2",
      hh_size == "3" ~ household_size_indiv_new_incl == "3",
      hh_size == "4" ~ household_size_indiv_new_incl == "4",
      hh_size == "5p" ~ as.integer(household_size_indiv_new_incl) >= 5,
      hh_size == "2or3" ~ household_size_indiv_new_incl %in% c("2", "3"),
      hh_size == "4p" ~ as.integer(household_size_indiv_new_incl) >= 4,
      TRUE ~ TRUE  # Keep all rows when hh_size does not match any of the above conditions
    )) %>%
    mutate(household_no_new_incl = dense_rank(household_no_new_incl))
  
  # Data of each household
  df_data_household <- data.frame(
    household_sizes_old_incl,
    household_sizes_new_incl
  ) %>%
    # specific household sizes
    filter(case_when(
      hh_size == "2" ~ household_sizes_new_incl == "2",
      hh_size == "3" ~ household_sizes_new_incl == "3",
      hh_size == "4" ~ household_sizes_new_incl == "4",
      hh_size == "5p" ~ as.integer(household_sizes_new_incl) >= 5,
      hh_size == "2or3" ~ household_sizes_new_incl %in% c("2", "3"),
      hh_size == "4p" ~ as.integer(household_sizes_new_incl) >= 4,
      TRUE ~ TRUE  # Keep all rows when hh_size does not match any of the above conditions
    ))
  
  # Save data
  save(df_data_individual,
       df_data_household,
       household_months_incl,
       file = paste0("../R/Results-temp/", scn, "/RData/data_initial.RData"))
}
