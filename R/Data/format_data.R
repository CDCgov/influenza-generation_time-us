## format_data.m ###############################################################
format_data <- function(i, param_seq, scn) {
  
  ## format_data.m ###############################################################
  # Format imported data into a structure array.
  
  # Throughout, arrays with the suffix "_dir" are ordered according to
  # (i) the household number assigned to each household, and
  # (ii) the (unknown) order in which household members became infected.
  # Otherwise, arrays are ordered according to
  # (i) household number, and
  # (ii) the (known) order in which household members developed symptoms.
  
  # Load imported data from "data_initial.RData"
  # (this file is created by running "import_data.R")
  load(paste0("../R/Results-temp/", scn, "/RData/data_initial.RData"))
  
  # Data of each household
  household_sizes_old_incl <- df_data_household$household_sizes_old_incl
  household_sizes_new_incl <- df_data_household$household_sizes_new_incl
  # Data of each individual
  household_no_new_incl         <- df_data_individual$household_no_new_incl
  household_size_indiv_old_incl <- df_data_individual$household_size_indiv_old_incl
  household_size_indiv_new_incl <- df_data_individual$household_size_indiv_new_incl
  household_months_incl         <- df_data_individual$household_months_incl
  d_s_incl                      <- df_data_individual$d_s_incl
  d_iR_incl                     <- df_data_individual$d_iR_incl
  symp_incl                     <- df_data_individual$symp_incl
  asymp_incl                    <- df_data_individual$asymp_incl
  uninf_incl                    <- df_data_individual$uninf_incl
  
  # Number of hosts and households
  no_hosts      <- length(d_s_incl)
  no_households <- length(household_sizes_new_incl)
  
  # Simplify names of data arrays
  household_sizes_incl      <- household_sizes_new_incl
  household_sizes_full      <- household_sizes_old_incl
  household_size_indiv_incl <- household_size_indiv_new_incl
  household_size_indiv_full <- household_size_indiv_old_incl
  household_no              <- household_no_new_incl
  household_months          <- household_months_incl
  
  d_s            <- d_s_incl
  d_iR           <- d_iR_incl
  uninfected_dir <- uninf_incl
  infected_dir   <- !uninfected_dir
  symp           <- symp_incl
  asymp          <- asymp_incl
  
  # Lower and upper bounds for symptom onset and infection times
  t_sL <- d_s - 0.5
  t_sR <- d_s + 0.5
  t_iL <- rep(-Inf, no_hosts)
  t_iR <- d_iR + 0.5
  
  # Indicator matrix showing which household each individual belongs to
  blocks <- lapply(household_sizes_incl, function(size) matrix(1, nrow = size, ncol = 1))
  household_indicator_mat <- as.matrix(Matrix::bdiag(blocks))
  
  # Calculate the number of infected and symptomatic/asymptomatic hosts in
  # each household, and working in the (unknown) order of infection,
  # determine the possible (household) infectors for each individual (i.e.,
  # household members who were infected before that individual).
  no_infected_in_household <- numeric(no_households)   # number of infected     in each household
  no_symp_in_household     <- numeric(no_households)   # number of symptomatic  in each household
  no_asymp_in_household    <- numeric(no_households)   # number of asymptomatic in each household
  no_poss_infectors_dir    <- numeric(no_hosts)        # number of possible infectors for each individual
  poss_infectors_dir_cell  <- vector('list', no_hosts) # cell array to populate with indices of possible infectors for each individual
  primary_dir              <- logical(no_hosts)        # logical array to indicate which individual is the primary case in the household
  
  for (i in 1:no_households) {
    in_household                  <- household_no == i
    infected_hosts_in_household   <- which(in_household & infected_dir)   # individual id of infected
    symp_hosts_in_household       <- which(in_household & symp)           # individual id of symptomatic
    asymp_hosts_in_household      <- which(in_household & asymp)          # individual id of asymptomatic
    uninfected_hosts_in_household <- which(in_household & uninfected_dir) # individual id of uninfected
    
    no_infected_in_household[i] <- length(infected_hosts_in_household)    # number of infected     in each household
    no_symp_in_household[i]     <- length(symp_hosts_in_household)        # number of symptomatic  in each household
    no_asymp_in_household[i]    <- length(asymp_hosts_in_household)       # number of asymptomatic in each household
    
    no_poss_infectors_dir[infected_hosts_in_household[1]]     <- 1    # number of possible infectors for each individual
    poss_infectors_dir_cell[[infected_hosts_in_household[1]]] <- 0    # possible infectors (use 0 to denote infection from outside the household)
    primary_dir[infected_hosts_in_household[1]]               <- TRUE # primary index case
    
    if (length(infected_hosts_in_household) >= 2) {
      for (j in 2:length(infected_hosts_in_household)) { # except index case
        poss_infectors_dir_cell[[infected_hosts_in_household[j]]] <- infected_hosts_in_household[1:(j-1)]
        no_poss_infectors_dir[infected_hosts_in_household[j]]     <- length(poss_infectors_dir_cell[[infected_hosts_in_household[j]]])
      }
    }
    
    if (length(uninfected_hosts_in_household) >= 1) {
      for (j in 1:length(uninfected_hosts_in_household)) {
        poss_infectors_dir_cell[[uninfected_hosts_in_household[j]]] <- infected_hosts_in_household
        no_poss_infectors_dir[uninfected_hosts_in_household[j]]     <- length(infected_hosts_in_household)
      }
    }
  }
  
  # Logical arrays indicating whether or not each household contained any symptomatic or asymptomatic infected hosts
  symp_in_household  <- (no_symp_in_household  > 0)
  asymp_in_household <- (no_asymp_in_household > 0)
  
  # Create a vector v, enumerating each value within the cell array listing the possible infectors for each individual (in the order of infection)
  v <- unlist(poss_infectors_dir_cell)
  
  # Indicator matrix giving the index of the potential infector corresponding to each entry of v
  M1 <- matrix(0, nrow = length(v), ncol = no_hosts)
  for (i in 1:length(v)) {
    j <- v[i]
    
    if (j > 0) {
      M1[i, j] <- 1
    }
  }
  
  # Indicator matrix giving the index of the potential infectee corresponding to each entry of v
  blocks <- lapply(no_poss_infectors_dir, function(size) matrix(1, nrow = size, ncol = 1))
  M2 <- as.matrix(Matrix::bdiag(blocks))
  
  # Logical vectors indicating, respectively, whether or not each entry of v corresponds to:
  # (i) the potential infector of a host who actually became infected;
  # (ii) the infection of the primary case (denoted by a 0 in v);
  # (iii) the potential infector of a non-primary case who became infected
  v_to_infected_indicator  <- as.logical(M2 %*% infected_dir)
  v_to_primary_indicator   <- (v == 0)
  v_to_recipient_indicator <- (v_to_infected_indicator & !v_to_primary_indicator)
  
  # The size of the household that each entry of v corresponds to
  household_size_v <- as.numeric(M2 %*% household_size_indiv_full)
  
  # Structure list containing information about possible infectors
  poss_infectors_dir <- list(
    cell                   = poss_infectors_dir_cell,
    all                    = v,
    from_indicator_mat     = M1,
    to_indicator_mat       = M2,
    to_primary_indicator   = v_to_primary_indicator,
    to_infected_indicator  = v_to_infected_indicator,
    to_recipient_indicator = v_to_recipient_indicator,
    household_size         = household_size_v
  )
  
  # Structure list containing all observed data
  data_struct_observed <- list(
    t_iL                         = t_iL,
    t_iR                         = t_iR,
    t_sL                         = t_sL,
    t_sR                         = t_sR,
    household_sizes_incl         = household_sizes_incl,
    household_sizes_full         = household_sizes_full,
    household_size_indiv_incl    = household_size_indiv_incl,
    household_size_indiv_full    = household_size_indiv_full,
    household_no                 = household_no,
    household_indicator_mat      = household_indicator_mat,
    no_infected_in_household     = no_infected_in_household,
    no_symp_in_household         = no_symp_in_household,
    symp_in_household            = symp_in_household,
    no_asymp_in_household        = no_asymp_in_household,
    asymp_in_household           = asymp_in_household,
    primary_dir                  = primary_dir,
    infected_dir                 = infected_dir,
    symp                         = symp,
    asymp                        = asymp,
    no_poss_infectors_dir        = no_poss_infectors_dir,
    poss_infectors_dir           = poss_infectors_dir,
    household_months             = household_months
  )
  
  # Save data to .RData file
  save(data_struct_observed, file = paste0("../R/Results-temp/", scn, "/RData/data.RData"))
}
