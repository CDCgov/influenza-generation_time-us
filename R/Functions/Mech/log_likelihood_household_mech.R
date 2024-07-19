## log_likelihood_household_mech.m ##########################################################
log_likelihood_household_mech <- function(f_inc, b_cond, B_cond, mean_transmissions, data_struct_augmented) {
  # Calculate the contribution from each household to the
  # log-likelihood of the augmented data, data_struct_augmented, for the
  # mechanistic model.
  
  # Throughout, arrays with the suffix "_dir" are ordered according to
  # (i) the household number assigned to each household, and (ii) the
  # (unknown) order in which household members became infected.
  
  # Augmented infection and symptom onset times
  t_i_dir <- data_struct_augmented$t_i_dir
  t_s_dir <- data_struct_augmented$t_s_dir
  
  # Number of hosts
  no_hosts <- length(t_i_dir)
  
  # Logical vectors indicating infection/symptom/primary status
  primary_dir    <- data_struct_augmented$primary_dir
  infected_dir   <- data_struct_augmented$infected_dir
  uninfected_dir <- !infected_dir
  asymp_dir      <- data_struct_augmented$asymp_dir
  
  # Indicator matrix showing which household each individual belongs to
  household_indicator_mat <- data_struct_augmented$household_indicator_mat
  
  # Information about possible infectors
  poss_infectors_dir        <- data_struct_augmented$poss_infectors_dir
  v                         <- poss_infectors_dir$all # list of all possible infectors for each host in order
  M_from                    <- poss_infectors_dir$from_indicator_mat
  M_to                      <- poss_infectors_dir$to_indicator_mat
  to_uninfected_indicator_v <- !poss_infectors_dir$to_infected_indicator
  to_recipient_indicator_v  <- poss_infectors_dir$to_recipient_indicator
  household_size_v          <- poss_infectors_dir$household_size
  
  from_asymp_indicator_v    <- as.logical(M_from %*% asymp_dir)
  
  # Contribution to likelihood from incubation periods
  t_i_dir_inf <- t_i_dir[infected_dir]
  t_s_dir_inf <- t_s_dir[infected_dir]
  
  t_inc_inf <- t_s_dir_inf - t_i_dir_inf
  l1_inf    <- log(f_inc(t_inc_inf))
  
  l1_indiv               <- numeric(no_hosts)
  l1_indiv[infected_dir] <- l1_inf
  
  # Contribution to likelihood from transmissions occurring
  t_i_dir[t_i_dir == Inf]                      <- 999999999
  t_s_dir[t_s_dir == Inf]                      <- 999999999
  t_tost_contribs                              <- M_to %*% t_i_dir - M_from %*% t_s_dir # vector of every possible generation time for each infected host
  t_inc_contribs                               <- M_from %*% (t_s_dir - t_i_dir)
  t_tost_contribs[t_tost_contribs > 900000000] <- Inf
  t_inc_contribs[t_inc_contribs > 900000000]   <- Inf
  
  t_tost_recipient_contribs         <- t_tost_contribs[to_recipient_indicator_v]
  t_inc_recipient_contribs          <- t_inc_contribs[to_recipient_indicator_v]
  household_size_recipient_contribs <- household_size_v[to_recipient_indicator_v]
  from_asymp_recipient_contribs     <- from_asymp_indicator_v[to_recipient_indicator_v]
  
  t_inc_uninfected_contribs          <- t_inc_contribs[to_uninfected_indicator_v]
  household_size_uninfected_contribs <- household_size_v[to_uninfected_indicator_v]
  from_asymp_uninfected_contribs     <- from_asymp_indicator_v[to_uninfected_indicator_v]
  
  L2a_contribs                           <- numeric(length(v))
  L2a_contribs[to_recipient_indicator_v] <- b_cond(t_tost_recipient_contribs, t_inc_recipient_contribs, household_size_recipient_contribs, from_asymp_recipient_contribs)
  
  L2a                               <- t(M_to) %*% L2a_contribs
  L2a[primary_dir | uninfected_dir] <- 1
  
  l2a                               <- log(L2a)
  
  # Contribution to likelihood from evasion of infection up to time of infection (or for all time for individuals who remained uninfected)
  l2b_contribs                            <- numeric(length(v))
  l2b_contribs[to_recipient_indicator_v]  <- B_cond(t_tost_recipient_contribs, t_inc_recipient_contribs, household_size_recipient_contribs, from_asymp_recipient_contribs)
  l2b_contribs[to_uninfected_indicator_v] <- mean_transmissions(t_inc_uninfected_contribs, household_size_uninfected_contribs, from_asymp_uninfected_contribs)
  
  l2b                                     <- -t(M_to) %*% l2b_contribs
  
  # Calculate overall likelihood contribution from each individual
  l2_indiv <- l2a + l2b
  l_indiv  <- l1_indiv + l2_indiv
  
  # Likelihood contribution from each household
  l_indiv[l_indiv == -Inf]              <- -999999999
  l_household                           <- t(household_indicator_mat) %*% l_indiv
  l_household[l_household < -900000000] <- -Inf
  
  return(l_household)
}
