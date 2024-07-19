## log_likelihood_household_indep.m ############################################
log_likelihood_household_indep <- function(f_inc, beta0, rho, x_A, f_gen, F_gen, data_struct_augmented) {
  # Calculate the contribution from each household to the
  # log-likelihood of the augmented data, data_struct_augmented, for the
  # independent transmission and symptoms model.
  
  # Throughout, arrays with the suffix "_dir" are ordered according to
  # (i) the household number assigned to each household, and (ii) the
  # (unknown) order in which household members became infected.
  
  # Augmented infection and symptom onset times
  t_i_dir <- data_struct_augmented$t_i_dir
  t_s_dir <- data_struct_augmented$t_s_dir
  
  # Number of hosts
  no_hosts <- length(t_i_dir)
  
  # Logical vectors indicating infection/symptom/primary status
  infected_dir   <- data_struct_augmented$infected_dir
  uninfected_dir <- !infected_dir
  symp_dir       <- data_struct_augmented$symp_dir
  asymp_dir      <- data_struct_augmented$asymp_dir
  primary_dir    <- data_struct_augmented$primary_dir
  
  # Indicator matrix showing which household each individual belongs to
  household_indicator_mat <- data_struct_augmented$household_indicator_mat
  
  # Information about possible infectors
  poss_infectors_dir        <- data_struct_augmented$poss_infectors_dir
  v                         <- poss_infectors_dir$all # list of all possible infectors for each host in order
  M_from                    <- poss_infectors_dir$from_indicator_mat
  M_to                      <- poss_infectors_dir$to_indicator_mat
  to_uninfected_indicator_v <- !poss_infectors_dir$to_infected_indicator
  to_primary_indicator_v    <- poss_infectors_dir$to_primary_indicator
  to_recipient_indicator_v  <- poss_infectors_dir$to_recipient_indicator
  household_size_v          <- poss_infectors_dir$household_size
  
  from_asymp_indicator_v    <- as.logical(M_from %*% asymp_dir)
  
  beta_v                         <- beta0 / (household_size_v^rho)
  beta_v[to_primary_indicator_v] <- 0
  beta_v[from_asymp_indicator_v] <- x_A * beta_v[from_asymp_indicator_v]
  
  # Contribution to likelihood from incubation periods
  t_i_dir_symp <- t_i_dir[symp_dir]
  t_s_dir_symp <- t_s_dir[symp_dir]
  t_inc_symp   <- t_s_dir_symp - t_i_dir_symp
  l1_symp      <- log(f_inc(t_inc_symp))
  
  l1_indiv           <- numeric(no_hosts)
  l1_indiv[symp_dir] <- l1_symp
  
  # Contribution to likelihood from transmissions occurring
  t_i_dir[t_i_dir == Inf]                    <- 999999999
  t_gen_contribs                             <- (M_to - M_from) %*% t_i_dir # vector of every possible generation time for each infected host
  t_gen_contribs[t_gen_contribs > 900000000] <- Inf
  
  t_gen_recipient_contribs <- t_gen_contribs[to_recipient_indicator_v]
  beta_recipient_contribs  <- beta_v[to_recipient_indicator_v]
  
  beta_uninfected_contribs <- beta_v[to_uninfected_indicator_v]
  
  L2a_contribs                            <- numeric(length(v))
  L2a_contribs[to_recipient_indicator_v]  <- beta_recipient_contribs * f_gen(t_gen_recipient_contribs)
  
  L2a                                     <- t(M_to) %*% L2a_contribs
  L2a[primary_dir | uninfected_dir]       <- 1
  
  l2a                                     <- log(L2a)
  
  # Contribution to likelihood from evasion of infection up to time of
  # infection (or for all time for individuals who remained uninfected)
  l2b_contribs                            <- numeric(length(v))
  l2b_contribs[to_recipient_indicator_v]  <- beta_recipient_contribs * F_gen(t_gen_recipient_contribs)
  l2b_contribs[to_uninfected_indicator_v] <- beta_uninfected_contribs
  
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
