## empirical_summary_mech.m ####################################################
empirical_summary_mech <- function(b_cond, data_struct_augmented) {
  # Estimate the mean and standard deviation of realised household
  # generation times, in addition to the proportion of presymptomatic
  # transmissions, for the mechanistic model with augmented data
  # data_struct_augmented.
  
  # Throughout, arrays with the suffix "_dir" are ordered according to
  # (i) the household number assigned to each household, and (ii) the
  # (unknown) order in which household members became infected.
  
  # Augmented infection and symptom onset times
  t_i_dir <- data_struct_augmented$t_i_dir
  t_s_dir <- data_struct_augmented$t_s_dir
  
  # Logical vectors indicating symptom status
  symp_dir  <- data_struct_augmented$symp_dir
  asymp_dir <- data_struct_augmented$asymp_dir
  
  # Information about possible infectors
  poss_infectors_dir       <- data_struct_augmented$poss_infectors_dir
  v                        <- poss_infectors_dir$all # list of all possible infectors for each host in order
  M_from                   <- poss_infectors_dir$from_indicator_mat
  M_to                     <- poss_infectors_dir$to_indicator_mat
  to_recipient_indicator_v <- poss_infectors_dir$to_recipient_indicator
  household_size_v         <- poss_infectors_dir$household_size
  
  from_symp_indicator_v    <- as.logical(M_from %*% symp_dir)
  from_asymp_indicator_v   <- as.logical(M_from %*% asymp_dir)
  
  # Vectors of every possible generation time for each infectee (with
  # entries corresponding to infection by every possible infector), and
  # of the incubation period (which is infinite for asymptomatic hosts)
  # of the corresponding infectee.
  t_i_dir[t_i_dir == Inf]                      <- 999999999
  t_s_dir[t_s_dir == Inf]                      <- 999999999
  t_gen_contribs                               <- (M_to - M_from) %*% t_i_dir # vector of every possible generation time for each infected host
  t_tost_contribs                              <- M_to %*% t_i_dir - M_from %*% t_s_dir
  t_inc_contribs                               <- M_from %*% (t_s_dir - t_i_dir)
  t_gen_contribs[t_gen_contribs > 900000000]   <- Inf
  t_tost_contribs[t_tost_contribs > 900000000] <- Inf
  t_inc_contribs[t_inc_contribs > 900000000]   <- Inf
  
  # Calculate likelihood contributions corresponding to infection of each individual
  t_gen_recipient_contribs          <- t_gen_contribs[to_recipient_indicator_v]
  t_tost_recipient_contribs         <- t_tost_contribs[to_recipient_indicator_v]
  t_inc_recipient_contribs          <- t_inc_contribs[to_recipient_indicator_v]
  household_size_recipient_contribs <- household_size_v[to_recipient_indicator_v]
  from_asymp_recipient_contribs     <- from_asymp_indicator_v[to_recipient_indicator_v]
  
  L2a_contribs                           <- rep(0, length(v))
  L2a_contribs[to_recipient_indicator_v] <- b_cond(t_tost_recipient_contribs, t_inc_recipient_contribs, household_size_recipient_contribs, from_asymp_recipient_contribs)
  L2a                                    <- t(M_to) %*% L2a_contribs
  
  # Non-trivial possible generation times and times from symptom onset to
  # transmission (TOST)
  t_gen  <- t_gen_recipient_contribs
  t_tost <- t_tost_contribs[from_symp_indicator_v & to_recipient_indicator_v]
  
  # Use likelihood contributions to calculate weights indicating the
  # relative probabilities of infection by different possible infectors
  weights_all                      <- L2a_contribs / (M_to %*% L2a)
  weights_all[is.nan(weights_all)] <- 0
  
  weights_gen  <- weights_all[to_recipient_indicator_v]
  weights_tost <- weights_all[from_symp_indicator_v & to_recipient_indicator_v]
  
  weights_gen  <- weights_gen  / sum(weights_gen)
  weights_tost <- weights_tost / sum(weights_tost)
  
  # Calculate the mean and standard deviation of realized generation
  # times, using the computed weights to account for different
  # probabilities that individuals were infected by different infectors.
  m <- sum(t_gen * weights_gen)
  s <- sqrt(sum(t_gen * t_gen * weights_gen) - m^2)
  p <- sum((t_tost < 0) * weights_tost)
  
  empirical_summary <- c(m, s, p)
  return(empirical_summary)
}
