## initialise_augmented_data_mech.m ############################################
initialise_augmented_data_mech <- function(data_struct_observed) {
  # Initialise the structure array of augmented data,
  # data_struct_augmented, from the structure array of observed data,
  # data_struct_observed, for the mechanistic model
  
  # Throughout, arrays with the suffix "_dir" are ordered according to
  # (i) the household number assigned to each household, and (ii) the
  # (unknown) order in which household members became infected.
  # Otherwise, arrays are ordered according to (i) household number, and
  # (ii) the (known) order in which household members developed symptoms.
  
  # Household number that each individual belongs to, and the total
  # number of individuals
  household_no <- data_struct_observed$household_no
  no_hosts     <- length(household_no)
  
  # Left and right bounds for infection and symptom onset times
  t_iL <- data_struct_observed$t_iL
  t_iR <- data_struct_observed$t_iR
  t_sL <- data_struct_observed$t_sL
  t_sR <- data_struct_observed$t_sR
  
  # Logical arrays indicating symptom status
  symp  <- data_struct_observed$symp
  asymp <- data_struct_observed$asymp
  
  # Initialise vectors t_i and t_s containing infection and symptom onset
  # times for each individual
  t_s <- rep(Inf, no_hosts)
  t_i <- rep(Inf, no_hosts)
  
  t_s[symp] <- t_sL[symp] + (t_sR[symp] - t_sL[symp]) * runif(sum(symp))
  t_i[symp] <- t_s[symp] - 6
  
  if (sum(symp) > 0) {
    t_s_min <- min(t_s[symp])
    t_s_max <- max(t_s[symp])
    t_i[asymp] <- t_s_min + (t_s_max - t_s_min) * runif(sum(asymp))
  } else {
    t_i[asymp] <- runif(sum(asymp))
  }
  
  t_i <- pmax(pmin(t_i, t_iR), t_iL)  # ensure initial values lie within bounds
  
  # Assign values to "onset" times of asymptomatic hosts (i.e., time of
  # entry into the I stage)
  t_s[asymp] <- t_i[asymp] + 6
  
  # Augmented data ordered by (estimated) time of infection within the
  # household
  t_dir_host_inds <- order(household_no, t_i)
  t_i_dir         <- t_i[t_dir_host_inds]
  t_s_dir         <- t_s[t_dir_host_inds]
  
  symp_dir  <- symp[t_dir_host_inds]
  asymp_dir <- asymp[t_dir_host_inds]
  
  # Store augmented data in a structure array
  data_struct_augmented                 <- data_struct_observed
  data_struct_augmented$t_i             <- t_i
  data_struct_augmented$t_s             <- t_s
  data_struct_augmented$t_i_dir         <- t_i_dir
  data_struct_augmented$t_s_dir         <- t_s_dir
  data_struct_augmented$t_dir_host_inds <- t_dir_host_inds
  data_struct_augmented$symp_dir        <- symp_dir
  data_struct_augmented$asymp_dir       <- asymp_dir
  
  return(data_struct_augmented)
}
