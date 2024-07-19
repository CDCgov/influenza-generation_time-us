## update_asymp_fun_indep.m ####################################################
update_asymp_fun_indep <- function(theta, data_struct_augmented_old, ll_household_old, ll_household_form, t_i_prop_sd_asymp) {
  # Update the infection times of asymptomatic infected hosts in the
  # parameter fitting procedure for the independent transmission and
  # symptoms model
  
  # Throughout, arrays with the suffix "_dir" are ordered according to
  # (i) the household number assigned to each household, and (ii) the
  # (unknown) order in which household members became infected.
  # Otherwise, arrays are ordered according to (i) household number, and
  # (ii) the (known) order in which household members developed symptoms
  
  # Left and right bounds for times of infection
  t_iL <- data_struct_augmented_old$t_iL
  t_iR <- data_struct_augmented_old$t_iR
  
  # Old augmented data
  t_i_old             <- data_struct_augmented_old$t_i
  t_s                 <- data_struct_augmented_old$t_s
  t_i_dir_old         <- data_struct_augmented_old$t_i_dir
  t_s_dir_old         <- data_struct_augmented_old$t_s_dir
  t_dir_host_inds_old <- data_struct_augmented_old$t_dir_host_inds
  symp                <- data_struct_augmented_old$symp
  asymp               <- data_struct_augmented_old$asymp
  symp_dir_old        <- data_struct_augmented_old$symp_dir
  asymp_dir_old       <- data_struct_augmented_old$asymp_dir
  
  # Household details
  household_no             <- data_struct_augmented_old$household_no
  household_sizes_incl     <- data_struct_augmented_old$household_sizes_incl
  household_indicator_mat  <- data_struct_augmented_old$household_indicator_mat
  no_infected_in_household <- data_struct_augmented_old$no_infected_in_household
  no_symp_in_household     <- data_struct_augmented_old$no_symp_in_household
  no_asymp_in_household    <- data_struct_augmented_old$no_asymp_in_household
  asymp_in_household       <- data_struct_augmented_old$asymp_in_household
  no_households            <- length(household_sizes_incl)
  
  # Randomly select one asymptomatic infected host in every household
  # containing at least one such individual
  r1                            <- runif(sum(asymp_in_household))
  update_hosts_in_household     <- ceiling(no_asymp_in_household[asymp_in_household] * r1)
  household_asymp_start_indices <- cumsum(household_sizes_incl) - household_sizes_incl + no_symp_in_household
  update_hosts                  <- household_asymp_start_indices[asymp_in_household] + update_hosts_in_household
  
  # Update the infection times of the chosen hosts
  r2                     <- rnorm(sum(asymp_in_household))
  t_i_updates            <- t_i_prop_sd_asymp * r2
  t_i_prop               <- t_i_old
  t_i_prop[update_hosts] <- t_i_prop[update_hosts] + t_i_updates # must work with undirected vectors so asymptomatics placed after symptomatics
  
  # Determine occasions when the proposed infection times lie outside the allowed bounds
  t_i_bdry_hosts_prop      <- apply(cbind(t_i_prop < t_iL, t_i_prop > t_iR), 1, any)
  t_i_bdry_households_prop <- as.logical(t(household_indicator_mat) %*% t_i_bdry_hosts_prop)
  
  # Proposed augmented data ordered according to times of infection within each household
  sort_indices   <- order(household_no, t_i_prop)
  t_i_dir_prop   <- t_i_prop[sort_indices]
  t_s_dir_prop   <- t_s[sort_indices]
  symp_dir_prop  <- symp[sort_indices]
  asymp_dir_prop <- asymp[sort_indices]
  
  # Populate the structure array data_struct_augmented_prop with the proposed infection times
  data_struct_augmented_prop           <- data_struct_augmented_old
  data_struct_augmented_prop$t_i       <- t_i_prop
  data_struct_augmented_prop$t_i_dir   <- t_i_dir_prop
  data_struct_augmented_prop$t_s_dir   <- t_s_dir_prop
  data_struct_augmented_prop$symp_dir  <- symp_dir_prop
  data_struct_augmented_prop$asymp_dir <- asymp_dir_prop
  
  # Calculate the log-likelihood contribution from each household using the proposed parameters and onset times
  ll_household_prop <- ll_household_form(theta, data_struct_augmented_prop)
  
  # If proposed times for a household lie outside the allowed bounds, set the likelihood contribution to 0
  ll_household_prop[t_i_bdry_households_prop] <- -Inf
  
  # Calculate the logarithm of the ratio between the proposed and previous likelihood contributions from each household
  la_vec                 <- ll_household_prop - ll_household_old
  la_vec[is.nan(la_vec)] <- -Inf
  
  # Accept the changes for each household with probabilities given by the entries in exp(la_vec)
  accept_households <- (log(runif(no_households)) < la_vec)
  accept_hosts      <- as.logical(household_indicator_mat %*% accept_households)
  
  # New augmented data after accepting/rejecting proposals
  t_i_new               <- t_i_old
  t_i_new[accept_hosts] <- t_i_prop[accept_hosts]
  
  t_i_dir_new               <- t_i_dir_old
  t_i_dir_new[accept_hosts] <- t_i_dir_prop[accept_hosts]
  
  t_s_dir_new               <- t_s_dir_old
  t_s_dir_new[accept_hosts] <- t_s_dir_prop[accept_hosts]
  
  t_dir_host_inds_new               <- t_dir_host_inds_old
  t_dir_host_inds_new[accept_hosts] <- sort_indices[accept_hosts]
  
  symp_dir_new               <- symp_dir_old
  symp_dir_new[accept_hosts] <- symp_dir_prop[accept_hosts]
  
  asymp_dir_new               <- asymp_dir_old
  asymp_dir_new[accept_hosts] <- asymp_dir_prop[accept_hosts]
  
  # Populate the structure array data_struct_augmented_new with the new infection times
  data_struct_augmented_new                 <- data_struct_augmented_old
  data_struct_augmented_new$t_i             <- t_i_new
  data_struct_augmented_new$t_i_dir         <- t_i_dir_new
  data_struct_augmented_new$t_s_dir         <- t_s_dir_new
  data_struct_augmented_new$t_dir_host_inds <- t_dir_host_inds_new
  data_struct_augmented_new$symp_dir        <- symp_dir_new
  data_struct_augmented_new$asymp_dir       <- asymp_dir_new
  
  # New likelihood contributions from each household
  ll_household_new                    <- ll_household_old
  ll_household_new[accept_households] <- ll_household_prop[accept_households]
  
  # Information about acceptance (by household and overall)
  acceptance <- list(
    household = accept_households,
    overall   = mean(accept_households[asymp_in_household & (no_infected_in_household > 1)]) # Only include households with at least one asymptomatic and one other infection so this step is non-trivial
  )
  acceptance$infection <- acceptance$overall
  
  return(list(data_struct_augmented = data_struct_augmented_new,
              ll_household          = ll_household_new,
              acceptance            = acceptance))
}
