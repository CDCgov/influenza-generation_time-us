## update_onset_fun.m ##########################################################
update_onset_fun <- function(theta, data_struct_augmented_old, ll_household_old, ll_household_form) {
  # Update the symptom onset times of symptomatic infected hosts in the
  # parameter fitting procedure for either model
  
  # Throughout, arrays with the suffix "_dir" are ordered according to
  # (i) the household number assigned to each household, and (ii) the
  # (unknown) order in which household members became infected.
  # Otherwise, arrays are ordered according to (i) household number, and
  # (ii) the (known) order in which household members developed symptoms
  
  # Left and right bounds for times of infection
  t_sL <- data_struct_augmented_old$t_sL
  t_sR <- data_struct_augmented_old$t_sR
  
  # Old augmented data
  t_s_old         <- data_struct_augmented_old$t_s
  t_s_dir_old     <- data_struct_augmented_old$t_s_dir
  t_dir_host_inds <- data_struct_augmented_old$t_dir_host_inds
  symp            <- data_struct_augmented_old$symp
  
  # Household details
  household_sizes_incl    <- data_struct_augmented_old$household_sizes_incl
  household_indicator_mat <- data_struct_augmented_old$household_indicator_mat
  symp_in_household       <- data_struct_augmented_old$symp_in_household
  no_households           <- length(household_sizes_incl)
  
  # Resample the symptom onset times of all symptomatic infected hosts
  t_s_prop       <- t_s_old
  r              <- runif(sum(symp))
  t_s_prop[symp] <- t_sL[symp] + (t_sR[symp] - t_sL[symp]) * r
  
  # Proposed augmented data ordered according to times of infection within each household
  t_s_dir_prop <- t_s_prop[t_dir_host_inds]
  
  # Populate the structure array data_struct_augmented_prop with the proposed onset times
  data_struct_augmented_prop         <- data_struct_augmented_old
  data_struct_augmented_prop$t_s     <- t_s_prop
  data_struct_augmented_prop$t_s_dir <- t_s_dir_prop
  
  # Calculate the log-likelihood contribution from each household using the proposed onset times
  ll_household_prop <- ll_household_form(theta, data_struct_augmented_prop)
  
  # Calculate the logarithm of the ratio between the proposed and previous likelihood contributions from each household
  la_vec                 <- ll_household_prop - ll_household_old
  la_vec[is.nan(la_vec)] <- -Inf
  
  # Accept the changes for each household with probabilities given by the entries in exp(la_vec)
  accept_households <- (log(runif(no_households)) < la_vec)
  accept_hosts      <- as.logical(household_indicator_mat %*% accept_households)
  
  # New augmented data after accepting/rejecting proposals
  t_s_new               <- t_s_old
  t_s_new[accept_hosts] <- t_s_prop[accept_hosts]
  
  t_s_dir_new               <- t_s_dir_old
  t_s_dir_new[accept_hosts] <- t_s_dir_prop[accept_hosts]
  
  # Populate the structure array data_struct_augmented_new with the new infection times
  data_struct_augmented_new         <- data_struct_augmented_old
  data_struct_augmented_new$t_s     <- t_s_new
  data_struct_augmented_new$t_s_dir <- t_s_dir_new
  
  # New likelihood contributions from each household
  ll_household_new                    <- ll_household_old
  ll_household_new[accept_households] <- ll_household_prop[accept_households]
  
  # Information about acceptance (by household and overall)
  acceptance <- list(
    household = accept_households,
    overall   = mean(accept_households[symp_in_household])
  )
  
  return(list(data_struct_augmented = data_struct_augmented_new,
              ll_household          = ll_household_new,
              acceptance            = acceptance))
}
