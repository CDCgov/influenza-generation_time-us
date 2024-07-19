## scenario_name.m #############################################################
scenario_name <- function(i, param_seq) {
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
  
  # scenario text
  virus_text     <- virus
  period_text    <- as.character(period)
  coprimary_text <- as.character(coprimary)
  hh_size_text   <- paste0("hh_", as.character(hh_size))
  type_text      <- as.character(type)
  no_steps_text <- format(no_steps, scientific = FALSE)
  inc_mean_text <- stringr::str_replace(sprintf("%.2f", inc_mean), "[.]", "p")
  inc_sd_text   <- stringr::str_replace(sprintf("%.2f", inc_sd),   "[.]", "p")
  x_A_text      <- stringr::str_replace(sprintf("%.2f", x_A),      "[.]", "p")
  
  # scenario name
  scn <- paste(index,
               virus_text, period_text, coprimary_text, hh_size_text, type_text,
               no_steps_text, inc_mean_text, inc_sd_text, x_A_text, sep = "-")
  
  return(scn)
}
