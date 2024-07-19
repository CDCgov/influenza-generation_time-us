## f_tost_form_mech.m ##########################################################
f_tost_form_mech <- function(params, t) {
  # TOST distribution for our mechanistic approach with parameters given
  # by params.
  
  gamma <- params[1]
  mu    <- params[2]
  k_inc <- params[3]
  k_E   <- params[4]
  k_I   <- params[5]
  alpha <- params[6]
  k_P   <- k_inc - k_E
  
  # parameter (page 23)
  C <- k_inc * gamma * mu / (alpha * k_P * mu + k_inc * gamma)
  
  # TOST distribution (page 27)
  fm <- alpha * C * (1 - pgamma(-t, shape = k_P, scale = 1 / (k_inc * gamma)))
  fp <-         C * (1 - pgamma( t, shape = k_I, scale = 1 / (k_I   * mu   )))
  
  indm <- (t <  0)
  ind0 <- (t == 0)
  indp <- !(indm | ind0)
  
  f_tost_mech       <- rep(NA, length(t))
  f_tost_mech[indm] <- fm[indm]
  f_tost_mech[indp] <- fp[indp]
  f_tost_mech[ind0] <- fp[ind0]
  
  return(f_tost_mech)
}
