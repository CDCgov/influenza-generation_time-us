## get_gen_mean_sd_mech.m ######################################################
get_gen_mean_sd_mech <- function(params_mat) {
  # Calculate the mean and standard deviation of the generation time
  # distribution for our mechanistic approach, at each set of parameters
  # given by the params rows of params_mat.
  
  gamma <- params_mat[,1]
  mu    <- params_mat[,2]
  k_inc <- params_mat[,3]
  k_E   <- params_mat[,4]
  k_I   <- params_mat[,5]
  alpha <- params_mat[,6]
  k_P   <- k_inc - k_E
  
  C     <- k_inc * gamma * mu / (alpha * k_P * mu + k_inc * gamma)
  
  m_E <- k_E / (k_inc * gamma)
  v_E <- k_E / (k_inc * gamma)^2
  
  m_P   <- k_P / (k_inc * gamma)
  m_PP  <- k_P * (k_P + 1) / (k_inc * gamma)^2
  m_PPP <- k_P * (k_P + 1) * (k_P + 2) / (k_inc * gamma)^3
  
  m_I   <- 1 / mu
  m_II  <- (k_I + 1) / (k_I * mu^2)
  m_III <- (k_I + 1) * (k_I + 2) / (k_I^2 * mu^3)
  
  m_PI  <- m_P * m_I
  m_PPI <- m_PP * m_I
  m_PII <- m_P * m_II
  
  m_star <- (C / 2) * (alpha * m_PP + 2 * m_PI + m_II)
  v_star <- (C / 3) * (alpha * m_PPP + 3 * m_PPI + 3 * m_PII + m_III) - m_star^2
  
  m_gen <- m_E + m_star
  v_gen <- v_E + v_star
  s_gen <- sqrt(v_gen)
  
  return(list(m_gen = m_gen,
              s_gen = s_gen))
}
