calc_jmax_temp <- function(acclim_jmax, tg_c, R){
  # tg_k <- tg_c + 273 # Kelvin
  # t_ref <- 25 + 273 # Kelvin
  # t_acclim <- 25 + 273 # Kelvin
  # 
  # Ha <- 49.884
  # Hd <- 200
  
  # a <- 659.70
  # b <- -0.75
  
  # activation_term <- Ha * (tg_k - t_ref) / (t_ref * R * tg_k)
  # deactivation_term_num <- 1 + exp((t_ref * (a + b * t_acclim) - Hd) /
  #                                    (t_ref * R))
  # deactivation_term_denom <- 1 + exp((tg_k * (a + b * t_acclim) - Hd) /
  #                                      (tg_k * R))
  # 
  # jmax_temp <- acclim_jmax * exp(activation_term) * 
  #   (deactivation_term_num / deactivation_term_denom)
  tref <- 25
  tmean <- 25
  temp <- tg_c + 273.15
  Ha <- 49884 # Activation Energy (J mol^-1) Massad 2007
  Hd <- 200000 # Deactivaiton energy (J mol^-1) Massad 2007
  adelS <- 659.70
  bdelS <- -0.75
  
  trefK <- tref + 273.15
  R <- 8.314 # Universal Gas constant (J mol^-1 K^-1)
  
  kbeg <- exp(Ha*(temp-trefK)/(trefK*R*temp))
  kend <-((1+exp((trefK*(adelS+bdelS*tmean)-Hd)/(trefK*R)))/(1+exp((temp*(adelS+bdelS*tmean)-Hd)/(temp*R))))
  
  ktotal <- kbeg*kend # Equation 20 in Smith 2019
  jmax_temp <- acclim_jmax * ktotal
  return(jmax_temp)
}