calc_kp_temp_pa <- function(temp, z){ # Boyd 2015
  patm = calc_patm(z) 
  rat = patm / calc_patm(0)
  
  R <- 8.31 # in Joules
  temp_k <- temp + 273
  kp_25 <- 13.9 # Pa Co2
  kp_25 <- kp_25 * rat
  Ea <- 36.3 # kJ mol-1
  Ea_j <- Ea * 1000
  kp <- kp_25 * exp(Ea_j * (temp_k - 298.15)/(298.25 * R * temp_k))
  return(kp)
}