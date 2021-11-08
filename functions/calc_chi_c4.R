# calc_chi_c4


calc_chi_c4 <- function(ca, temp, vpd, z){
  patm <- calc_patm(z)
  ca_pa <- ca * 1e-6 * patm
  
  R <- 8.3145 # Pa
  tempK <- temp + 273.15
  vpd_pa <- vpd * 1000
  beta <- 166

  Kp_25 <- 16 # Pa
  Ea_Kp <- 36300 # J mol^-1, for the Pa parameter ## Boyd et al 2015
  Kp <- Kp_25 * exp((Ea_Kp * (tempK - 298.15))/(298.15 * R * tempK))
  
  eta_star <- calc_nstar(temp, z)
  
  xi <- sqrt((beta * Kp) / (1.6 * eta_star))
  chi <- (xi / (xi + sqrt(vpd_pa)))
  # res <- c(chi, ci)
  res <- chi
  
  return(res)
}
