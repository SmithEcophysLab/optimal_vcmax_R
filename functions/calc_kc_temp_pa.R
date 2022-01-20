# modify of nick's function calc_km_pa in the optimal-vcmax repository

calc_kc_temp_pa <- function(temp, z) {
  temp_k <- 273.15 + temp
  
  patm <- calc_patm(z) 
  rat <- patm / calc_patm(0)
  
  R <- 8.314   # in joules    
  
  Kc25 <- 121 # Pa CO2 Boyd 2015 
  Kc25 <- Kc25 * rat
  Hkc <-  64200 # J mol^-1
  
  Kc_pa <- Kc25 * exp(Hkc * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  
  return(Kc_pa)
}