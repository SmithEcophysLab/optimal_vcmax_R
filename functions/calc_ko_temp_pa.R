# modify of nick's function calc_km_pa in the optimal-vcmax repository

calc_ko_temp_pa <- function(temp, z) {
  temp_k <- 273.15 + temp
  
  patm <- calc_patm(z) 
  rat <- patm / calc_patm(0)
  
  R <- 8.314   # in joules    
  
  Ko25_kpa <- 29.2 # Boyd 2015
  Ko25 <- Ko25_kpa * 1000
  Ko25 <- Ko25 * rat
  Hko <-  10500
 
  Ko_pa <- Ko25 * exp(Hko * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  
  return(Ko_pa)
}