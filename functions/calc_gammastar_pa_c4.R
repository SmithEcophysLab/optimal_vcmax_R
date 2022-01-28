# calculate gammastar (Pa)

calc_gammastar_pa_c4 = function(temp, z) {
  
  patm <- calc_patm(z)
  rat <- calc_patm(z) / calc_patm(0)
  temp_k = 273.15 + temp
  
  gammastar25 <- 2.6 * rat  # Pa
  Hgm <- 37830 # J mol-1
  R <- 8.314        # J K-1 mol-1

  gStar_pa <- gammastar25*exp((Hgm/R)*(1/298.15-1/temp_k))
  
  return(gStar_pa)
}
