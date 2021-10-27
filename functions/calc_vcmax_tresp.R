### function to scale vcmax to different temperatures
calc_vcmax_tresp_mult = function(tleaf, tmean, tref){
  
  temp <- tleaf + 273.15
  Ha <- 71513
  Hd <- 200000
  adelS <- 668.39
  bdelS <- -1.07
  tmeanK <- tmean + 273.15
  trefK <- tref + 273.15
  R <- 8.314
  kbeg <- exp(Ha*(temp - trefK) / (trefK * R * temp))
  kend <- ((1 + exp((trefK * (adelS + bdelS * tmean) - Hd) / 
                      (trefK * R))) / (1 + exp((temp * (adelS + bdelS * tmean) - Hd) / (temp * R))))
  kbeg * kend
  
}

### function to scale jmax to different temperatures
calc_jmax_tresp_mult = function(tleaf, tmean, tref){
  
  temp <- tleaf + 273.15
  Ha <- 49884
  Hd <- 200000
  adelS <- 659.7
  bdelS <- -0.75
  tmeanK <- tmean + 273.15
  trefK <- tref + 273.15
  R <- 8.314
  kbeg <- exp(Ha * (temp - trefK) / (trefK * R * temp))
  kend <- ((1 + exp((trefK * (adelS + bdelS * tmean) - Hd) / (trefK * R))) / 
             (1 + exp((temp * (adelS + bdelS * tmean) - Hd) / (temp * R))))
  kbeg * kend
  
}
