# function to calculature the temperature response multiplier for phi give a temperature in degrees Celsius

phi_ftemp = function(temperature){
  
  temperatureK = temperature + 273.15
  Ha = 4.178162e+04
  delS = 3.570332e+02
  Hd = 1.096021e+05
  R = 8.314
  trefK = 298.15
  
  kbeg = exp(Ha*(temperatureK-trefK)/(trefK*R*temperatureK))
  kend = (1+exp((trefK*delS-Hd)/(trefK*R)))/(1+exp((temperatureK*delS-Hd)/(temperatureK*R)))
  
  multiplier = kbeg * kend
  
  multiplier
  
}