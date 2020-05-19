# calculate nstar (unitless relative viscosity of h2o at tempeature relative to 25°C)

calc_nstar = function(temp, z){ # temp in °C and z in m
	
	patm = calc_patm(z)
	
	# viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa) 
    ns      = calc_viscosity_h2o( temp, z )  # Pa s 
    ns25    = calc_viscosity_h2o( 25, z )  # Pa s 
    nstar = ns / ns25                       # (unitless)
    
    nstar
	
}
