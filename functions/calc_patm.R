# calculate atmospheric pressure (Pa) from elevation (m)

calc_patm = function(z) {
	
	kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
	kTo = 298.15   # base temperature, K (Prentice, unpublished)
	kL = 0.0065    # temperature lapse rate, K/m (Allen, 1973)
	kG = 9.80665   # gravitational acceleration, m/s**2 (Allen, 1973)
	kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
	kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
	
	patm = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
	
	patm
}


