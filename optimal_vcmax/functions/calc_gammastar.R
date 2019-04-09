# calculate gammastar (Pa)

calc_gammastar_pa = function(temp, z) {
	
	patm = calc_patm(z)
	rat = calc_patm(z) / calc_patm(0)
	
	#gammastar25 = 42.75  # ppm
	gammastar25 = 4.332 * rat  # Pa
	Hgm=37830 # J mol-1
	R = 8.314        # J K-1 mol-1
	O2 = 2.09476e5 # ppm
	O2_0 = O2 * 1e-6 * calc_patm(0)
	O2_z = O2 * 1e-6 * calc_patm(z)
	
	temp_k = 273.15+ temp
	
	gStar_pa = gammastar25*exp((Hgm/R)*(1/298.15-1/temp_k))
	
	gStar_pa
	
}
