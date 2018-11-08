# calculate omega for use in calculating c

calc_omega_c = function(temp, z, vpdo, cao, jv, theta){
	
	patm = calc_patm(z)
	vpd = calc_vpd(temp, z, vpdo)
	vpd_pa = vpd * 1000
	ca_pa = cao * 1e-6 * patm
	
	chi = calc_chi(temp, z, vpdo, cao)
	ci = chi * ca_pa
	gammastar = calc_gammastar_pa(temp, z)
	km_pa = calc_km_pa(temp, z) 
	
	mc = (ci - gammastar) / (ci + km_pa) 
	m = (ci - gammastar) / (ci + 2*gammastar)
	
	v = (jv * m) / (8*theta*mc)
	
	omega_c = (2 - (4 * theta * v)) / ((1 / v) - 2)
	
	omega_c	
	
}


