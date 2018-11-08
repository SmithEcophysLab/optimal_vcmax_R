# calculate c for Vcmax prediction

calc_c = function(temp, z, vpdo, cao, jv25, theta){
	
	#theta = 0.9
	
	patm = calc_patm(z)
	ca_pa = cao * 1e-6 * patm
	omega_c = calc_omega_c(temp, z, vpdo, cao, jv25, theta)
	gammastar_pa = calc_gammastar_pa(temp, z) # gamma star Pa
	chi = calc_chi(temp, z, vpdo, cao) # ci/ca		
	ci = chi * ca_pa # Pa
	m = ((ci - gammastar_pa)/(ci + (2 * gammastar_pa)))
	c_quad = -(theta)
	b_quad = 1
	a_quad1 = (omega_c + (1 - 2*theta))^2
	a_quad2 = a_quad1 / (1-theta)
	a_quad3 = a_quad2 + 4*theta
	a_quad = -1/a_quad3
	c_star = (-b_quad + sqrt(b_quad^2 - (4 * a_quad * c_quad)))/(2 * c_quad)
	
	c = (c_star * m) / 4
	
}


