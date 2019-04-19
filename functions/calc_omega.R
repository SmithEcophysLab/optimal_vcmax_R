# calculate omega for vcmax calculation

calc_omega = function(theta, c, m){ # paro in µmol m-2 s-1 and z in m
	
	cm = ((4 * c) / m) # simplification term for omega calculation
	v = 1/(cm * (1 - (theta * cm))) - 4 * theta # simplification term for omega calculation
	
	# account for non-linearities at low m values
	cquad = -1
	bquad = 3.4
	aquad = -(3.4 * theta)
	quad_root = (-bquad + sqrt(bquad^2 - (4 * aquad * cquad))) / (2 * aquad)
	m_star = (4 * c) / quad_root

	omega = ifelse(m < m_star, -(1 - (2 * theta)) - sqrt((1 - theta) * v), -(1 - (2 * theta)) + sqrt((1 - theta) * v))
	  
	omega

}