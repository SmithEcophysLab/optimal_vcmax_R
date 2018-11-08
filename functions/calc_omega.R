# calculate omega for vcmax calculation

calc_omega = function(theta, c, m){ # paro in µmol m-2 s-1 and z in m
	
	cm = ((4 * c) / m) # simplification term for omega calculation
	v = 1/(cm * (1 - (theta * cm))) - 4 * theta # simplification term for omega calculation
	
	# account for non-linearities at low m values
	capP = (((1/1.4) - 0.7)^2 / (1-theta)) + 3.4
	aquad = -1
	bquad = capP
	cquad = -(capP * theta)
	m_star = (4 * c) / polyroot(c(aquad, bquad, cquad))
	
	omega = ifelse(m < Re(m_star[1]), -(1 - (2 * theta)) - sqrt((1 - theta) * v), -(1 - (2 * theta)) + sqrt((1 - theta) * v))
	
	omega

}