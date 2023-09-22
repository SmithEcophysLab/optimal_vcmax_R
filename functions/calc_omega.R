# calculate omega for vcmax calculation

calc_omega = function(theta, c, m){ # paro in µmol m-2 s-1 and z in m
	
	cm = ((4 * c) / m) # simplification term for omega calculation
	v = 1/(cm * (1 - (theta * cm))) - 4 * theta # simplification term for omega calculation
	
	# account for non-linearities at low m values

	omega = ifelse(8*c > m & 4*c < m & m/c < 8*theta, -(1 - (2 * theta)) - sqrt((1 - theta) * v), -(1 - (2 * theta)) + sqrt((1 - theta) * v))

	omega

}