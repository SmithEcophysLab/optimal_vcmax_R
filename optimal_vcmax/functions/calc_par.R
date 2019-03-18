# calculate par in µmol m-2 s-1, from Wang (2016) New Phytologist

calc_par = function(paro, z){ # paro in µmol m-2 s-1 and z in m
	
	z_km = z/1000 
	
	par = paro * (1 + 0.027*z_km)
	
	par

}

