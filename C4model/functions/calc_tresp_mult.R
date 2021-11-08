# calculate the temperature response multiplier for Vcmax and Jmax (Kattge and Knorr 2007)

calc_tresp_mult = function(tleaf, tmean, tref){
	
	temp <- tleaf + 273.15
	Ha <- 71513 # Activation Energy (J mol^-1) Massad 2007
	Hd <- 200000 # Deactivaiton energy (J mol^-1) Massad 2007
	adelS <- 668.39
	bdelS <- -1.07

	trefK <- tref + 273.15
	R <- 8.314 # Universal Gas constant (J mol^-1 K^-1)
	
	kbeg <- exp(Ha*(temp-trefK)/(trefK*R*temp))
	kend <-((1+exp((trefK*(adelS+bdelS*tmean)-Hd)/(trefK*R)))/(1+exp((temp*(adelS+bdelS*tmean)-Hd)/(temp*R))))
	
	ktotal <- kbeg*kend # Equation 20 in Smith 2019
	return(ktotal)
}
