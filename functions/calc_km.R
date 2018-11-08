# calcualte the Michaelis-Menton coefficient (Pa) for Rubisco from temperature

calc_km_pa = function(temp, z) {
	
	patm = calc_patm(z) 
	rat = patm / calc_patm(0)
	
	R = 8.314        
    O2 = 2.09476e5      
    Kc25 = 41.03 * rat 
    Ko25 = 28210 * rat 
    Hkc = 79430  
    Hko = 36380 
    
    temp_k = 273.15 + temp

    Kc_pa =Kc25 * exp(Hkc * ((temp_k - 298.15) / (298.15 * R * temp_k)))
    Ko_pa =Ko25* exp(Hko * ((temp_k - 298.15) / (298.15 * R * temp_k)))
    
    O2_pa = O2 * (1e-6) * patm 
    
    Km_pa = Kc_pa * (1 + O2_pa/Ko_pa)
    
    Km_pa 
	
}



