# calcualte viscosity of water
# from Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V. Sengers, M. J. Assael, ..., K. Miyagawa (2009) New international formulation for the viscosity of H2O, J. Phys. Chem. Ref. Data, Vol. 38(2), pp. 101-125.

calc_viscosity_h2o = function(temp, z){ #temp in °C and z in m
	
	tk_ast  = 647.096    # Kelvin
    rho_ast = 322.0      # kg/m**3
    mu_ast  = 1e-6       # Pa s
	
	patm = calc_patm(z)
	
	rho = calc_density_h2o(temp, z) # density of water (kg m-3)
	
	# Calculate dimensionless parameters:
    tbar = (temp + 273.15)/tk_ast
    tbarx = tbar**(0.5)
    tbar2 = tbar**2
    tbar3 = tbar**3
    rbar = rho/rho_ast
    
    # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
    mu0 = 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3
    mu0 = 1e2*tbarx/mu0
    
    # Create Table 3, Huber et al. (2009):

    h_array1 = c(0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0)
    h_array2 = c(0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573)
    h_array3 = c(-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0)
    h_array4 = c(0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0)
    h_array5 = c(-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0)
    h_array6 = c(0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0)
    h_array7 = c(0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264)
    
    h_array = rbind(h_array1, h_array2, h_array3, h_array4, h_array5, h_array6, h_array7)
    
    # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
    
    mu1 = 0.0
    ctbar = (1.0/tbar) - 1.0
    for (i in 1:6){
    	
    		coef1 = ctbar**(i-1)
    		coef2 = 0.0
    	
    		for (j in 1:7){
    			
    			coef2 = coef2 + h_array[j,i] * (rbar - 1.0)**(j-1)
    			    		
    		}
    		
    		mu1 = mu1 + coef1 * coef2 
    	
    }
    
    mu1 = exp( rbar * mu1 )

    # Calculate mu_bar (Eq. 2, Huber et al., 2009)
    #   assumes mu2 = 1
    mu_bar = mu0 * mu1

    # Calculate mu (Eq. 1, Huber et al., 2009)
    viscosity_h2o = mu_bar * mu_ast    # Pa s
    viscosity_h2o
	
}


