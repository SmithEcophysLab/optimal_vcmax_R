# calculate the density of h2o (kg m-3) given tmeperature and pressure

calc_density_h2o = function(temp, z){ # temp in °C and z in m
	
	patm = calc_patm(z)
	
	# Calculate lambda, (bar cm**3)/g:
    my_lambda = 1788.316 + 1.55053*temp + -0.4695911*temp*temp + (3.096363e-3)*temp*temp*temp + -(7.341182e-6)*temp*temp*temp*temp

    # Calculate po, bar
    po = 5918.499 + 58.05267*temp + -1.1253317*temp*temp + (6.6123869e-3)*temp*temp*temp + -(1.4661625e-5)*temp*temp*temp*temp

    # Calculate vinf, cm**3/g
    vinf = 0.6980547 + -(7.435626e-4)*temp + (3.704258e-5)*temp*temp + -(6.315724e-7)*temp*temp*temp + (9.829576e-9)*temp*temp*temp*temp + -(1.197269e-10)*temp*temp*temp*temp*temp + (1.005461e-12)*temp*temp*temp*temp*temp*temp + -(5.437898e-15)*temp*temp*temp*temp*temp*temp*temp + (1.69946e-17)*temp*temp*temp*temp*temp*temp*temp*temp + -(2.295063e-20)*temp*temp*temp*temp*temp*temp*temp*temp*temp

    # Convert pressure to bars (1 bar = 100000 Pa)
    pbar = (1e-5)*patm
    
    # Calculate the specific volume (cm**3 g**-1):
    vau = vinf + my_lambda/(po + pbar)

    # Convert to density (g cm**-3) -> 1000 g/kg; 1000000 cm**3/m**3 -> kg/m**3:
    density_h2o = (1e3/vau)
    
    density_h2o
	
}
