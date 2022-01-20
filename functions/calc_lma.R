# From Wang et al. 2021: https://doi.org/10.1101/2021.02.07.430028 

#### 1. define function and parameters in photosynthesis model
## functions
{
  # calculate air pressure in Pa
  calc_patm <- function( elv ){
    #-----------------------------------------------------------------------
    # Input:    - elevation, m (elv)
    # Output:   - float, atmospheric pressure at elevation 'elv', Pa (patm)
    # Features: Returns the atmospheric pressure as a function of elevation
    #           and standard atmosphere (1013.25 hPa)
    # Depends:  - connect_sql
    #           - flux_to_grid
    #           - get_data_point
    #           - get_msvidx
    # Ref:      Allen et al. (1998)
    #-----------------------------------------------------------------------
    
    # Define constants:
    kPo <- 101325   # standard atmosphere, Pa (Allen, 1973)
    kTo <- 298.15   # base temperature, K (Prentice, unpublished)
    kL <- 0.0065    # temperature lapse rate, K/m (Allen, 1973)
    kG <- 9.80665   # gravitational acceleration, m/s^2 (Allen, 1973)
    kR <- 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
    kMa <- 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
    
    # Convert elevation to pressure, Pa:
    patm <- kPo*(1.0 - kL*elv/kTo)**(kG*kMa/(kR*kL))
    
    return (patm)
  }
  
  # calculate K (MM coefficient of Rubisco) in Pa
  calc_k <- function(temp, patm) {
    #-----------------------------------------------------------------------
    # Input:    - float, air temperature, deg C (temp)
    #           - float, atmospheric pressure, Pa (patm)
    # Output:   float, Pa (mmk)
    # Features: Returns the temperature & pressure dependent Michaelis-Menten
    #           coefficient, K (Pa).
    # Ref:      Bernacchi et al. (2001), Improved temperature response 
    #           functions for models of Rubisco-limited photosynthesis, 
    #           Plant, Cell and Environment, 24, 253--259.
    #-----------------------------------------------------------------------
    
    # Define constants
    kc25 <- 39.97      # Pa, assuming 25 deg C & 98.716 kPa
    ko25 <- 2.748e4    # Pa, assuming 25 deg C & 98.716 kPa
    dhac <- 79430      # J/mol
    dhao <- 36380      # J/mol
    kR   <- 8.3145     # J/mol/K
    kco  <- 2.09476e5  # ppm, US Standard Atmosphere
    
    vc <- kc25*exp(dhac*(temp - 25.0)/(298.15*kR*(temp + 273.15)))
    vo <- ko25*exp(dhao*(temp - 25.0)/(298.15*kR*(temp + 273.15)))
    k  <- vc*(1 + kco*(1e-6)*patm/vo)
    
    return(k)
    
  }
  
  # calculate Gstar (CO2 compensation point) in Pa
  calc_gstar_gepisat <- function( temp ) {
    #-----------------------------------------------------------------------
    # Input:    float, air temperature, degrees C (tc)
    # Output:   float, gamma-star, Pa (gs)
    # Features: Returns the temperature-dependent photorespiratory 
    #           compensation point, Gamma star (Pascals), based on constants 
    #           derived from Bernacchi et al. (2001) study.
    # Ref:      Bernacchi et al. (2001), Improved temperature response 
    #           functions for models of Rubisco-limited photosynthesis, 
    #           Plant, Cell and Environment, 24, 253--259.
    #-----------------------------------------------------------------------
    
    # Define constants
    gs25 <- 4.220    # Pa, assuming 25 deg C & 98.716 kPa)
    dha  <- 37830    # J/mol
    kR   <- 8.3145   # J/mol/K
    
    gs <- gs25 * exp( dha * ( temp - 25.0 ) / ( 298.15 * kR * ( temp + 273.15 ) ) )
    
    return( gs )
    
  }
  
  # conver CO2 from ppm to Pa
  co2_to_ca <- function( co2, patm ){
    #-----------------------------------------------------------------------
    # Input:    - float, annual atm. CO2, ppm (co2)
    #           - float, monthly atm. pressure, Pa (patm)
    # Output:   - ca in units of Pa
    # Features: Converts ca (ambient CO2) from ppm to Pa.
    #-----------------------------------------------------------------------
    ca   <- ( 1.e-6 ) * co2 * patm         # Pa, atms. CO2
    return( ca )
  }
}

calc_lma <- function(deciduous = 'yes', # yes or no for deciduous species
                     f = 0.5, # fraction of time in growing season
                     par = 300, # light absorbed in µmol m-2 s-1
                     temperature = 25, # temperature in C
                     k1 = 30, # scaling factor, unit: g biomass / molC
                     u = 768, # Xu, X. et al. Variations of leaf longevity in tropical moist forests predicted by a trait-driven carbon optimality model. Ecology letters 20, 1097-1106 (2017).
                     vpd_kpa = 1, # vpd in kpa
                     z = 0, # elevation in m
                     co2 = 400, # co2 in ppm
                     phi_intercept = 0.352, # intercept of the phi0 temperature response
                     alpha = 0.8
                     ){
  
  # constants
  Ha.v <- 72000
  R <- 8.3145 # J mol-1 K-1
  beta <- 146 
  
  # calculate pressure
  pressure = calc_patm(z)
  
  # calculate ca
  ca <- co2_to_ca(co2, pressure)
  
  # calculate Gstar
  Gstar <- calc_gstar_gepisat(temperature)
  
  # calculate M-M terms
  K <- calc_k(temperature, pressure)
  f1 <- exp(-0.0227*(temperature-25)) # the viscosity of water relative to its value at 25˚C
  g1 <- sqrt(beta*(Gstar+K)/(1.6*f1))
  x <- Gstar/ca + (1-Gstar/ca)*g1/(g1+sqrt(vpd_kpa*1000))
  ci <- x*ca
  mc <- (ci - Gstar)/(ci + K)
  m <- (ci - Gstar)/(ci + 2*Gstar)
  
  # calculate phi0
  # phi0 <- (phi_intercept+0.021*temperature-3.4*10^(-4)*(temperature)^2)/8
  phi0 = 0.25
  
  # calculate hT
  dS.v<-668.39-1.07*(temperature + 273.15)
  hT <- ( exp(Ha.v*(temperature - 25)/(298.15*R*(temperature + 273.15))) ) * 
    ( (1 + exp((298.15*dS.v - 200000)/(298.15*R)))/(1 + exp(((temperature + 273.15)*dS.v - 200000)/((temperature + 273.15)*R))) )
  
  if(deciduous == 'yes'){
    # deciduous LMA
    xTde <- hT*mc/(phi0*m)
    Cde <- log(k1) + log(365) - log(u) - log(xTde)
    # LMA <- exp(log(par) - 0.052*temperature - 0.6*log(alpha) + log(f) + 1.7)
    LMA <- exp(log(f) + log(par) - 0.052*temperature + Cde)
  }
  
  else{
    # evergreen LMA
    LMA <- exp(0.5*log(par) - 0.013*temperature - 0.12*log(alpha) + 0.25*log(f) + 3.5)
    #LMA <- exp(0.5*log(par) - 0.013*temperature + 0.25*log(f) + 3.5)
  }

  return(LMA)
  
}

return_lma <- function(given_lma, lma, f, par, temperature, vpd_kpa, z, co2){
  ifelse(given_lma == "yes",
         return(lma),
         return(calc_lma(f = f, par = par, temperature = temperature, vpd_kpa = vpd_kpa, z = z, co2 = co2)))
}

# given_lma <- function(lma){
#   ifelse(!is.na(lma), "yes", "no")
# }

## to do: need par adjustments based on elevation
## need to integrate with other code
## need to figure out how to estimate growing season length (f)
