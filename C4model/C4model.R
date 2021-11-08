# C4model: Predicts acclimated 
# Helen Scott
# Last Updated: 12/18/2020
#
# Arguments
## tg_c: acclimated temperature (degC)
## z: elevation (m)
## vpdo: vapor pressure deficit at sea level (kPa)
## cao: atmospheric CO2 at sea level (ppm)
## oao: atmospheric O2 at sea level (ppm)
## paro: photosynthetically active radiation at sea level (?mol m-2 s-1)
## theta: curvature of the light response of electron transport (unitless)
## R: universal gas constant (J mol-1 K-1)
#
# Returns
## tg_c: acclimated temperature (degC)
## par: photosynthetically active radiation at z (?mol m-2 s-1)
## ca: atmospheric CO2 at z (Pa)
## z: elevation (m)
## vpd: vapor pressure deficit at z (kPa)
## q0: quantum efficiency of photosynthetic electron transport at tg_c (mol/mol) # NOT AN INPUT ANYMORE
## kp: michaelis menten coefficient for PEPc (Pa)
## kr: Michaelis-Menten constant for Rubisco carboxylation (Pa)
## chi: leaf intercellular to atmospheric CO2 ratio (ci/ca) (unitless)
## ci: leaf intercellular CO2 concentation (Pa)
## Leakage: Rate of leakage of CO2 from bundle sheath to mesophyll (?mol m-2 s-1)
## cbs: bundle sheath CO2 pressure (Pa)
## obs: bundle sheath O2 pressure (Pa)
## jmax: optimal maximum rate of electron transport at tg (?mol m-2 s-1)
## vpmax: maximum rate of PEPc (?mol m-2 s-1)
## vcmax: maximum rate of Rubisco carboxylation (?mol m-2 s-1)
## jv ratio: ratio of jmax to vcmax (unitless)
## Al: light-limited photosynthesis (?mol m-2 s-1)
## Ap: PEPc-limited photosynthesis (?mol m-2 s-1)
## Ac: Rubisco-limited photosynthesis (?mol m-2 s-1)

C4model <- function(tg_c = 25, z = 0, vpdo = 1, cao = 400, oao = 209460, 
                  paro = 800, theta = 0.85, leakiness = 0.01, 
                  R = 8.314){
  
  # environmental terms
  patm <- calc_patm(z)
  par <- calc_par(paro, z)
  vpd <- calc_vpd(tg_c, z, vpdo)
  ca <- cao * 1e-6 * patm
  oa <- oao * 1e-6 * patm
  
  # Calculate Gamma star
  
  gamma_star <- calc_gammastar_pa_c4(tg_c, z) # pa
  
  # calc chi
  chi_m <- calc_chi_c4(cao, tg_c, vpd, z)
  # calc ci ( = cm)
  ci <- ca * chi_m
  cm <- ci
  oi <- oa * chi_m
    
  # Light Limited Photosynthesis
  m <- (ci - gamma_star) / (ci + 2 * gamma_star)
  omega <- calc_omega(theta = theta, c = 0.01, m = m) # Eq. S4
  omega_star <- (1 + (omega) - sqrt((1 + (omega))^2 - (4 * theta * omega)))  # Eq. 18
  # calculate q0 using Bernacchi et al. (2003) temperature response (set to 0.257 at 25C)
  q0 <- -0.0805 + (0.022 * tg_c) - (0.00034 * tg_c * tg_c)
  Al <- q0 * par * m * omega_star / (8 * theta) # Eqn. 2.2
  jmax <- q0 * par * omega
    
  # calc kp
  kp <- calc_kp_temp_pa(tg_c, z) # Eqn. 2.43
  # calc kr
  kr <- calc_kc_temp_pa(tg_c, z) # Eqn. 2.48
  # calc ko
  ko <- calc_ko_temp_pa(tg_c, z)
   
  # calc vpmax
  vpmax <- ((kp + cm)/cm) * (q0 * par * m * omega_star / (8 * theta)) # Eqn. 2.42
  Ap <- vpmax * (cm / (cm + kp))
  
  # calc cbs
  leakage <- leakiness * Al
  cbs <- calc_cbs(cm, leakage) # Eqn. 2.41
  chi_bs <- cbs / ca
  # calc obs
  obs <- oi
  
  # calc vcmax
  vcmax <- (q0 * par * m * omega_star / (8 * theta)) * ((cbs + kr * (1 + obs/ko)) / (cbs - gamma_star)) # Eqn. 2.47
  Ac <- vcmax * ((cbs - gamma_star) / (kr * (1 + obs/ko) + cbs)) # Eqn. 2.4
  
  results <- data.frame("tg_c" = tg_c,
                        "par" = par,
                        "cao" = cao,
                        "ca" = ca,
                        "oa" = oa,
                        "z" = z,
                        "vpd" = vpd,
                        "q0" = q0,
                        "kp" = kp,
                        "kr" = kr,
                        "chi_m" = chi_m,
                        "ci" = ci,
                        "Leakage" = leakage,
                        "cbs" = cbs,
                        "chi_bs" = chi_bs,
                        "obs" = obs,
                        "jmax" = jmax,
                        "vpmax" = vpmax,
                        "vcmax" = vcmax,
                        "jvc_ratio" = jmax/vcmax,
                        "jvp_ratio" = jmax/vpmax,
                        "vcvp_ratio" = vcmax/vpmax,
                        "Al" = Al,
                        "Ap" = Ap,
                        "Ac" = Ac,
                        "Gamma" = gamma_star)
  return(results)
}
