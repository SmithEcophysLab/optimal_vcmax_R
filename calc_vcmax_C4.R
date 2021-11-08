# calculate optimal vcmax for C4 plants

########################
## variable key
########################

# tg_c: acclimated temperature (degC)
# z: elevation (m)
# vpdo: vapor pressure deficit at sea level (kPa)
# cao: atmospheric CO2 at sea level (umol mol-1)
# oao: atmospheric O2 at sea level (ppm)
# paro: photosynthetically active radiation at sea level (µmol m-2 s-1)
# q0: quantum efficiency of photosynthetic electron transport (mol/mol)
# theta: curvature of the light response of electron transport (unitless)
# R: universal gas constant (J mol-1 K-1)
# patm: atmospheric pressure (Pa)
# ca: atmospheric CO2 at z (Pa)
# km: Michaelis-Menten constant for Rubisco (Pa)
# gammastar: CO2 compensation point (Pa)
# chi: leaf intercellular to atmospheric CO2 ratio (ci/ca) (unitless)
# vpd: vapor pressure deficit at z (kPa)
# par: photosynthetically active radiation at z (µmol m-2 s-1)
# ci: leaf intercellular CO2 concentation (Pa)
# m: CO2 limiation of electron transport rate limited photosynthesis (Pa)
# mc: CO2 limiation of Rubisco carboxylation rate limited photosynthesis (Pa)
# c: constant describing cost of maintaining electron transport (unitless)
# omega: omega term from Smith et al.
# omega_star: omega_star term from Smith et al.
# vcmax_star: maximum rate of Rubisco carboxylation without temperature correction (µmol m-2 s-1) 
# vcmax_prime: optimal maximum rate of Rubisco carboxylation at tg (µmol m-2 s-1)
# jvrat: ratio of the maximum rate of electron transport to the maximum rate of Rubisco carbocylation at tg (unitless)
# jmax_prime: optimal maximum rate of electron transport at tg (µmol m-2 s-1)
# vcmax25: optimal maximum rate of Rubisco carboxylation at a temperature of 25 C (µmol m-2 s-1) 
# vpmax25: ...at a temperature of 25 C (µmol m-2 s-1)
# jmax25: optimal maximum rate of electron transport at a temperature of 25 C (µmol m-2 s-1)
# lma: leaf mass area
# grossphoto: gross photosynthesis
# resp: respiration
# netphoto: net photosynthesis
# nrubisco: leaf N in rubisco
# nbioe: leaf N in bioenergetics
# npep: leaf N in rubisco from predicted vpmax with PEP-specific constants
# nstructure: nitrogen in structural tissue from lma
# nall: all leaf N predictions
# nphoto: leaf N used for photosynthesis
# nrubisco_frac: fraction of leaf N in rubisco out of all leaf N
# nphoto_frac: fraction of leaf N for photosynthesis out of all leaf N


# All equation numbers refer to equations presented in Smith et al., 2019 

# libraries
# install.packages('R.utils')
library(R.utils)

# load necessary functions
sourceDirectory('functions')

calc_vcmax_C4 <- function(tg_c = 25, z = 0, vpdo = 1, cao = 400, oao = 209460,
                               paro = 800, q0 = 0.257, theta = 0.85, f = 0.5){
  
  # constants
  R <- 8.314
  c <- 0.05336251
  leakiness <- 0.01
  
  # environmental terms
  patm <- calc_patm(z)
  par <- calc_par(paro, z)
  vpd <- calc_vpd(tg_c, z, vpdo)
  ca <- cao * 1e-6 * patm
  oa <- oao * 1e-6 * patm
  
  # K and Gamma* model terms
  km <- calc_km_pa(tg_c, z)                
  gammastar <- calc_gammastar_pa(tg_c, z)   
  
  # calc chi
  chi <- calc_chi_c4(cao, tg_c, vpd, z)
  # calc ci ( = cm)
  ci <- ca * chi
  cm <- ci
  oi <- oa * chi
    
  # Light Limited Photosynthesis
  m <- (ci - gammastar) / (ci + 2 * gammastar)
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
  vcmax <- (q0 * par * m * omega_star / (8 * theta)) * ((cbs + kr * (1 + obs/ko)) / (cbs - gammastar)) # Eqn. 2.47
  mc <- ((cbs + kr * (1 + obs/ko)) / (cbs - gammastar))
  Ac <- vcmax * ((cbs - gammastar) / (kr * (1 + obs/ko) + cbs)) # Eqn. 2.4
    
  # calculate vcmax and jmax at temperature = 25 C
  vcmax25 <- vcmax / calc_tresp_mult(tg_c, tg_c, 25)
  vpmax25 <- vpmax / calc_tresp_mult(tg_c, tg_c, 25) # using vcmax parameters - find better solution
  jmax25 <- jmax / calc_jmax_tresp_mult(tg_c, tg_c, 25)
  
  # estimate LMA 
  lma <- calc_lma(f = f, par = paro, temperature = tg_c, vpd_kpa = vpdo, z = z, co2 = 400)
  
  # calculate leaf N in rubisco from predicted vcmax
  nrubisco <- fvcmax25_nrubisco(vcmax25)
  
  # calculate leaf N in bioenergetics from predicted jmax
  nbioe <- fjmax25_nbioe(jmax25)
  
  # calculate leaf N in rubisco from predicted vpmax with PEP-specific constants
  npep <- fvpmax25_npep(vpmax25)
  
  # calculate nitrogen in structural tissue from lma 
  nstructure <- flma_nstructure(lma)
  
  # sum all leaf N predictions
  nall <- nrubisco + nbioe + nstructure + npep
  
  # calculate leaf N used for photosynthesis
  nphoto <- nrubisco + nbioe + npep
  
  # calculate the fraction of leaf N in rubisco out of all leaf N
  nrubisco_frac <- nrubisco / nall
  
  # calculate the fraction of leaf N for photosynthesis out of all leaf N
  nphoto_frac <- nphoto / nall
  
  # output
  results <- data.frame("tg_c" = tg_c,
                        "z" = z,
                        "vpdo" = vpdo,
                        "cao" = cao,
                        "paro" = paro,
                        "q0" = q0,
                        "theta" = theta,
                        "c" = c,
                        "par" = par,
                        "patm" = patm,
                        "ca" = ca,
                        "oa" = oa,
                        "vpd" = vpd,
                        "kp" = kp,
                        "kr" = kr,
                        "chi" = chi,
                        "ci" = ci,
                        "Leakage" = leakage,
                        "cbs" = cbs,
                        "chi_bs" = chi_bs,
                        "obs" = obs,
                        "Al" = Al,
                        "Ap" = Ap,
                        "Ac" = Ac,
                        "km" = km,
                        "gammastar" = gammastar,
                        "omega" = omega,
                        "m" = m,
                        "mc" = mc,
                        "omega_star" = omega_star,
                        "vcmax" = vcmax,
                        "vpmax" = vpmax,
                        "jmax" = jmax,
                        "vcmax25" = vcmax25,
                        "vpmax25" = vpmax25,
                        "jmax25" = jmax25,
                        "lma" = lma,
                        "nrubisco" = nrubisco,
                        "nbioe" = nbioe,
                        "npep" = npep,
                        "nstructure" = nstructure,
                        "nall" = nall,
                        "nphoto" = nphoto,
                        "nrubisco_frac" = nrubisco_frac,
                        "nphoto_frac" = nphoto_frac)

  return(results)
  # return(list(tg_c, z, vpdo, cao, paro, q0, theta, c, par, patm, ca, oa, vpd, kp, kr, chi, ci, leakage, cbs, 
  #             chi_bs, obs, Al, Ap, Ac, km, gammastar, omega, m, mc, omega_star, vcmax, vpmax, jmax, vcmax25, 
  #             vpmax25, jmax25, lma, nrubisco, nbioe, npep, nstructure, nall, nphoto, nrubisco_frac, nphoto_frac))
}
