# calculate optimal vcmax as in Smith et al. Photosynthetic capacity is optimized to the 
# environment. Ecology Letters. 

########################
## variable key
########################

# pathway: photosynthetic pathway, either "c3" or "c4"
# tg_c: acclimated temperature (degC)
# z: elevation (m)
# vpdo: vapor pressure deficit at sea level (kPa)
# cao: atmospheric CO2 at sea level (umol mol-1)
# oao: atmospheric O2 at sea level (ppm)
# paro: photosynthetically active radiation at sea level (µmol m-2 s-1)
# q0_resp: yes or no, use the q0 response curve calculation
# q0_int: intercept for the q0 response curve calculation
# q0: quantum efficiency of photosynthetic electron transport (mol/mol)
# theta: curvature of the light response of electron transport (unitless)
# given_chi: yes or no based on if chi is known (yes) or calculated (no)
# chi: leaf intercellular to atmospheric CO2 ratio (ci/ca) (unitless)
# f: fraction of year in growing season
# lma: leaf mass area (g m-2)
# given_lma: yes or no based on if lma is known (yes) or calculated (no)
# R: universal gas constant (J mol-1 K-1)
# patm: atmospheric pressure (Pa)
# ca: atmospheric CO2 at z (Pa)
# km: Michaelis-Menten constant for Rubisco (Pa)
# gammastar: CO2 compensation point (Pa)
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
# parea: area-based leaf photosynthetic phosphorus (g m-2)
# nue = nitrogen use efficiency
# wue = water use efficiency
# pue = phosphorus use efficiency


# All equation numbers refer to equations presented in Smith et al., 2019 

# libraries
# install.packages('R.utils')
library(R.utils)

# load necessary functions
sourceDirectory('functions', modifiedOnly = FALSE)

calc_optimal_vcmax <- function(pathway = "C3", deciduous = "yes", tg_c = 25, z = 0, vpdo = 1, cao = 400, oao = 209460,
                               paro = 800, beta = 146, theta = 0.85, chi = NA, q0 = 0.257, q0_resp = "yes", 
                               q0_int = -0.0805, lma = NA, f = 0.5){
  
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

  if(pathway == "C3"){
    # C3

    gammastar <- calc_gammastar_pa(tg_c, z)
    
    # Coordination and least-cost hypothesis model terms
    given_chi <- given_chi(chi) # returns yes or no based on presence of chi value
    # calculate chi if unknown
    chi <- return_chi(given_chi = given_chi, chi = chi, temp = tg_c, z = z, vpdo = vpdo, cao = cao, beta = beta)
    ci <- chi * ca # Pa
    mc <- ((ci - gammastar) / (ci + km))
    m <- ((ci - gammastar)/(ci + (2 * gammastar)))
    omega <- calc_omega(theta = theta, c = c, m = m)
    omega_star <- (1 + (omega) - sqrt((1 + (omega))^2 - (4 * theta * omega)))
    
    # calculate q0
    if(q0_resp == "yes"){
      # Bernacchi et al. (2003) temperature response (set to 0.257 at 25C)
      q0 = q0_int + (0.022 * tg_c) - (0.00034 * tg_c * tg_c)
    }else{
      q0
    }

    # calculate vcmax and jmax
    vcmax <- ((q0 * par * m) / mc) * (omega_star / (8 * theta))
    jvrat <- ((8 * theta * mc * omega) / (m * omega_star))
    jmax <- jvrat * vcmax

    vpmax <- 0
    vpmax25 <- 0
    leakage <- NA
    chi_bs <- NA
    kp <- NA
    kr <- NA
    ko <- NA
    cbs <- NA
    chi_bs <- NA
    obs <- NA
    Al <- q0 * par * m * omega_star / (8 * theta)
    Ap <- 0
    Ac <- vcmax * mc
    
    # calc rd
    rd <- vcmax * 0.018 # Ren et al. (2024)
    
    # calc anet
    Anet <- Ac - rd
    
    }else{
      # C4
      
      gammastar <- calc_gammastar_pa_c4(tg_c, z)
      
      # Coordination and least-cost hypothesis model terms
      given_chi <- given_chi(chi) # returns yes or no based on presence of chi value
      # calculate chi if unknown
      chi <- return_chi_c4(given_chi = given_chi, chi = chi, ca = ca, temp = tg_c, vpd = vpdo, z = z, beta = beta)
      ci <- ca * chi
      cm <- ci
      oi <- oa * chi
      m <- (ci - gammastar) / (ci + 2 * gammastar)
      omega <- calc_omega(theta = theta, c = 0.01, m = m) # Eq. S4
      omega_star <- (1 + (omega) - sqrt((1 + (omega))^2 - (4 * theta * omega)))  # Eq. 18
      
      # calculate q0
      if(q0_resp == "yes"){
        # Bernacchi et al. (2003) temperature response (set to 0.257 at 25C)
        q0 = q0_int + (0.022 * tg_c) - (0.00034 * tg_c * tg_c) 
      }else{
        q0
      }
      Al <- q0 * par * m * omega_star / (8 * theta) # Eqn. 2.2
      jmax <- q0 * par * omega

      # calc kp, kr, and ko
      kp <- calc_kp_temp_pa(tg_c, z) # Eqn. 2.43
      kr <- calc_kc_temp_pa(tg_c, z) # Eqn. 2.48
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
      
      # calc rd
      rd <- vcmax * 0.018 # Ren et al. (2024)
      
      # calc anet
      Anet <- Ac - rd

    }
  
  # calculate vcmax and jmax at temperature = 25 C
  vcmax25 <- vcmax / calc_tresp_mult(tg_c, tg_c, 25)
  jmax25 <- jmax / calc_jmax_tresp_mult(tg_c, tg_c, 25)
  
  # calculate vpmax at temperature = 25 C
  vpmax25 <- vpmax / calc_tresp_mult(tg_c, tg_c, 25) # using vcmax parameters - find better solution
  
  # calculate rd at temperature = 25 C
  rd25 <- rd /calc_rd_tresp_mult(tg_c, tg_c, 25)
  
  # LMA
  given_lma <- ifelse(!is.na(lma), "yes", "no") # returns yes or no based on presence of LMA value
  # calculate LMA if unknown
  #lma <- return_lma(given_lma, lma, deciduous = deciduous, f = f, par = paro, temperature = tg_c, vpd_kpa = vpdo, z = z, co2 = cao)
  if(given_lma == "yes"){
    lma <- lma
  }else{
    lma <- calc_lma(deciduous = deciduous, f = f, par = par, temperature = tg_c, 
                    beta = beta, R = R, m = m, mc = mc)
  }
  
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
  
  # calculate stomatal conductance
  gsc <- ((Ac / cao) / (1 - chi))
  gsw <- gsc * 1.6
  
  # calculate leaf n metrics
  narea <- nall
  nmass <- narea * (1/lma)
  
  # calculate leaf P (g m-2)
  parea <- Al * (1/1e6) * (1/11.48) * 5 * 92.91  
  
  # calculate efficiency metrics
  wue <- Anet/gsw
  nue <- Anet/nall
  pue <- Anet/parea

  
	# output
  results <- data.frame("pathway" = pathway,
                        "deciduous" = deciduous,
	                      "tg_c" = tg_c,
	                      "z" = z,
	                      "vpdo" = vpdo,
	                      "cao" = cao,
                        "oao" = oao,
	                      "paro" = paro,
                        "q0_response" = q0_resp,
                        "q0_intercept" = q0_int,
	                      "q0" = q0,
                        "beta" = beta,
	                      "theta" = theta,
                        "lma" = lma,
                        "given_lma" = given_lma,
                        "f" = f,
                        "patm" = patm,
                        "par" = par,
                        "vpd" = vpd,
                        "ca" = ca,
                        "oa" = oa,
                        "km" = km,
                        "gammastar" = gammastar,
                        "chi" = chi,
                        "given_chi" = given_chi,
                        "ci" = ci,
                        "mc" = mc,
                        "m" = m,
                        "omega" = omega,
                        "omega_star" = omega_star,
                        "Al" = Al,
                        "kp" = kp,
                        "kr" = kr,
                        "ko" = ko,
                        "Ap" = Ap,
                        "leakage" = leakage,
                        "cbs" = cbs,
                        "chi_bs" = chi_bs,
                        "obs" = obs,
                        "Ac" = Ac,
                        "rd" = rd,
                        "Anet" = Anet,
                        "jmax" = jmax,
                        "vpmax" = vpmax,
	                      "vcmax" = vcmax,
                        "jmax25" = jmax25,
                        "vpmax25" = vpmax25,
                        "vcmax25" = vcmax25,
                        "rd25" = rd25,
	                      "nrubisco" = nrubisco, 
	                      "nbioe" = nbioe, 
	                      "npep" = npep, 
	                      "nstructure" = nstructure,
	                      "nall" = nall,
	                      "nphoto" = nphoto, 
	                      "nrubisco_frac" = nrubisco_frac, 
	                      "nphoto_frac" = nphoto_frac,
                        "gsw"=gsw,
                        "gsc"=gsc,
                        "narea"=narea,
                        "nmass"=nmass,
                        "parea"=parea,
                        "wue"=wue,
                        "nue"=nue,
                        "pue"=pue)

	return(results)
}
