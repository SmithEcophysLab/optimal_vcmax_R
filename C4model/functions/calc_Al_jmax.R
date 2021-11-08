calc_Al_jmax <- function(jmax, tg_c, par, q0, theta, m){
  # q0 <- q025 * phi_ftemp(tg_c)
  
  omega <- jmax / (q0 * par)
  omega_star <- (1 + (omega) - sqrt((1 + (omega))^2 - (4 * theta * omega)))
  Al <- q0 * par * m * omega_star / (8 * theta)
  
  return(Al)
}