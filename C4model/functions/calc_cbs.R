# Bundle Sheath Conductance to Plug into whole model

calc_cbs <- function(cm, leakage) { # "leakiness" is the same as phi in Von Caemmerer (eqn. 4.5)

  gbs <- 0.003 # 3 mmol m^-2 s^-1 * 1000 to convert to micromol

  Cbs <- (leakage + (gbs * cm)) /gbs # rearrangement of eqn. 4.4 in Von Caemmerer
  
  return(Cbs)
}