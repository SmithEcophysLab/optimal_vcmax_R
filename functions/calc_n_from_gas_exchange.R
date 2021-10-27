# leaf Nrubisco from vcmax25
## calculate N from vcmax
### from Ning Dong and Harrison et al. (2009): Nrubisco (g per m2) = ((vcmax25 * Mr * Nr) / (kcat * nr)) * Mn
### vcmax25 = vcmax at 25C (umol m-2 s-1)
### Nr = N concentration of Rubisco = 0.0114 mol N per g Rubisco
### kcat = catalytic turnover at 25C = 3.5 mol CO2 per mol Rubisco site per second (3500000 µmol CO2 per mol Rubisco site per second)
### Mr = molecular mass of Rubisco = 0.55 g Rubisco per umol Rubisco (550000 g Rubisco per mol Rubisco)
### nr = number of catalytic sites per mole of Rubisco = 8 mol Rubisco site per mol Rubisco
### Mn = molecular mass of N = 14 gN per mol N
# https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-3040.2008.01918.x
fvcmax25_nrubisco = function(vcmax25){
  
  Mr = 550000 # g Rubisco per mol Rubisco
  Nr = 0.0114 # mol N per g Rubisco
  Mn = 14 # gN per mol N
  nr =  8 # mol Rubisco site per mol Rubisco
  kcat = 3500000 # µmol CO2 per mol Rubisco site per second
  
  nrubisco = vcmax25 * Mr * Nr * Mn * (1/nr) * (1/kcat)
  # nrubisco = vcmax25 * (1/3500000) * (8) * 550000 * (1/0.0114) * 14
  nrubisco # leaf n in rubisco gN/m2
  
}

# leaf Nbioenergetics from jmax25
## calculate N from jmax
### from Niinemets and Tenhunen (1997): Nbioe (g per m2) = (jmax25 * Ncyt) / jmc
### jmax25 = jmax at 25C (umol m-2 s-1)
### jmc = activity of electron transport at 25C = 156 µmol e- (µmol cyt f)-1 * s-1
### Ncyt = investment in bioenergetics = 0.124 gN (µmol cyt f)-1
# https://onlinelibrary.wiley.com/doi/epdf/10.1046/j.1365-3040.1997.d01-133.x
fjmax25_nbioe = function(jmax25){
  
  nbioe = jmax25 * 0.124 *(1/156)
  nbioe # leaf n in bioenergetics gN/m2
  
}

# structural N
### lma conversion from Dong et al. (2017) Biogeosciences
flma_nstructure = function(lma){
  
  nstructure = (10^-2.67) * (lma ^ 0.99)
  nstructure
  
}

fvpmax25_npep = function(vpmax25){ # assuming smaller size, but similar catalytic functions to rubisco
  
  Mr = 410000 # https://www.jbc.org/content/278/14/11867
  Nr = 0.0114 # mol N per g Rubisco (similar for PEP? https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1054259/pdf/plntphys00619-0049.pdf)
  Mn = 14 # gN per mol N
  nr =  2 # mol Rubisco site per mol Rubisco
  kcat = 5440000 # µmol CO2 per mol Rubisco site per second # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4634053/
  
  npep = vpmax25 * Mr * Nr * Mn * (1/nr) * (1/kcat)
  # nrubisco = vcmax25 * (1/3500000) * (8) * 550000 * (1/0.0114) * 14
  npep # leaf n in rubisco gN/m2
  
}

