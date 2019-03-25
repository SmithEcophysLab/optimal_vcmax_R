# code to examine different ways to calculate vcmax25

# setwd("/Users/nicksmith/Documents/Git/vcmax25")

## load libraries
library(R.utils)

## load optimal vcmax script (Smith et al., Ecology Letters)
sourceDirectory('optimal_vcmax/functions')
source('optimal_vcmax/calc_optimal_vcmax.R')

## data for temperature and light responses
temperature = calc_optimal_vcmax(tg_c = seq(5, 35, 5))
head(temperature)

light = calc_optimal_vcmax(paro = seq(100, 1500, 50))
head(light)

## some plots
plot((temperature$vcmax_25 / temperature$vcmax_25[5]) ~ temperature$tg_c, type = 'l', lwd = 4, col = 'blue', ylim = c(0, 4))
lines((temperature$vcmax_25_alt / temperature$vcmax_25_alt[5]) ~ temperature$tg_c, type = 'l', lwd = 4, col = 'red')
lines((temperature$vcmax_star / temperature$vcmax_star[5]) ~ temperature$tg_c, type = 'l', lwd = 4, col = 'black')
lines((temperature$vcmax_prime / temperature$vcmax_prime[5]) ~ temperature$tg_c, type = 'l', lwd = 4, col = 'grey')

plot((light$vcmax_25 / light$vcmax_25[5]) ~ light$par, type = 'l', lwd = 4, col = 'blue', ylim = c(0, 5))
lines((light$vcmax_25_alt / light$vcmax_25_alt[5]) ~ light$par, type = 'l', lwd = 4, col = 'red')
