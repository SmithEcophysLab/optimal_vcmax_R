# test optimal vcmax predicitons

source('calc_optimal_vcmax.R')

par_seq = calc_optimal_vcmax(paro = seq(100, 1000, 100))
t_seq = calc_optimal_vcmax(tg_c = seq(10, 30, 5))
vpd_seq = calc_optimal_vcmax(vpdo = seq(0.5, 4, 0.5))
z_seq = calc_optimal_vcmax(z = seq(0, 1000, 100))
ca_seq = calc_optimal_vcmax(cao = seq(250, 1000, 50))

