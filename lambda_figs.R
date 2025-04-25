# lambda_figs.R
## make some hypothesis figures for the NutNet lambda storyboard/manuscript

## load libraries
library(ggplot2)

## source model
source('calc_optimal_vcmax.R')

## run sequences
par_seq = calc_optimal_vcmax(paro = seq(100, 1000, 100))
t_seq = calc_optimal_vcmax(tg_c = seq(10, 30, 5))
vpd_seq = calc_optimal_vcmax(vpdo = seq(0.5, 4, 0.5))
z_seq = calc_optimal_vcmax(z = seq(0, 1000, 100))
ca_seq = calc_optimal_vcmax(cao = seq(250, 1000, 50))
beta_seq = calc_optimal_vcmax(beta = seq(20, 500, 20))

par_seq_c4 = calc_optimal_vcmax(paro = seq(100, 1000, 100), pathway = "C4")
t_seq_c4 = calc_optimal_vcmax(tg_c = seq(10, 30, 5), pathway = "C4")
vpd_seq_c4 = calc_optimal_vcmax(vpdo = seq(0.5, 4, 0.5), pathway = "C4")
z_seq_c4 = calc_optimal_vcmax(z = seq(0, 1000, 100), pathway = "C4")
ca_seq_c4 = calc_optimal_vcmax(cao = seq(250, 1000, 50), pathway = "C4")
beta_seq_c4 = calc_optimal_vcmax(beta = seq(20, 500, 20), pathway = "C4")

par_seq_c4$c4_adv <- ((par_seq_c4$Al - par_seq$Al) / par_seq$Al) * 100
t_seq_c4$c4_adv <- ((t_seq_c4$Al - t_seq$Al) / t_seq$Al) * 100
ca_seq_c4$c4_adv <- ((ca_seq_c4$Al - ca_seq$Al) / ca_seq$Al) * 100
beta_seq_c4$c4_adv <- ((beta_seq_c4$Al - beta_seq$Al) / beta_seq$Al) * 100

## make figures
narea_soiln_fig <- ggplot(data = beta_seq, aes(x = 1/beta, y = nall)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', linewidth = 3) +
  #geom_line(data = beta_seq_c4, color = 'brown', linewidth = 3) +
  ylab(expression(italic('N')['area'] * ' (gN m' ^ '-2' * ')')) +
  xlab(expression('Soil nutrient availability (1/β)')) +
  ylim(c(2.49, 2.54)) +
  xlim(c(0, 0.05))

narea_vpd_fig <- ggplot(data = vpd_seq, aes(x = vpd, y = nall)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', linewidth = 3) +
  ylab(expression(italic('N')['area'] * ' (gN m' ^ '-2' * ')')) +
  xlab(expression('VPD (kPa)')) +
  ylim(c(2.4, 2.8)) +
  xlim(c(0, 4))


## save figures

jpeg(filename = "lambda_figs/narea_soiln_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(narea_soiln_fig)
dev.off()

jpeg(filename = "lambda_figs/narea_vpd_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(narea_vpd_fig)
dev.off()