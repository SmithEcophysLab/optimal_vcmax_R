# test optimal vcmax predicitons

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
par_fig <- ggplot(data = par_seq, aes(x = paro, y = vcmax)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression(italic('V')['cmax'] * ' (µmol m' ^ '-2 ' * 's' ^ '-1' * ')')) +
  xlab(expression('PAR' * ' (µmol m' ^ '-2 ' * 's' ^ '-1' * ')')) +
  ylim(c(0, 150)) +
  xlim(c(0, 1000))

t_fig <- ggplot(data = t_seq, aes(x = tg_c, y = vcmax)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression(italic('V')['cmax'] * ' (µmol m' ^ '-2 ' * 's' ^ '-1' * ')')) +
  xlab(expression('Temperature (°C)')) +
  ylim(c(0, 150)) +
  xlim(c(10, 30))

ca_fig <- ggplot(data = ca_seq, aes(x = cao, y = vcmax)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression(italic('V')['cmax'] * ' (µmol m' ^ '-2 ' * 's' ^ '-1' * ')')) +
  xlab(expression('CO'['2'] * ' (µmol mol' ^ '-1' * ')')) +
  ylim(c(50, 150)) +
  xlim(c(200, 1000))

beta_fig <- ggplot(data = beta_seq, aes(x = beta, y = vcmax)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression(italic('V')['cmax'] * ' (µmol m' ^ '-2 ' * 's' ^ '-1' * ')')) +
  xlab(expression('Soil N limitation')) +
  ylim(c(80, 120)) +
  xlim(c(0, 500))

t_fig_c4_adv <- ggplot(data = t_seq_c4, aes(x = tg_c, y = c4_adv)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression('C'['4'] * ' advantage (%)')) +
  xlab(expression('Temperature (°C)')) +
  ylim(c(20, 40)) +
  xlim(c(10, 30))

ca_fig_c4_adv <- ggplot(data = ca_seq_c4, aes(x = cao, y = c4_adv)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression('C'['4'] * ' advantage (%)')) +
  xlab(expression('CO'['2'] * ' (µmol mol' ^ '-1' * ')')) +
  ylim(c(20, 50)) +
  xlim(c(200, 1000))

## save figures

jpeg(filename = "talk_figs/par_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(par_fig)
dev.off()

jpeg(filename = "talk_figs/t_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(t_fig)
dev.off()

jpeg(filename = "talk_figs/ca_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(ca_fig)
dev.off()

jpeg(filename = "talk_figs/beta_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(beta_fig)
dev.off()

jpeg(filename = "talk_figs/t_fig_c4_adv.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(t_fig_c4_adv)
dev.off()

jpeg(filename = "talk_figs/ca_fig_c4_adv.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(ca_fig_c4_adv)
dev.off()
