# make some nitrogen figures

library(ggplot2)

source('calc_optimal_vcmax.R')

par_seq = calc_optimal_vcmax(paro = seq(100, 1000, 100))
t_seq = calc_optimal_vcmax(tg_c = seq(15, 30, 1))
vpd_seq = calc_optimal_vcmax(vpdo = seq(0.5, 4, 0.5))
z_seq = calc_optimal_vcmax(z = seq(0, 1000, 100))
ca_seq = calc_optimal_vcmax(cao = seq(250, 1000, 50))
beta_seq = calc_optimal_vcmax(beta = seq(20, 500, 20))

par_fig <- ggplot(data = par_seq, aes(x = paro, y = nphoto)) +
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
  ylab(expression(italic('N')['photo'] * ' (g m' ^ '-2' * ')')) +
  xlab(expression('PAR' * ' (µmol m' ^ '-2 ' * 's' ^ '-1' * ')')) +
  ylim(c(0, 0.8)) +
  xlim(c(0, 1000))

t_fig <- ggplot(data = t_seq, aes(x = tg_c, y = nphoto)) +
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
  ylab(expression(italic('N')['photo'] * ' (g m' ^ '-2' * ')')) +
  xlab(expression('Temperature (°C)')) +
  ylim(c(0.3, 0.7)) +
  xlim(c(15, 30))

vpd_fig <- ggplot(data = vpd_seq, aes(x = vpdo, y = nphoto)) +
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
  ylab(expression(italic('N')['photo'] * ' (g m' ^ '-2' * ')')) +
  xlab(expression('VPD (kPa)')) +
  ylim(c(0.45, 0.55)) +
  xlim(c(0.5, 4))

ca_fig <- ggplot(data = ca_seq, aes(x = cao, y = nphoto)) +
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
  ylab(expression(italic('N')['photo'] * ' (g m' ^ '-2' * ')')) +
  xlab(expression('CO'['2'] * ' (µmol mol' ^ '-1' * ')')) +
  ylim(c(0.4, 0.55)) +
  xlim(c(200, 1000))

beta_fig <- ggplot(data = beta_seq, aes(x = beta, y = nphoto)) +
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
  ylab(expression(italic('N')['photo'] * ' (g m' ^ '-2' * ')')) +
  xlab(expression('Soil N limitation')) +
  ylim(c(0.45, 0.55)) +
  xlim(c(0, 500))

beta_fig_gsw <- ggplot(data = beta_seq, aes(x = beta, y = gsw)) +
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
  ylab(expression(italic('g')['sw'] * ' (mmol m' ^ '-2 ' * 's' ^'-1' * ')')) +
  xlab(expression('Soil N limitation')) +
  ylim(c(0, 0.8)) +
  xlim(c(0, 500))

## save figures

jpeg(filename = "nitrogen_figures/par_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(par_fig)
dev.off()

jpeg(filename = "nitrogen_figures/t_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(t_fig)
dev.off()

jpeg(filename = "nitrogen_figures/ca_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(ca_fig)
dev.off()

jpeg(filename = "nitrogen_figures/vpd_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(vpd_fig)
dev.off()

jpeg(filename = "nitrogen_figures/beta_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(beta_fig)
dev.off()

jpeg(filename = "nitrogen_figures/beta_fig_gsw.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(beta_fig_gsw)
dev.off()
