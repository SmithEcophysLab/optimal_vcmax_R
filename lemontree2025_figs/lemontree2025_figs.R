# test optimal vcmax predicitons

## load libraries
library(ggplot2)

## source and test functions
source('../calc_optimal_vcmax.R')
sourceDirectory('../functions')

par_seq = calc_optimal_vcmax(paro = seq(100, 1000, 10))
t_seq = calc_optimal_vcmax(tg_c = seq(10, 30, 1))
vpd_seq = calc_optimal_vcmax(vpdo = seq(0.5, 4, 0.1))
z_seq = calc_optimal_vcmax(z = seq(0, 1000, 100))
ca_seq = calc_optimal_vcmax(cao = seq(250, 1000, 10))
beta_seq = calc_optimal_vcmax(beta = seq(10, 500, 10))

## make lemontree2025 figs
par_nphoto_fig <- ggplot(data = par_seq, aes(x = paro, y = nphoto)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(4), colour = 'black'),
        axis.title.x = element_text(size = rel(4), colour = 'black'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression(italic('N')['photo'])) +
  xlab(expression('PAR')) +
  ylim(c(0, 0.7)) +
  xlim(c(0, 1000))

jpeg(filename = "figures/par_nphoto_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(par_nphoto_fig)
dev.off()

t_nphoto_fig <- ggplot(data = t_seq, aes(x = tg_c, y = nphoto)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(4), colour = 'black'),
        axis.title.x = element_text(size = rel(4), colour = 'black'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression(italic('N')['photo'])) +
  xlab(expression('Temperature')) +
  ylim(c(0.3, 0.7)) +
  xlim(c(20, 30))

jpeg(filename = "figures/t_nphoto_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(t_nphoto_fig)
dev.off()

beta_nphoto_fig <- ggplot(data = beta_seq, aes(x = beta, y = nphoto)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(4), colour = 'black'),
        axis.title.x = element_text(size = rel(4), colour = 'black'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression(italic('N')['photo'])) +
  xlab(expression('Water avail.')) +
  ylim(c(0.45, 0.55)) +
  xlim(c(10, 500))

jpeg(filename = "figures/beta_nphoto_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(beta_nphoto_fig)
dev.off()

par_vcmax25_fig <- ggplot(data = par_seq, aes(x = paro, y = vcmax25)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(4), colour = 'black'),
        axis.title.x = element_text(size = rel(4), colour = 'black'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression(italic('V')['cmax25'])) +
  xlab(expression('PAR')) +
  ylim(c(0, 125)) +
  xlim(c(0, 1000))

jpeg(filename = "figures/par_vcmax25_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(par_vcmax25_fig)
dev.off()

beta_vcmax25_fig <- ggplot(data = beta_seq, aes(x = 1/beta, y = vcmax25)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(4), colour = 'black'),
        axis.title.x = element_text(size = rel(4), colour = 'black'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression(italic('V')['cmax25'])) +
  xlab(expression('N avail.')) +
  ylim(c(0, 125)) +
  xlim(c(0, 0.1))

jpeg(filename = "figures/beta_vcmax25_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(beta_vcmax25_fig)
dev.off()

ca_vcmax25_fig <- ggplot(data = ca_seq, aes(x = ca, y = vcmax25)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(4), colour = 'black'),
        axis.title.x = element_text(size = rel(4), colour = 'black'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(color = 'darkblue', lwd = 3) +
  ylab(expression(italic('V')['cmax25'])) +
  xlab(expression('CO'[2])) +
  ylim(c(70, 125)) +
  xlim(c(20, 110))

jpeg(filename = "figures/ca_vcmax25_fig.jpeg", width = 5.5, height = 5.5, units = 'in', res = 600)
plot(ca_vcmax25_fig)
dev.off()
