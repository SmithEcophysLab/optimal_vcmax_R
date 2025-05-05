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

## c3c4 state figures
library(raster)
library(RColorBrewer)
library(maps)
library(mapdata)
library(gridBase)
library(mapproj)
library(grDevices)
library(geodata)

tmp_globe_4model = read.csv('/Users/nicksmith/Documents/Research/Colimitation/Spatial_Maps/cru_tmp_climExtract_growingseason_globe.csv')
par_globe_4model = read.csv('/Users/nicksmith/Documents/Research/Colimitation/Spatial_Maps/cru_par_climExtract_growingseason_globe.csv')
vpd_globe_4model = read.csv('/Users/nicksmith/Documents/Research/Colimitation/Spatial_Maps/cru_vpd_climExtract_growingseason_globe.csv')
z_globe_4model =  read.csv('/Users/nicksmith/Documents/Research/Colimitation/Spatial_Maps/z_globe.csv')

c3_globe_pres = calc_optimal_vcmax(tg_c = tmp_globe_4model$tmp, z = z_globe_4model$z, vpdo = vpd_globe_4model$vpd, 
                       cao = 400, paro = par_globe_4model$par, pathway = "C3")
c4_globe_pres = calc_optimal_vcmax(tg_c = tmp_globe_4model$tmp, z = z_globe_4model$z, vpdo = vpd_globe_4model$vpd, 
                       cao = 400, paro = par_globe_4model$par, pathway = "C4")
c3_globe_fut = calc_optimal_vcmax(tg_c = tmp_globe_4model$tmp + 4, z = z_globe_4model$z, vpdo = vpd_globe_4model$vpd, 
                      cao = 1000, paro = par_globe_4model$par, pathway = "C3")
c4_globe_fut = calc_optimal_vcmax(tg_c = tmp_globe_4model$tmp + 4, z = z_globe_4model$z, vpdo = vpd_globe_4model$vpd, 
                      cao = 1000, paro = par_globe_4model$par, pathway = "C4")

c3_globe_pres_A = cbind(tmp_globe_4model$lon, tmp_globe_4model$lat, c3_globe_pres$Ac)
c4_globe_pres_A = cbind(tmp_globe_4model$lon, tmp_globe_4model$lat, c4_globe_pres$Ac)
c3_globe_fut_A = cbind(tmp_globe_4model$lon, tmp_globe_4model$lat, c3_globe_fut$Ac)
c4_globe_fut_A = cbind(tmp_globe_4model$lon, tmp_globe_4model$lat, c4_globe_fut$Ac)

globe_pres_A_adv = cbind(tmp_globe_4model$lon, tmp_globe_4model$lat, ((c4_globe_pres$Ac - c3_globe_pres$Ac)/c3_globe_pres$Ac) * 100)
globe_fut_A_adv = cbind(tmp_globe_4model$lon, tmp_globe_4model$lat, ((c4_globe_fut$Ac - c3_globe_fut$Ac)/c3_globe_fut$Ac) * 100)

## plot
us <- gadm(country="USA", level=1, path = '~/Documents/Data/geodata')
neb_usa = c("Nebraska", "Kansas", "Oklahoma", "Texas")
neb = us[match(toupper(neb_usa),toupper(us$NAME_1)),]

globe_pres_A_adv_ras = rast(globe_pres_A_adv, type ='xyz')
globe_pres_A_adv_ras[globe_pres_A_adv_ras > 30] <- 30
globe_pres_A_adv_ras_crop <- crop(globe_pres_A_adv_ras, neb)
mask_pres = mask(globe_pres_A_adv_ras_crop, neb)

globe_fut_A_adv_ras = rast(globe_fut_A_adv, type ='xyz')
globe_fut_A_adv_ras[globe_fut_A_adv_ras > 30] <- 30
globe_fut_A_adv_ras_crop <- crop(globe_fut_A_adv_ras, neb)
mask_fut = mask(globe_fut_A_adv_ras_crop, neb)

mean(as.matrix(mask_pres), na.rm = T)
mean(as.matrix(mask_fut), na.rm = T)

pale = colorRampPalette(c(brewer.pal(9,'Reds')))
cols = pale(21)
#arg = list(at = seq(0, 30, 2), labels = seq(0, 30, 2))

par(mfrow = c(2, 1), mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 3))
plot(neb, axes = F)
plot(mask_pres, add = T, col = cols, 
     legend=T,
     type = 'continuous', range = c(0, 30),
     axes = F, smooth = T)
plot(neb, add = T, axes = F)
title('Present day')

plot(neb, axes = F)
plot(mask_fut,add = T, col = cols, 
     legend=T,
     type = 'continuous', range = c(0, 30),
     axes = F, smooth = T)
plot(neb, add = T)
title('Future')


