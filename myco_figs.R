# make some figures for lalasia's myco paper

## load libraries
library(ggplot2)
library(RColorBrewer)
library(grid)

## load model
source('calc_optimal_vcmax.R')

## make model for gs-N tradeoff
beta_seq_par400 = calc_optimal_vcmax(beta = seq(9, 2000, 1), par = 400)
beta_seq_par500 = calc_optimal_vcmax(beta = seq(8, 1800, 1), par = 500)
beta_seq_par600 = calc_optimal_vcmax(beta = seq(7, 1600, 1), par = 600)
beta_seq_par700 = calc_optimal_vcmax(beta = seq(6, 1400, 1), par = 700)
beta_seq_par800 = calc_optimal_vcmax(beta = seq(5, 1200, 1), par = 800)
beta_seq_par900 = calc_optimal_vcmax(beta = seq(4, 1000, 1), par = 900)
beta_seq_par1000 = calc_optimal_vcmax(beta = seq(3, 800, 1), par = 1000)
beta_seq_par1100 = calc_optimal_vcmax(beta = seq(2, 575, 1), par = 1100)
beta_seq_par1200 = calc_optimal_vcmax(beta = seq(5, 900, 1), par = 1200)
beta_seq_par1300 = calc_optimal_vcmax(beta = seq(5, 900, 1), par = 1300)
beta_seq_par1400 = calc_optimal_vcmax(beta = seq(5, 900, 1), par = 1400)
beta_seq_par1500 = calc_optimal_vcmax(beta = seq(5, 900, 1), par = 1500)
beta_seq_par1600 = calc_optimal_vcmax(beta = seq(5, 900, 1), par = 1600)
beta_seq_par1700 = calc_optimal_vcmax(beta = seq(5, 900, 1), par = 1700)
beta_seq_par1800 = calc_optimal_vcmax(beta = seq(5, 900, 1), par = 1800)

## make color palette
fig_colors <- brewer.pal(9, "Greens")

## make plots
vcmax_gs_tradeoff_plot <- ggplot(aes(x = gsw, y = vcmax), data = beta_seq_par1200) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_blank()) +
  geom_line(color = fig_colors[9]) +
  geom_line(data = beta_seq_par1300, color = fig_colors[8]) +
  geom_line(data = beta_seq_par1400, color = fig_colors[7]) +
  geom_line(data = beta_seq_par1500, color = fig_colors[6]) +
  geom_line(data = beta_seq_par1600, color = fig_colors[5]) +
  geom_line(data = beta_seq_par1700, color = fig_colors[4]) +
  geom_line(data = beta_seq_par1800, color = fig_colors[3]) +
  ylab(expression(italic('V')[cmax])) +
  xlab(expression(italic('g')[sw]))

jpeg(filename = "figs/vcmax_gs_tradeoff_plot.jpeg", 
    width = 4, height = 4, units = 'in', res = 300)
grid.newpage()
grid.draw(vcmax_gs_tradeoff_plot)
dev.off()
