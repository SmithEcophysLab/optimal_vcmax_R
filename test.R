# test LMA predictions on C3 and C4 NutNet data

library(dplyr)
library(ggplot2)

nutnet <- read.csv('NutNet_traits.csv')

## load model
source('calc_optimal_vcmax.R')
source('calc_vcmax_C3.R')
source('calc_vcmax_C4.R')
sourceDirectory('functions', modifiedOnly = FALSE)

nutnet_C3 <- dplyr::filter(nutnet, photosynthetic_pathway == "C3")
nutnet_C4 <- dplyr::filter(nutnet, photosynthetic_pathway == "C4")

model_output_C3 <- calc_optimal_vcmax(pathway = "C3",
                                     tg_c = nutnet_C3$tmp, 
                                     z = nutnet_C3$z, 
                                     vpdo = nutnet_C3$vpd,
                                     paro = nutnet_C3$par,
                                     cao = 391)

model_output_C4 <- calc_optimal_vcmax(pathway = "C4",
                                      tg_c = nutnet_C4$tmp, 
                                      z = nutnet_C4$z, 
                                      vpdo = nutnet_C4$vpd,
                                      paro = nutnet_C4$par,
                                      cao = 391)

nutnet_C3$lma_pred <- model_output_C3$lma
nutnet_C4$lma_pred <- model_output_C4$lma

ggplot(nutnet_C3, aes(lma_pred, lma)) + geom_point()
ggplot(nutnet_C4, aes(lma_pred, lma)) + geom_point()

