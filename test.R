# test LMA predictions on C3 and C4 NutNet data

library(dplyr)
library(ggplot2)

nutnet <- read.csv('NutNet_traits.csv')

## load model
source('calc_optimal_vcmax.R')
sourceDirectory('functions', modifiedOnly = FALSE)

nutnet$chi <- rep(NA)
nutnet_C3 <- dplyr::filter(nutnet, photosynthetic_pathway == "C3")
nutnet_C4 <- dplyr::filter(nutnet, photosynthetic_pathway == "C4")

nutnet_C3$lma_test <- rep(NA) 

model_output_C3 <- calc_optimal_vcmax(pathway = "C3",
                                     tg_c = nutnet_C3$tmp, 
                                     z = nutnet_C3$z, 
                                     vpdo = nutnet_C3$vpd,
                                     paro = nutnet_C3$par,
                                     cao = 391, lma = nutnet_C3$lma_test)

model_output_C4 <- calc_optimal_vcmax(pathway = "C4",
                                      tg_c = nutnet_C4$tmp, 
                                      z = nutnet_C4$z, 
                                      vpdo = nutnet_C4$vpd,
                                      paro = nutnet_C4$par,
                                      cao = 391, lma = nutnet_C4$lma)

nutnet_C3$model_lma <- model_output_C3$lma
ggplot(nutnet_C3, aes(x = model_lma, y = lma)) + geom_point()
