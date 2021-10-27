# test LMA predictions on C3 and C4 NutNet data

library(ggplot2)

nutnet <- read.csv('NutNet_traits.csv')

## load model
source('../calc_optimal_vcmax.R')
sourceDirectory('../functions', modifiedOnly = FALSE)

model_output <- calc_optimal_vcmax(tg_c = nutnet$tmp, 
                                   z = nutnet$z, 
                                   vpdo = nutnet$vpd,
                                   paro = nutnet$par,
                                   cao = 391)

nutnet$lma_pred <- model_output$lma

ggplot(nutnet, aes(lma_pred, lma)) + geom_point(aes(color = photosynthetic_pathway))

