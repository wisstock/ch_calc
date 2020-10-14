# Flyorescent spectra anslysis
# Copyright © 2020 Borys Olifirov

require(ggplot2)
require(gridExtra)
require(dplyr)
require(magrittr)
require(wavelength2colour)

setwd('/home/astria/Bio/Note/diff/fluo')

fluo.1 <- read.csv('fluo_4.csv')
fluo.2 <- read.csv('mTFP1.csv')
pass.band.1 <- c(475, 500) 
pass.band.2 <- c(540, 600)


fluo.1[is.na(fluo.1)] <- 0
fluo.1 <- subset(fluo.1, w >= 300)
if (max(fluo.1$ex) < 90) {
  fluo.1$ex <- fluo.1$ex * 100
  fluo.1$em <- fluo.1$em * 100
}

fluo.2[is.na(fluo.2)] <- 0
fluo.2 <- subset(fluo.2, w >= 300)
if (max(fluo.2$ex) < 90) {
  fluo.2$ex <- fluo.2$ex * 100
  fluo.2$em <- fluo.2$em * 100
}

fluo.1.ex.col <- wavelength2hex(fluo.1$w[fluo.1$ex == max(fluo.1$ex)])
fluo.1.em.col <- wavelength2hex(530) # fluo.1$w[fluo.1$em == max(fluo.1$em)])

fluo.2.ex.col <- wavelength2hex(fluo.2$w[fluo.2$ex == max(fluo.2$ex)])
fluo.2.em.col <- wavelength2hex(fluo.2$w[fluo.2$em == max(fluo.2$em)])


ggplot() +
  geom_ribbon(data = fluo.1, mapping = aes(x = w,
                                        ymin = 0,
                                        ymax = ex),
              colour = fluo.1.ex.col,
              fill = fluo.1.ex.col,
              size = 0.2,
              alpha = 0.25) +
  geom_ribbon(data = fluo.1, mapping = aes(x = w,
                                        ymin = 0,
                                        ymax = em),
              colour = fluo.1.em.col,
              fill = fluo.1.em.col,
              size = 1,
              alpha = 0.25) +
  geom_ribbon(data = fluo.2, mapping = aes(x = w,
                                      ymin = 0,
                                      ymax = ex),
              colour = fluo.2.ex.col,
              fill = fluo.2.ex.col,
              size = 0.2,
              alpha = 0.25) +
  geom_ribbon(data = fluo.2, mapping = aes(x = w,
                                       ymin = 0,
                                       ymax = em),
              colour = fluo.2.em.col,  # '#808080',
              fill = fluo.2.em.col,
              size = 1,
              alpha = 0.25) +
  geom_rect(aes(xmin = pass.band.1[1],
                xmax = pass.band.1[2],
                ymin = 0, ymax = 110),
            fill = 'red', 
            alpha = 0.3) +
  geom_rect(aes(xmin = pass.band.2[1],
                xmax = pass.band.2[2],
                ymin = 0, ymax = 110),
            fill = 'red', 
            alpha = 0.3) +
  scale_x_continuous(limits = c(300, 700),
                     breaks = seq(300, 700, 50)) +
  scale_y_continuous(limits = c(0, 110),
                     breaks = seq(0, 100, 20)) +
  labs(y ='Intensity (%)', x = 'Wavelength (nm)') +
  theme_minimal(base_size = 18,
                base_family = 'oswald')

sum(fluo.2$em[fluo.2$w >= pass.band.1[1] & fluo.2$w <= pass.band.1[2]]) / sum(fluo.2$em)

