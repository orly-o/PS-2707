#####################################
## Class work for Application Day:
## unsupervised
#####################################

#####################################
library(tidyverse)
library(MCMCpack)
library(scales)
library(emIRT)
library(dplyr)
#####################################

# remotes::install_github("xmarquez/democracyData")

setwd('C:/Users/orlyo/OneDrive/Desktop/Grad School/Spring 2022/1. PS 2707')
# load('democracy1946.2000.rda')
load('democracy1946.2008.rda')

# set NAs to zero
democracy[is.na(democracy)] <- 0

# rescale different democracy scores to 0/1
democracy$arat_scale = rescale(democracy$arat)
democracy$blm_scale = rescale(democracy$blm)
democracy$bollen_scale = rescale(democracy$bollen)
democracy$freedomhouse_scale = rescale(democracy$freedomhouse)
democracy$hadenius_scale = rescale(democracy$hadenius)
democracy$mainwaring_scale = rescale(democracy$mainwaring)
democracy$munck_scale = rescale(democracy$munck)
democracy$pacl_scale = rescale(democracy$pacl)
democracy$polity_scale = rescale(democracy$polity)
democracy$polyarchy_scale = rescale(democracy$polyarchy)
democracy$prc_scale = rescale(democracy$prc)
democracy$vanhanen_scale = rescale(democracy$vanhanen)

# subset to include scaled variables
democracy_sub = democracy[,c('cowcode', 'country', 'year', 'arat_scale', 'blm_scale', 'bollen_scale', 'freedomhouse_scale', 
                             'hadenius_scale', 'mainwaring_scale', 'munck_scale', 'pacl_scale', 'polity_scale', 
                             'polyarchy_scale', 'prc_scale', 'vanhanen_scale')]

# prep for emIRT
# combine
democracy_sub$countryyear = paste(democracy_sub$country, democracy_sub$year, sep = ' ')
rownames(democracy_sub) = democracy_sub$countryyear
democracy_sub2 = subset(democracy_sub, select = -c(cowcode, country, year, countryyear))

# make binary
democracy_sub2[democracy_sub2 <= 0.5 & democracy_sub2 > 0] = -1
democracy_sub2[democracy_sub2 > 0.5] = 1

# matrix
democracy_mat = as.matrix(democracy_sub2)

# model
# starting
start_dec <- list(
  alpha = matrix(rnorm(ncol(democracy_mat))),
  beta = matrix(rnorm(ncol(democracy_mat))),
  x = matrix(rnorm(nrow(democracy_mat)))
)

# setup prior
prior_dec <- list(
  x = list(mu = matrix(0), sigma = matrix(1)),
  beta = list(mu = matrix(c(0,0)), sigma = 25 * diag(2))
)

# run the model
dec_bin_emIRT <- binIRT(.rc = list(votes = democracy_mat),
                        .starts = start_dec,
                        .priors = prior_dec,
                        .control = list(verbose = FALSE, thresh = 1e-6)
)


###

# means from binIRT model come out in the same order - just make it a new column
democracy_sub$ideal = dec_bin_emIRT$means$x

# plot country and ideal point
ggplot(democracy_sub) +
  geom_line(aes(y = country, x = ideal))

# hard to see!

# see a few countries separately
country = 'United States'
dem = democracy_sub[democracy_sub$country == country,]
ggplot(dem_USA) +
  geom_line(aes(x = year, y = ideal)) +
  ggtitle(country)

country = 'Spain'
dem = democracy_sub[democracy_sub$country == country,]
ggplot(dem) +
  geom_line(aes(x = year, y = ideal)) +
  ggtitle(country)

country = 'Israel'
dem = democracy_sub[democracy_sub$country == country,]
ggplot(dem) +
  geom_line(aes(x = year, y = ideal)) +
  ggtitle(country)

country = 'Fiji'
dem = democracy_sub[democracy_sub$country == country,]
ggplot(dem) +
  geom_line(aes(x = year, y = ideal)) +
  ggtitle(country)

# next: plot ideal  points with uncertainty interval for countries




