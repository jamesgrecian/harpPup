#######################################
### Harp seal morphometric analysis ###
#######################################

# load libraries
require(tidyverse)

# load data
dat <- readRDS("data/mophometrics.rds")

# what is the range of weights of tagged animals? Need this for opening of methods section...
dat %>% pull(weight_kg) %>% median
dat %>% pull(weight_kg) %>% range

# do the two populations differ in length?
m0 <- lm(length_cm ~ 1, data = dat)
m1 <- lm(length_cm ~ population, data = dat)
anova(m0, m1)
summary(m1)

# do the two populations differ in girth?
m0 <- lm(girth_cm ~ 1, data = dat)
m1 <- lm(girth_cm ~ population, data = dat)
anova(m0, m1)
summary(m1)

# do the two populationa differ in weight
m0 <- lm(weight_kg ~ 1, data = dat)
m1 <- lm(weight_kg ~ population, dat = dat)
anova(m0, m1)
summary(m1)
# Canadian seals are heavier by ~ 2.2 kg

# How does body condition differ between the two populations?
# Hammill 1995
# volume index = length * girth^2
# take residuals of linear model
vi_m <- lm(weight_kg ~ length_cm * I(girth_cm^2), dat = dat)
summary(vi_m) # adj R squared = 0.487
dat$bmi = residuals(vi_m)

# plot to check
ggplot() +
  geom_histogram(aes(x = bmi), bins = 8, data = dat) +
  facet_wrap(~ population, ncol = 1)

# model to estimate population difference in bmi
m1 <- lm(bmi ~ population, data = dat)
m0 <- lm(bmi ~ 1, data = dat)
anova(m1, m0, test = "F")



