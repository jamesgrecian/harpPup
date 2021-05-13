###############################
### Dive behaviour analysis ###
###############################

# load libraries
require(tidyverse)
require(mgcv)
require(gratia)
require(cowplot)
require(ggExtra)

############################
### Individual dive data ###
############################

# load individual dive data
dat <- readRDS("data/individual_dives.rds")

##########################
### Average dive depth ###
##########################

# for a great post on nested interaction terms in GAM see:
# https://stats.stackexchange.com/questions/403772/different-ways-of-modelling-interactions-between-continuous-and-categorical-pred

# full model
# fixed effects of population
# spline interactions between time and population
# random intercept and random slope terms for individual
# AR1 correlation structure
# drop the three dives greater than 600 metres

m01 <- gamm(max.dep ~ population +
                      s(days, k = 6, bs = "cs") +
                      s(days, by = population, k = 6, bs = "cs", m = 1) +
                      s(id, bs = "re") +
                      s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dat %>% filter(max.dep < 600),
            family = Gamma(link = log),
            method = "ML")

m02 <- gamm(max.dep ~ population +
                      s(days, k = 6, bs = "cs") +
                      s(id, bs = "re") +
                      s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dat %>% filter(max.dep < 600),
            family = Gamma(link = log),
            method = "ML")

AIC(m01$lme, m02$lme)
AIC(m01$lme) - AIC(m02$lme)
plot(acf(resid(m01$lme, type = "normalized")))
gratia::appraise(m01$gam)

# refit with REML for prediction
m01a <- gamm(max.dep ~ population +
                       s(days, k = 6, bs = "cs") +
                       s(days, by = population, k = 6, bs = "cs", m = 1) +
                       s(id, bs = "re") +
                       s(days, id, bs = "re", k = c(6, 10)),
             correlation = corAR1(form = ~ 1 | id),
             data = dat %>% filter(max.dep < 600),
             family = Gamma(link = log),
             method = "REML")

#############################
### Average dive duration ###
#############################

# full model
# fixed effects of population
# spline interactions between time and population
# random intercept and random slope terms for individual
# AR1 correlation structure
# drop the dives longer than 1000 seconds

m03 <- gamm(dive.dur ~ population +
                       s(days, k = 6, bs = "cs") +
                       s(days, by = population, k = 6, bs = "cs", m = 1) +
                       s(id, bs = "re") +
                       s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dat %>% filter(dive.dur < 1000),
            family = tw,
            method = "ML")

m04 <- gamm(dive.dur ~ population +
                       s(days, k = 6, bs = "cs") +
                       s(id, bs = "re") +
                       s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dat %>% filter(dive.dur < 1000),
            family = tw,
            method = "ML")

AIC(m03$lme, m04$lme)
AIC(m03$lme) - AIC(m04$lme)
plot(acf(resid(m03$lme, type = "normalized")))
gratia::appraise(m03$gam)

# refit with REML for prediction
m03a <- gamm(dive.dur ~ population +
                        s(days, k = 6, bs = "cs") +
                        s(days, by = population, k = 6, bs = "cs", m = 1) +
                        s(id, bs = "re") +
                        s(days, id, bs = "re", k = c(6, 10)),
             correlation = corAR1(form = ~ 1 | id),
             data = dat %>% filter(dive.dur < 1000),
             family = tw,
             method = "REML")

################################
### Average surface duration ###
################################

# full model
# fixed effects of population
# spline interactions between time and population
# random intercept and random slope terms for individual
# AR1 correlation structure

m05 <- gamm(surf.dur ~ population +
                       s(days, k = 6, bs = "cs") +
                       s(days, by = population, k = 6, bs = "cs", m = 1) +
                       s(id, bs = "re") +
                       s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dat,
            family = tw,
            method = "ML") 

# drop interaction between population and time
m06 <- gamm(surf.dur ~ population +
                       s(days, k = 6, bs = "cs") +
                       s(id, bs = "re") +
                       s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dat,
            family = tw,
            method = "ML") 
AIC(m05$lme, m06$lme)
AIC(m05$lme) - AIC(m06$lme) # no support for population interaction

# drop fixed effect of population
m07 <- gamm(surf.dur ~ s(days, k = 6, bs = "cs") +
                       s(id, bs = "re") +
                       s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dat,
            family = tw,
            method = "ML")
AIC(m06$lme) - AIC(m07$lme) # no support for population effect

AIC(m05$lme, m06$lme, m07$lme)

# model 7 is final model
plot(acf(resid(m07$lme, type = "normalized")))
gratia::appraise(m07$gam)

# refit with REML for prediction
m07a <- gamm(surf.dur ~ s(days, k = 6, bs = "cs") +
                        s(id, bs = "re") +
                        s(days, id, bs = "re", k = c(6, 10)),
             correlation = corAR1(form = ~ 1 | id),
             data = dat,
             family = tw,
             method = "REML") 

#########################
### Dive summary data ###
#########################

# load individual dive data
dive_sum <- readRDS("data/dive_summaries.rds")

#################
### Dive rate ###
#################

m08 <- gamm(n.cycles ~ population + 
                       s(days, k = 6, bs = "cs") +
                       s(days, by = population, k = 6, bs = "cs", m = 1) +
                       s(id, bs = "re") +
                       s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dive_sum,
            method = "ML") 

m09 <- gamm(n.cycles ~ population + 
                       s(days, k = 6, bs = "cs") +
                       s(id, bs = "re") +
                       s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dive_sum,
            method = "ML") 

AIC(m08$lme, m09$lme)
AIC(m08$lme) - AIC(m09$lme) # no support for difference between populations

m10 <- gamm(n.cycles ~ s(days, k = 6, bs = "cs") +
                       s(id, bs = "re") +
                       s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dive_sum,
            method = "ML") 
AIC(m08$lme, m09$lme, m10$lme) # no support for retaining population fixed effect
AIC(m09$lme) - AIC(m10$lme) # no support for difference between populations

plot(acf(resid(m08$lme, type = "normalized")))
gratia::appraise(m08$gam)

# refit with REML for prediction
m10a <- gamm(n.cycles ~ s(days, k = 6, bs = "cs") +
                        s(id, bs = "re") +
                        s(days, id, bs = "re", k = c(6, 10)),
             correlation = corAR1(form = ~ 1 | id),
             data = dive_sum,
             method = "REML") 

######################
### Max dive depth ###
######################

m11 <- gamm(max.depth ~ population +
                        s(days, k = 6, bs = "cs") +
                        s(days, by = population, k = 6, bs = "cs", m = 1) +
                        s(id, bs = "re") +
                        s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dive_sum,
            method = "ML") 

m12 <- gamm(max.depth ~ population +
                        s(days, k = 6, bs = "cs") +
                        s(id, bs = "re") +
                        s(days, id, bs = "re", k = c(6, 10)),
            correlation = corAR1(form = ~ 1 | id),
            data = dive_sum,
            method = "ML") 

AIC(m11$lme, m12$lme) # support for including population effect
AIC(m11$lme) - AIC(m12$lme)
gratia::appraise(m11$gam)

# refit with REML for prediction
m11a <- gamm(max.depth ~ population +
                         s(days, k = 6, bs = "cs") +
                         s(days, by = population, k = 6, bs = "cs", m = 1) +
                         s(id, bs = "re") +
                         s(days, id, bs = "re", k = c(6, 10)),
             correlation = corAR1(form = ~ 1 | id),
             data = dive_sum,
             method = "REML") 

#########################
### Max dive duration ###
#########################

# drop 7 dive summary periods where max duration exceeded 1000 seconds
# filter summary periods where max duration is 0

m13 <- gamm(max.dur ~ population +
                      s(days, k = 6, bs = "cs") +
                      s(days, by = population, k = 6, bs = "cs", m = 1) +
                      s(id, bs = "re") +
                      s(days, id, bs = "re", k = c(6, 10)),
            data = dive_sum %>% filter(max.dur < 1000) %>% filter(max.dur > 0),
            family = tw,
            method = "ML") 

m14 <- gamm(max.dur ~ population +
                      s(days, k = 6, bs = "cs") +
                      s(id, bs = "re") +
                      s(days, id, bs = "re", k = c(6, 10)),
            data = dive_sum %>% filter(max.dur < 1000) %>% filter(max.dur > 0),
            family = tw,
            method = "ML") 

AIC(m13$lme, m14$lme)
AIC(m13$lme) - AIC(m14$lme)
gratia::appraise(m13$gam)

# refit with REML for prediction
m13a <- gamm(max.dur ~ population +
                       s(days, k = 6, bs = "cs") +
                       s(days, by = population, k = 6, bs = "cs", m = 1) +
                       s(id, bs = "re") +
                       s(days, id, bs = "re", k = c(6, 10)),
             data = dive_sum %>% filter(max.dur < 1000) %>% filter(max.dur > 0),
             family = tw,
             method = "REML") 

#########################################################
### Generate predictions from each model for plotting ###
#########################################################

### predictions for dive depth ###

# generate prediction dataframe
newdata = expand.grid(days = 0:100,
                      id = unique(dat$id))
ids <- dat %>% group_by(id) %>% slice(1) %>% dplyr::select(id, population) %>% ungroup()
newdata <- newdata %>% left_join(ids, by = c("id" = "id"))
# append model predictions
newdata$pop <- predict.gam(m01a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"))
newdata$se.fit <- predict.gam(m01a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"), se.fit = T)$se.fit
newdata$fit <- predict.gam(m01a$gam, newdata, type = "response")

max_depths <- dat %>% mutate(day = round(days)) %>% group_by(day) %>% summarise(depth95 = quantile(max.dep, probs = 0.95))
max_depths <- max_depths %>% mutate(max = cummax(depth95))

wpal <- wesanderson::wes_palette("Zissou1", n = 5, "discrete")
wpal[c(1, 4)]

p1 <- ggplot() +
  theme_bw(base_size = 10) +
  #  theme(text = element_text(family = "serif")) +
  geom_line(aes(x = days, y = fit, group = id, colour = population), data = newdata, alpha = 0.4) +
  geom_ribbon(aes(x = days, ymin = pop - (1.96* se.fit), ymax = pop + (1.96* se.fit), group = population, fill = population), data = newdata, alpha = 0.5) +
  geom_line(aes(x = days, y = pop, group = population, colour = population), data = newdata, size = 1) +
  geom_path(aes(x = day, y = max), data = max_depths) +
  ylab("Dive depth (m)") + xlab("Days since start of diving") +
  scale_colour_manual(values = wpal[c(4, 1)]) + scale_fill_manual(values = wpal[c(4, 1)]) +
  geom_point(aes(x = days, y = max.dep, colour = population), data = dat, alpha = 0) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 100), oob = scales::squish, expand = c(0, 0)) + # ggMarginal doesn't work with coord_cartesian
  scale_y_continuous(limits = c(0, 200), oob = scales::squish, expand = c(0, 0))


### predictions for dive duration ###

# generate prediction dataframe
newdata = expand.grid(days = 0:100,
                      id = unique(dat$id))
ids <- dat %>% group_by(id) %>% slice(1) %>% dplyr::select(id, population) %>% ungroup()
newdata <- newdata %>% left_join(ids, by = c("id" = "id"))
# append model predictions
newdata$pop <- predict.gam(m03a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"))/60
newdata$se.fit <- predict.gam(m03a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"), se.fit = T)$se.fit/60
newdata$fit <- predict.gam(m03a$gam, newdata, type = "response")/60

max_dur <- dat %>% mutate(day = round(days)) %>% group_by(day) %>% summarise(dur95 = quantile(dive.dur, probs = 0.95))
max_dur <- max_dur %>% mutate(max = cummax(dur95))

p2 <- ggplot() +
  theme_bw(base_size = 10) +
  #  theme(text = element_text(family = "serif")) +
  geom_point(aes(x = days, y = dive.dur/60, colour = population), data = dat, alpha = 0) +
  theme(legend.position = "none") +
  geom_line(aes(x = days, y = fit, group = id, colour = population), data = newdata, alpha = 0.4) +
  geom_ribbon(aes(x = days, ymin = pop - (1.96* se.fit), ymax = pop + (1.96* se.fit), group = population, fill = population), data = newdata, alpha = 0.5) +
  geom_line(aes(x = days, y = pop, group = population, colour = population), data = newdata, size = 1) +
  geom_hline(yintercept = 4.3, linetype = "dashed") +
  geom_line(aes(x = day, y = max/60), data = max_dur) +
  ylab("Dive duration (min)") + xlab("Days since start of diving") +
  scale_colour_manual(values = wpal[c(4, 1)]) + scale_fill_manual(values = wpal[c(4, 1)]) +
  scale_x_continuous(limits = c(0, 100), oob=scales::squish, expand = c(0,0)) + # ggMarginal doesn't work with coord_cartesian
  scale_y_continuous(limits = c(0, 8), oob=scales::squish, expand = c(0,0))


### predictions for surface duration

# generate prediction dataframe
newdata = expand.grid(days = 0:100,
                      id = unique(dat$id))
ids <- dat %>% group_by(id) %>% slice(1) %>% dplyr::select(id, population) %>% ungroup()
newdata <- newdata %>% left_join(ids, by = c("id" = "id"))

# append model predictions
newdata$pop <- predict.gam(m07a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"))/60
newdata$se.fit <- predict.gam(m07a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"), se.fit = T)$se.fit/60
newdata$fit <- predict.gam(m07a$gam, newdata, type = "response")/60

p3 <- ggplot() +
  theme_bw(base_size = 10) +
  #  theme(text = element_text(family = "serif")) +
  geom_point(aes(x = days, y = surf.dur/60, colour = population), data = dat, alpha = 0) +
  theme(legend.position = "none") +
  geom_line(aes(x = days, y = fit, group = id, colour = population), data = newdata, alpha = 0.4) +
  geom_ribbon(aes(x = days, ymin = pop - (1.96* se.fit), ymax = pop + (1.96* se.fit)), data = newdata, alpha = 0.5) +
  geom_line(aes(x = days, y = pop), data = newdata, size = 1) +
  ylab("Surface duration (min)") + xlab("Days since start of diving") +
  scale_colour_manual(values = wpal[c(4, 1)]) + scale_fill_manual(values = wpal[c(4, 1)]) +
  scale_x_continuous(limits = c(0, 100), oob=scales::squish, expand = c(0,0)) + # ggMarginal doesn't work with coord_cartesian
  scale_y_continuous(limits = c(0, 8), oob=scales::squish, expand = c(0,0))


### predictions for dive rate

# generate prediction dataframe
newdata = expand.grid(days = 0:100,
                      id = unique(dat$id))
ids <- dive_sum %>% group_by(id) %>% slice(1) %>% dplyr::select(id, population) %>% ungroup()
newdata <- newdata %>% left_join(ids, by = c("id" = "id"))

# only predict to day 80 for Greenland Sea
newdata <- newdata %>% filter(case_when(population == "Greenland Sea" ~ days < 81,
                                        population == "Northwest Atlantic" ~ days < 101))
# append model predictions
newdata$pop <- predict.gam(m10a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"))
newdata$se.fit <- predict.gam(m10a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"), se.fit = T)$se.fit
newdata$fit <- predict.gam(m10a$gam, newdata, type = "response")

p4 <- ggplot() +
  theme_bw(base_size = 10) +
  #  theme(text = element_text(family = "serif")) +
  geom_point(aes(x = days, y = n.cycles/6, colour = population), data = dive_sum, alpha = 0) +
  geom_line(aes(x = days, y = fit/6, group = id, colour = population), data = newdata, alpha = 0.4) +
  geom_ribbon(aes(x = days, ymin = (pop - (1.96* se.fit))/6, ymax = (pop + (1.96* se.fit))/6), data = newdata, alpha = 0.5) +
  geom_line(aes(x = days, y = pop/6, group = population), data = newdata, size = 1) +
  theme(legend.position = "none") +
  ylab("Dive rate (dives per h)") + xlab("Days since start of diving") +
  scale_colour_manual(values = wpal[c(4, 1)]) + scale_fill_manual(values = wpal[c(4, 1)]) +
  scale_x_continuous(limits = c(0, 100), oob = scales::squish, expand = c(0, 0)) + # ggMarginal doesn't work with coord_cartesian
  scale_y_continuous(limits = c(0, 20), oob = scales::squish, expand = c(0, 0))


### predictions for max depth

# generate prediction dataframe
newdata = expand.grid(days = 0:100,
                      id = unique(dat$id))
ids <- dive_sum %>% group_by(id) %>% slice(1) %>% dplyr::select(id, population) %>% ungroup()
newdata <- newdata %>% left_join(ids, by = c("id" = "id"))

# only predict to day 80 for Greenland Sea
newdata <- newdata %>% filter(case_when(population == "Greenland Sea" ~ days < 81,
                                        population == "Northwest Atlantic" ~ days < 101))
# append model predictions
newdata$pop <- predict.gam(m11a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"))
newdata$se.fit <- predict.gam(m11a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"), se.fit = T)$se.fit
newdata$fit <- predict.gam(m11a$gam, newdata, type = "response")

max_depths2 <- dive_sum %>% mutate(day = round(days)) %>% group_by(day) %>% summarise(depth95 = quantile(max.depth, probs = 0.95))
max_depths2 <- max_depths2 %>% mutate(max = cummax(depth95))

p5 <- ggplot() +
  theme_bw(base_size = 10) +
  #  theme(text = element_text(family = "serif")) +
  geom_point(aes(x = days, y = max.depth, colour = population), data = dive_sum, alpha = 0) +
  geom_line(aes(x = days, y = fit, group = id, colour = population), data = newdata, alpha = 0.4) +
  geom_ribbon(aes(x = days, ymin = pop - (1.96* se.fit), ymax = pop + (1.96* se.fit), group = population, fill = population), data = newdata, alpha = 0.5) +
  geom_line(aes(x = days, y = pop, group = population, colour = population), data = newdata, size = 1) +
  theme(legend.position = "none") +
  geom_line(aes(x = day, y = max), data = max_depths2) +
  ylab("Max dive depth (m)") + xlab("Days since start of diving") +
  scale_colour_manual(values = wpal[c(4, 1)]) + scale_fill_manual(values = wpal[c(4, 1)]) +
  scale_x_continuous(limits = c(0, 100), oob = scales::squish, expand = c(0, 0)) + # ggMarginal doesn't work with coord_cartesian
  scale_y_continuous(limits = c(0, 200), oob = scales::squish, expand = c(0, 0))


### predictions for max duration

# generate prediction dataframe
newdata = expand.grid(days = 0:100,
                      id = unique(dat$id))
ids <- dive_sum %>% group_by(id) %>% slice(1) %>% dplyr::select(id, population) %>% ungroup()
newdata <- newdata %>% left_join(ids, by = c("id" = "id"))

# only predict to day 80 for Greenland Sea
newdata <- newdata %>% filter(case_when(population == "Greenland Sea" ~ days < 81,
                                        population == "Northwest Atlantic" ~ days < 101))
# append model predictions
newdata$pop <- predict.gam(m13a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"))/60
newdata$se.fit <- predict.gam(m13a$gam, newdata, type = "response", exclude = c("s(days,id)", "s(id)"), se.fit = T)$se.fit/60
newdata$fit <- predict.gam(m13a$gam, newdata, type = "response")/60

max_dur2 <- dive_sum %>% mutate(day = round(days)) %>% group_by(day) %>% summarise(dur95 = quantile(max.dur, probs = 0.95))
max_dur2 <- max_dur2 %>% mutate(max = cummax(dur95))

p6 <- ggplot() +
  theme_bw(base_size = 10) +
  #  theme(text = element_text(family = "serif")) +
  geom_point(aes(x = days, y = max.dur/60, colour = population), data = dive_sum %>% filter(max.dur > 0), alpha = 0) +
  geom_line(aes(x = days, y = fit, group = id, colour = population), data = newdata, alpha = 0.4) +
  geom_ribbon(aes(x = days, ymin = (pop - (1.96* se.fit)), ymax = (pop + (1.96* se.fit)), group = population, fill = population), data = newdata, alpha = 0.5) +
  geom_line(aes(x = days, y = pop, group = population, colour = population), data = newdata, size = 1) +
  geom_hline(yintercept = 4.3, linetype = "dashed") +
  theme(legend.position = "none") +
  geom_line(aes(x = day, y = max/60), data = max_dur2) +
  ylab("Max dive duration (min)") + xlab("Days since start of diving") +
  scale_colour_manual(values = wpal[c(4, 1)]) + scale_fill_manual(values = wpal[c(4, 1)]) +
  scale_x_continuous(limits = c(0, 100), oob = scales::squish, expand = c(0, 0)) + # ggMarginal doesn't work with coord_cartesian
  scale_y_continuous(limits = c(0, 8), oob = scales::squish, expand = c(0, 0))

### Plot layout

leg <- ggplot() +
  geom_point(aes(x = days, y = max.dur, colour = population), data = dive_sum) +
  scale_colour_manual("", values = wpal[c(4, 1)]) +
  theme(legend.key=element_blank()) 

# 2. Save the legend from the dot plot
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  leg + theme(legend.position = "bottom")
)

#https://wilkelab.org/cowplot/articles/shared_legends.html
quartz(title = "Panel Plot", width = 9, height = 6)
prow <- plot_grid(ggExtra::ggMarginal(p1,
                                      type = "density",
                                      margins = "both",
                                      groupFill = T,
                                      size = 5,
                                      xparams = list(size = 0.2),
                                      yparams = list(size = 0.2)),
                  ggExtra::ggMarginal(p2,
                                      type = "density",
                                      margins = "both",
                                      groupFill = T,
                                      size = 5,
                                      xparams = list(size = 0.2),
                                      yparams = list(size = 0.2)),
                  ggExtra::ggMarginal(p3,
                                      type = "density",
                                      margins = "both",
                                      groupFill = T,
                                      size = 5,
                                      xparams = list(size = 0.2),
                                      yparams = list(size = 0.2)),
                  ggExtra::ggMarginal(p5,
                                      type = "density",
                                      margins = "both",
                                      groupFill = T,
                                      size = 5,
                                      xparams = list(size = 0.2),
                                      yparams = list(size = 0.2)),
                  ggExtra::ggMarginal(p6,
                                      type = "density",
                                      margins = "both",
                                      groupFill = T,
                                      size = 5,
                                      xparams = list(size = 0.2),
                                      yparams = list(size = 0.2)),
                  ggExtra::ggMarginal(p4,
                                      type = "density",
                                      margins = "both",
                                      groupFill = T,
                                      size = 5,
                                      xparams = list(size = 0.2),
                                      yparams = list(size = 0.2)),
                  nrow = 2,
                  ncol = 3,
                  labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
                  label_fontface = "italic",
                  label_size = 12)

plot_grid(prow, legend, nrow = 2, rel_heights = c(1, .1))

quartz.save(file = "figures/combo dive stats 20200513.jpeg", type = "jpeg",
            dev  = dev.cur(), dpi = 500)
dev.off()
