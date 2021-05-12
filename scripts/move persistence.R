###########################################################################
### Estimate link between move persistence and environmental covariates ###
###########################################################################

# load libraries
require(tidyverse)
#remotes::install_github("ianjonsen/mpmm")
require(mpmm)
source("R/discrete_gradient.R")

# import harp seal migration tracks truncated if there were long gaps in transmission
dat <- readRDS("data/move_persistence.rds")
# force animal id to character from factor for mpmm model object
dat <- dat %>% mutate(id = as.character(id))

# run mpmm models seperately for two populations
dat_GS <- dat %>% filter(deployment == "hp4")
dat_NW <- dat %>% filter(deployment == "hp6")
dat_GS <- dat_GS %>% arrange(id, date)
dat_NW <- dat_NW %>% arrange(id, date)

# scale and centre covariates
dat_GS <- dat_GS %>%
  mutate(ice_c = scale(ice, center = T, scale = T),
         depth_c = scale(depth, center = T, scale = T))
dat_NW <- dat_NW %>%
  mutate(ice_c = scale(ice, center = T, scale = T),
         depth_c = scale(depth, center = T, scale = T))

# run mpmm on ice and depth
# allow random slopes on ice
# compare models with LRTs

# Greenland Sea animals
m01 <- mpmm(~ ice_c + depth_c + (ice_c | id),
                data = dat_GS,
                control = mpmm_control(REML = F)) 
m02 <- mpmm(~  ice_c + depth_c + (1 | id),
            data = dat_GS,
            control = mpmm_control(REML = F)) 
anova(m01, m02)

m03 <- mpmm(~  depth_c + (1 | id),
            data = dat_GS,
            control = mpmm_control(REML = F)) 
m04 <- mpmm(~  ice_c + (1 | id),
            data = dat_GS,
            control = mpmm_control(REML = F)) 
anova(m02, m03)
anova(m02, m04)

# Northwest Atlantic animals
m01 <- mpmm(~ ice_c + depth_c + (ice_c | id),
            data = dat_NW,
            control = mpmm_control(REML = F))
m02 <- mpmm(~ ice_c + depth_c + (1 | id),
            data = dat_NW,
            control = mpmm_control(REML = F))
anova(m01, m02)

m03 <- mpmm(~ ice_c + (ice_c | id),
            data = dat_NW,
            control = mpmm_control(REML = F))
anova(m01, m03)

# refit models with REML for prediction
m_GS <- mpmm(~  ice_c + depth_c + (1 | id),
             data = dat_GS,
             control = mpmm_control(REML = T)) 
m_NW <- mpmm(~ ice_c + depth_c + (ice_c | id),
             data = dat_NW,
             control = mpmm_control(REML = T))

####################################
### custom plots for reml models ###
####################################

# for GS
m <- m_GS
terms <- attr(terms(lme4::nobars(m$formula)), "term.labels")
n <- length(terms)
rng <- sapply(1:n, function(i) range(m$fr[, terms[i]]))
xt <- sapply(1:n, function(i) seq(rng[1,i], rng[2,i], l = 200))
xt.mn <- apply(xt, 2, mean)

f_int <- m$par["Intercept","Estimate"]
betas <- sapply(1:n, function(i) m$par[terms[i],"Estimate"])

fxd <- sapply(1:n, function(i) {
  if(n > 2) {
    plogis(f_int + betas[i] * xt[, i] + betas[-i] %*% xt.mn[-i])
  } else if(n > 1){
    plogis(f_int + betas[i] * xt[, i] + betas[-i] * xt.mn[-i])
  } else {
    plogis(f_int + betas * xt)
  }
})

## intercept only random effect
re_ints <- m$par["Intercept", "Estimate"] + m$re$`(Intercept)`
k <- length(re_ints)

re <- lapply(1:n, function(j) {
  if(n > 1) {
    plogis(outer(betas[j] * xt[, j] + betas[-j] * xt.mn[-j], re_ints, FUN = "+"))
  } else {
    plogis(outer(betas * xt, re_ints, FUN = "+"))
  }
})

p <- lapply(1:n, function(j) {
  pdat <- data.frame(x = xt[, j], g = re[[j]])
  pdat <-
    reshape2::melt(
      pdat,
      id.vars = "x",
      value.name = "g",
      variable.name = "Intercept"
    )
  pdat <- data.frame(id = rep(as.character(m$re$id), each = 200), pdat)
  pdat1 <- data.frame(x=xt[,j], y=fxd[,j])
})

rnms <- names(m$re)[!names(m$re) %in% c("id","(Intercept)")]
# check for fixed terms not in random terms
rmiss <- which(!terms %in% rnms)
rpos <- which(terms %in% rnms)

bs <- matrix(0, ncol = n, nrow = k)
bs[, rmiss] <- 0
bs[, rpos] <- unlist(m$re[, rnms])

re_betas <- sapply(1:n, function(i) {
  (betas[i] + bs[, i])
})

re <- lapply(1:n, function(j){
  if(n > 2) {
    plogis(re_ints + re_betas[,j] %o% xt[,j] + as.vector(re_betas[,-j] %*% xt.mn[-j]))
  } else if(n > 1) {
    plogis(re_ints + re_betas[,j] %o% xt[,j] + re_betas[,-j] * xt.mn[-j])
  } else {
    plogis(re_ints + re_betas[,j] %o% xt[,j])
  }
})

j = 1
pdat <- data.frame(x = xt[, j], g = t(re[[j]]))
pdat <- reshape2::melt(
  pdat,
  id.vars = "x",
  value.name = "g",
  variable.name = "re"
)
pdat1 <- data.frame(x=xt[,j], y=fxd[,j])

pdat$x <- (pdat$x * attr(dat_GS$ice_c, 'scaled:scale')) + attr(dat_GS$ice_c, 'scaled:center')
pdat1$x <- (pdat1$x * attr(dat_GS$ice_c, 'scaled:scale')) + attr(dat_GS$ice_c, 'scaled:center')

GS_dat1 <- pdat
GS_dat2 <- pdat1


j = 2
pdat <- data.frame(x = xt[, j], g = t(re[[j]]))
pdat <- reshape2::melt(
  pdat,
  id.vars = "x",
  value.name = "g",
  variable.name = "re"
)
pdat2 <- data.frame(x=xt[,j], y=fxd[,j])

pdat$x <- pdat$x * attr(dat_GS$depth_c, 'scaled:scale') + attr(dat_GS$depth_c, 'scaled:center')
pdat2$x <- pdat2$x * attr(dat_GS$depth_c, 'scaled:scale') + attr(dat_GS$depth_c, 'scaled:center')

GS_dat3 <- pdat #store predictions
GS_dat4 <- pdat2

# repeat for NW
m <- m_NW
terms <- attr(terms(lme4::nobars(m$formula)), "term.labels")
n <- length(terms)
rng <- sapply(1:n, function(i) range(m$fr[, terms[i]]))
xt <- sapply(1:n, function(i) seq(rng[1,i], rng[2,i], l = 200))
xt.mn <- apply(xt, 2, mean)

f_int <- m$par["Intercept","Estimate"]
betas <- sapply(1:n, function(i) m$par[terms[i],"Estimate"])

fxd <- sapply(1:n, function(i) {
  if(n > 2) {
    plogis(f_int + betas[i] * xt[, i] + betas[-i] %*% xt.mn[-i])
  } else if(n > 1){
    plogis(f_int + betas[i] * xt[, i] + betas[-i] * xt.mn[-i])
  } else {
    plogis(f_int + betas * xt)
  }
})

re_ints <- f_int + m$re$`(Intercept)`
k <- length(re_ints)

rnms <- names(m$re)[!names(m$re) %in% c("id","(Intercept)")]
# check for fixed terms not in random terms
rmiss <- which(!terms %in% rnms)
rpos <- which(terms %in% rnms)

bs <- matrix(0, ncol = n, nrow = k)
bs[, rmiss] <- 0
bs[, rpos] <- unlist(m$re[, rnms])

re_betas <- sapply(1:n, function(i) {
  (betas[i] + bs[, i])
})

re <- lapply(1:n, function(j){
  if(n > 2) {
    plogis(re_ints + re_betas[,j] %o% xt[,j] + as.vector(re_betas[,-j] %*% xt.mn[-j]))
  } else if(n > 1) {
    plogis(re_ints + re_betas[,j] %o% xt[,j] + re_betas[,-j] * xt.mn[-j])
  } else {
    plogis(re_ints + re_betas[,j] %o% xt[,j])
  }
})

j = 1
pdat <- data.frame(x = xt[, j], g = t(re[[j]]))
pdat <- reshape2::melt(
  pdat,
  id.vars = "x",
  value.name = "g",
  variable.name = "re"
)
pdat1 <- data.frame(x=xt[,j], y=fxd[,j])

pdat$x <- (pdat$x * attr(dat_NW$ice_c, 'scaled:scale')) + attr(dat_NW$ice_c, 'scaled:center')
pdat1$x <- (pdat1$x * attr(dat_NW$ice_c, 'scaled:scale')) + attr(dat_NW$ice_c, 'scaled:center')

NW_dat1 <- pdat
NW_dat2 <- pdat1

j = 2
pdat <- data.frame(x = xt[, j], g = t(re[[j]]))
pdat <- reshape2::melt(
  pdat,
  id.vars = "x",
  value.name = "g",
  variable.name = "re"
)
pdat2 <- data.frame(x=xt[,j], y=fxd[,j])

pdat$x <- pdat$x * attr(dat_NW$depth_c, 'scaled:scale') + attr(dat_NW$depth_c, 'scaled:center')
pdat2$x <- pdat2$x * attr(dat_NW$depth_c, 'scaled:scale') + attr(dat_NW$depth_c, 'scaled:center')

NW_dat3 <- pdat
NW_dat4 <- pdat2

# convert to dataframe for faceted ggplot
GS_dat1 <- GS_dat1 %>% mutate(Population = "Greenland Sea",
                              Covariate = "Sea Ice Concentration (%)")
GS_dat2 <- GS_dat2 %>% mutate(Population = "Greenland Sea",
                              Covariate = "Sea Ice Concentration (%)")
GS_dat3 <- GS_dat3 %>% mutate(Population = "Greenland Sea",
                              Covariate = "Depth (m)")
GS_dat4 <- GS_dat4 %>% mutate(Population = "Greenland Sea",
                              Covariate = "Depth (m)")
NW_dat1 <- NW_dat1 %>% mutate(Population = "Northwest Atlantic",
                              Covariate = "Sea Ice Concentration (%)")
NW_dat2 <- NW_dat2 %>% mutate(Population = "Northwest Atlantic",
                              Covariate = "Sea Ice Concentration (%)")
NW_dat3 <- NW_dat3 %>% mutate(Population = "Northwest Atlantic",
                              Covariate = "Depth (m)")
NW_dat4 <- NW_dat4 %>% mutate(Population = "Northwest Atlantic",
                              Covariate = "Depth (m)")
pop <- rbind(GS_dat2, GS_dat4, NW_dat2, NW_dat4)
res <- rbind(GS_dat1, GS_dat3, NW_dat1, NW_dat3)

# dummy to pad x axis to range of data
dummy <- tibble(x = c(100, 5000), y = c(1, 1), Covariate = c("Sea Ice Concentration (%)", "Depth (m)"))

# choose palette
wpal <- wesanderson::wes_palette("Zissou1", n = 5, "discrete")
wpal[c(1, 4)]

leg <- ggplot() +
  geom_density(aes(x = x, colour = Population, fill = Population), data = pop) +
  scale_colour_manual("", values = wpal[c(4, 1)]) +
  scale_fill_manual("", values = wpal[c(4, 1)]) +
  theme(legend.key=element_blank()) 

# 2. Save the legend from the dot plot
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  leg +  guides(color = guide_legend(nrow = 1)) + theme(legend.position = "bottom", legend.justification="center")
)

p1 <- ggplot() + 
  theme_bw(base_size = 10) +
  geom_line(aes(x = x, y = g, group = re), data = res, size = 0.2, colour = "grey") + 
  geom_line(aes(x = x, y = y, colour = Population), size = 1, data = pop) + 
  scale_colour_manual("", values = wpal[c(4, 1)]) +
  geom_blank(aes(x = x, y = y), data = dummy) +
  coord_cartesian(ylim = c(0, 1),expand = F) +
  xlab("") +
  ylab(expression(gamma[t])) +
  facet_grid(Population ~ Covariate,
             scales = "free_x",
             switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.spacing = unit(1, "lines"),
        legend.position = "none")

quartz(title = "Panel Plot", width = 6, height = 6)
cowplot::plot_grid(p1, legend, nrow = 2, rel_heights = c(2, .1))
quartz.save(file = "figures/mpmm covariate figure.jpeg", type = "jpeg",
            dev  = dev.cur(), dpi = 500)
dev.off()












