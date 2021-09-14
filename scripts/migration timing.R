##########################################################################
### Load regularised data and estimate timing of migration uing 1d HMM ###
##########################################################################

# load libraries
require(tidyverse)
require(lubridate)
require(depmixS4)

# load regularised location data for both deployments
dat <- readRDS("data/migration_timing.rds")

# calculate difference in latitude between successive steps
dat <- dat %>%
  group_by(id) %>%
  mutate(lat_diff = lat - lag(lat, default = first(lat), order_by = date))

# run HMM
mod = depmix(lat_diff ~ 1,
             nstates = 2,
             transition = ~ 1,
             ntimes = dat %>% group_by(id) %>% tally() %>% pull(n),
             family = gaussian(),
             data = dat)

fd <- fit(mod)
fdp <- posterior(fd)
dat$states <- fdp$state # append states

# identify switches between states
dat <- dat %>% group_by(id) %>% dplyr::select(deployment, id, date, lon, lat, lat_diff, states) %>% 
  mutate(select = states - lag(states, default = first(states), order_by = date))

# pull dates of switches and calculate length of time spent in each state
res <- dat %>% filter(select != 0) %>% mutate(lag = difftime(date, lag(date), units = "days")) 
res %>% filter(yday(date) > 150) %>% filter(lag > 5) %>% print(n = Inf)

# highlight range of departure dates
date_poly <- tibble(x1 = c(ymd_hms("2017-04-23 00:00:00"), ymd_hms("2019-06-05 00:00:00")),
                    x2 = c(ymd_hms("2017-06-16 00:00:00"), ymd_hms("2019-07-21 00:00:00")),
                    y1 = c(-Inf, -Inf),
                    y2 = c(Inf, Inf),
                    deployment = c("Greenland Sea",
                                   "Northwest Atlantic"))

p1 <- ggplot() +
  theme_bw(base_size = 14) +
  ylab("Latitude") + xlab("Month") +
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, group = deployment), alpha = 0.25, data = date_poly) +
  geom_point(aes(x = date, y = lat, colour = factor(states)), data = dat) + 
  scale_color_viridis_d("state") +
  facet_wrap(~deployment, scales = "free") +
  theme(legend.position = "none")

p2 <- ggplot() + 
  theme_bw(base_size = 14) +
  geom_density(aes(x = abs(lat_diff), fill = "state 2"), alpha = .75, data = dat %>% filter(states == 2)) +
  geom_density(aes(x = abs(lat_diff), fill = "state 1"), alpha = .75, data = dat %>% filter(states == 1)) +
  scale_fill_viridis_d("") +
  scale_x_continuous(limits = c(-.05, .6)) +
  ylab("density") + xlab(expression(Latitudinal~displacement~(degree*N~d^{-1}))) +
  facet_wrap(~deployment) +
  theme(legend.position = "bottom")

require(patchwork)
quartz(title = "Panel Plot", width = 9, height = 8)
p1 / p2 +  plot_layout(heights = c(2, 1)) + plot_annotation(tag_levels = "a")
quartz.save(file = "figures/supplementary migration timing.jpeg", type = "jpeg",
            dev  = dev.cur(), dpi = 500)
dev.off()

# ends