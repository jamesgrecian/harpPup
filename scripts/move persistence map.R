############################
### Map move persistence ###
############################

# example of how to map move persistence to tracks

# load libraries
require(tidyverse)
require(sf)
require(foieGras)
require(patchwork)
source("R/discrete_gradient.R")

# load data
dat <- readRDS("data/move_persistence.rds")

# estimate g from mpm model
fmpm <- foieGras::fit_mpm(dat[,2:5], model = "mpm")

# combine
dat <- dat %>% left_join(grab(fmpm, "fitted"))

# get land shapefile
world_shp <- rnaturalearth::ne_download(scale = 110, type = "land", category = "physical", returnclass = "sf")

# define projection for Northwest Atlantic
prj = "+proj=laea +lat_0=55 +lon_0=-50 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"

p1 <- ggplot() +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10)) +
  labs(title = "Northwest Atlantic", tag = expression(paste("(", italic(a), ")"))) + 
  geom_sf(aes(colour = g), size = 0.5,
          data = dat %>% filter(deployment == "hp6") %>%
            st_as_sf(coords = c('lon', 'lat')) %>% st_set_crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')) +
  geom_sf(aes(), data = world_shp, colour = "dark grey", fill = "dark grey") +
  scale_colour_discrete_gradient(expression(gamma[t]),
                                 colours = viridis::viridis(12),
                                 bins = 12,
                                 limits = c(0, 1),
                                 breaks = seq(0, 1, .25),
                                 guide = guide_colourbar(nbin = 500,
                                                         raster = T,
                                                         frame.colour = "black",
                                                         ticks.colour = "black",
                                                         frame.linewidth = 1,
                                                         barwidth = .75,
                                                         barheight = 15,
                                                         direction = "vertical",
                                                         title.position = "right", #or "right"
                                                         title.theme = element_text(angle = 0,
                                                                                    hjust = 0.5,
                                                                                    family = "serif"))) +
  
  theme(legend.position = "right") +
  coord_sf(xlim = c(-1500, 1500), ylim = c(-1500, 2500), crs = prj, expand = T) +
  scale_x_continuous(breaks = seq(from = -130, to = 130, by = 10)) +
  xlab("") +ylab("")

# define projection for Greenland Sea
prj = "+proj=laea +lat_0=75 +lon_0=-10 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"

p2 <- ggplot() +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10)) +
  labs(title = "Greenland Sea", tag = expression(paste("(", italic(b), ")"))) + 
  geom_sf(aes(colour = g), size = 0.5,
          data = dat %>% filter(deployment == "hp4") %>%
            st_as_sf(coords = c('lon', 'lat')) %>% st_set_crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')) +
  geom_sf(aes(), data = world_shp, colour = "dark grey", fill = "dark grey") +
  scale_colour_discrete_gradient(expression(gamma[t]),
                                 colours = viridis::viridis(12),
                                 bins = 12,
                                 limits = c(0, 1),
                                 breaks = seq(0, 1, .25),
                                 guide = guide_colourbar(nbin = 500,
                                                         raster = T,
                                                         frame.colour = "black",
                                                         ticks.colour = "black",
                                                         frame.linewidth = 1,
                                                         barwidth = .75,
                                                         barheight = 15,
                                                         direction = "vertical",
                                                         title.position = "right", #or "right"
                                                         title.theme = element_text(angle = 0,
                                                                                    hjust = 0.5,
                                                                                    family = "serif"))) +
  
  theme(legend.position = "right") +
  coord_sf(xlim = c(-1000, 1600), ylim = c(-1200, 1400), crs = prj, expand = T) +
  scale_x_continuous(breaks = seq(from = -180, to = 180, by = 10)) +
  xlab("") + ylab("")

quartz(title = "Panel Plot", width = 8, height = 4)
p1 + p2 + plot_layout(guides = "collect")
quartz.save(file = "figures/move persistence map.jpeg", type = "jpeg",
            dev  = dev.cur(), dpi = 500)
dev.off()
