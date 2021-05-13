#################################################################################
### Generate vertical histogram plots to highlight differences in water depth ###
#################################################################################

# load libraries
require(tidyverse)
require(cowplot)

# load data
dat <- readRDS("data/move_persistence.rds")

# append population name
dat <- dat %>% mutate(population = case_when(deployment == "hp4" ~ "Greenland Sea",
                                             deployment == "hp6" ~ "Northwest Atlantic"))
dat <- dat %>% mutate(population = factor(population))

# append days since commencement of diving
dat <- dat %>% group_by(id) %>% mutate(days = difftime(date, first(date), units = "days")) %>% mutate(days = as.numeric(days)) %>% ungroup()

# choose palette
wpal <- wesanderson::wes_palette("Zissou1", n = 5, "discrete")
wpal[c(1, 4)]

# create legend
leg <- ggplot() +
  geom_density(aes(x = depth, colour = population, fill = population), data = dat %>% filter(days < 25)) +
  scale_colour_manual("", values = wpal[c(4, 1)]) +
  scale_fill_manual("", values = wpal[c(4, 1)]) +
  theme(legend.key=element_blank()) 

# 2. Save the legend from the dot plot
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  leg + theme(legend.position = "bottom")
)

p1 <- ggplot() +
  theme_bw(base_size = 10) +
  #  theme(text = element_text(family = "serif")) +
  geom_density(aes(x = depth, y = -..scaled..), fill = wpal[4], dat %>% filter(population == "Greenland Sea") %>% filter(days < 25), outline.type = "full") +
  geom_density(aes(x = depth, y = ..scaled..), fill = wpal[1], dat %>% filter(population == "Northwest Atlantic") %>% filter(days < 25), outline.type = "full") +
  xlab("Depth (m)") + 
  ylab("Scaled density") +
  xlim(0, 5000) +
  ggtitle("Water depth during first 25 days") +
  theme(plot.title = element_text(size = 10)) +
  coord_flip()

p2 <- ggplot() +
  theme_bw(base_size = 10) +
  #  theme(text = element_text(family = "serif")) +
  geom_density(aes(x = depth, y = -..scaled..), fill = wpal[4], dat %>% filter(population == "Greenland Sea") %>% filter(days > 25), outline.type = "full") +
  geom_density(aes(x = depth, y = ..scaled..), fill = wpal[1], dat %>% filter(population == "Northwest Atlantic") %>% filter(days > 25), outline.type = "full") +
  xlab("Depth (m)") + 
  ylab("Scaled density") +
  xlim(0, 5000) +
  ggtitle("Water depth after first 25 days") +
  theme(plot.title = element_text(size = 10)) +
  coord_flip()

#https://wilkelab.org/cowplot/articles/shared_legends.html
quartz(title = "Panel Plot", width = 6, height = 4)
prow <- plot_grid(p1,
                  p2,
                  nrow = 1,
                  labels = c("(a)", "(b)"),
                  label_fontface = "italic",
                  label_size = 12)
plot_grid(prow, legend, nrow = 2, rel_heights = c(1, .1))
quartz.save(file = "figures/depth density plots.jpeg", type = "jpeg",
            dev  = dev.cur(), dpi = 500)
dev.off()

# ends