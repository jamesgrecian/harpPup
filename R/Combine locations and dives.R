############################################################################
### Understanding harp seal ontogeny from location and dive summary data ###
############################################################################

# Load libraries
require(Hmisc)
require(tidyverse)
require(foieGras)
require(lubridate)
require(sf)

# Download data from SMRU and unzip .mdb
download.file("http://gatty:SOI@www.smru.st-andrews.ac.uk/protected/hp6/db/hp6.zip", "~/hp6.zip")
unzip("~/hp6.zip", exdir = "~/", overwrite = T)

# Load in movement data from .mdb
locs <- Hmisc::mdb.get("~/hp6.mdb",  tables = "diag") %>% as_tibble()
locs <- locs %>% haven::zap_label()
locs <- locs %>% dplyr::select("REF", "D.DATE", "LQ", "LON", "LAT", "SEMI.MAJOR.AXIS", "SEMI.MINOR.AXIS", "ELLIPSE.ORIENTATION")

locs <- locs %>% rename(id = REF,
                      date = D.DATE,
                      lc = LQ,
                      lon = LON,
                      lat = LAT,
                      smaj = SEMI.MAJOR.AXIS,
                      smin = SEMI.MINOR.AXIS,
                      eor = ELLIPSE.ORIENTATION)

locs <- locs %>% mutate(date = mdy_hms(date, tz = "UTC"))

# Recode location class for prefilter algoritm
locs <- locs %>% mutate(lc = recode_factor(factor(lc),
                                         `3` = "3",
                                         `2` = "2",
                                         `1` = "1",
                                         `0` = "0",
                                         `-1` = "A",
                                         `-2` = "B",
                                         `-9` = "Z"))
# Convert names to character
locs <- locs %>% mutate(id = as.character(id))

locs <- locs[!locs$lc == "Z",] # drop the Z locations
locs <- locs[locs$lat > 40,] # drop locations below 30 N
locs <- locs[locs$lat < 85,] # drop locations above 85 N
locs <- locs[locs$lon > -100,] # drop locations below 100 W
locs <- locs[locs$lon < 100,] # drop locations above 100 E

# hp6-L752-19 was recovered dead
locs <- locs %>% filter(case_when(id == "hp6-L752-19" ~ date < "2019-07-29",
                                  id != "hp6-L752-19" ~ date > min(date))) # remove dead locations

# Find 6 hour intervals so interpolated locations match the dive data
source("R/SRDL_time.R")
times <- SRDL_time(locs)

# Regularise the location data to match the 6 hour dive summary data 
# Fit a continuous time random walk using the Argos least squares data
fit <- fit_ssm(locs, model = "rw", time.step = times)
fmap(fit, "predicted")

# Extract fitted values from model
plocs <- grab(fit, "predicted", as_sf = F)

# Load in dive data from .mdb
dives <- Hmisc::mdb.get("~/hp6.mdb",  tables = "summary") %>% as_tibble()
dives <- dives %>% haven::zap_label()

# Format date time
dives <- dives %>% mutate(S.DATE = mdy_hms(S.DATE, tz = "UTC"))
dives <- dives %>% mutate(E.DATE = mdy_hms(E.DATE, tz = "UTC"))

# Select the useful (non NA) columns
dives <- dives %>% dplyr::select("REF",
                             "PTT",
                             "S.DATE",
                             "SURF.TM",
                             "DIVE.TM",
                             "HAUL.TM",
                             "N.CYCLES",
                             "AV.DEPTH",
                             "MAX.DEPTH",
                             "AV.DUR",
                             "SD.DUR",
                             "MAX.DUR",
                             "SWIM.EFF.DESC",
                             "SWIM.EFF.ASC",
                             "SWIM.EFF.WHOLE",
                             "SECS.DESC",
                             "SECS.ASC",
                             "PITCH.DESC",
                             "PITCH.ASC")

# combine regularised locations with dive data
dat <- left_join(plocs, dives, by = c("id" = "REF", "date" = "S.DATE"))
saveRDS(dat, "data/dat.rds")

# ends
