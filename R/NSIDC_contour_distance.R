##################################################
### Extract  distance to NSIDC sea ice contour ###
##################################################

# Given a tibble of animal location data, a path to pull files for and the contour of interest
# Import NSIDC daily data, find contour and estimate distance to time matched points
# st_distance returns distance in metres so /1000 to convert to km

NSIDC_contour_distance <- function(df, path, contour){
  
  ## build a file-db
  files <- tibble(fullname = list.files(normalizePath(path), full.names = TRUE, pattern = "tif$"), 
                  date = as.POSIXct(strptime(basename(fullname), "N_%Y%m%d", tz = "UTC")))
  
  ## map our data to this file-db (we might have multiple points per file)
  df$fullname <- files$fullname[findInterval(df$date, files$date)]
  
  ## load a dummy raster to pull projection info from
  (rdummy <- raster(df$fullname[1]))
  
  ## project the query points
  df[c("X", "Y")] <- as_tibble(rgdal::project(as.matrix(df[c("lon", "lat")]), projection(rdummy)))
  
  ## set up progress bar based on number of files that will be loaded
  pb <- progress_estimated(length(unique(df$fullname)))

  df <- purrr::map_df(split(df, df$fullname)[unique(df$fullname)], 
                      function(.x) {
                        pb$tick()$print()
                        r <- raster(.x$fullname[1])
                        r[r == 2510] <- 1000 # make pole hole 100% ice cover
                        r[r > 1000] <- NA
                        r <- r/10
                        c <- rasterToContour(r, levels = contour) %>% st_as_sf()
                        .x["dist"] <- st_distance(.x %>% st_as_sf(coords = c("X", "Y")) %>% st_set_crs(projection(r)),
                                                  c %>% st_as_sf(), by_element = T) %>% as.vector()
                        .x["dist"] <- .x["dist"]/1000 # convert from m to km
                        .x
              })
  
  # check order of dataframe
  df <- df %>% arrange(id, date)
  df <- df %>% dplyr::select(-c("fullname", "X", "Y"))
  return(df)
}

# This version uses a loop rather than the map_df function
# Its slower for large datasets...

NSIDC_contour_distance_loop <- function(df, path, contour){
  
  ## build a file-db
  files <- tibble(fullname = list.files(normalizePath(path), full.names = TRUE, pattern = "tif$"), 
                  date = as.POSIXct(strptime(basename(fullname), "N_%Y%m%d", tz = "UTC")))
  
  ## map our data to this file-db (we might have multiple points per file)
  df$fullname <- files$fullname[findInterval(df$date, files$date)]
  
  ## load a dummy raster to pull projection info from
  (rdummy <- raster(df$fullname[1]))
  
  ## project the query points
  df[c("X", "Y")] <- as_tibble(rgdal::project(as.matrix(df[c("lon", "lat")]), projection(rdummy)))
  
  ## set up progress bar based on number of files that will be loaded
  #  pb <- progress_estimated(length(unique(df$fullname)))
  pb1 <- txtProgressBar(min = 1, max = length(unique(df$fullname)), style = 3) 
  
  for (i in 1:length(unique(df$fullname))){
    setTxtProgressBar(pb1, i) # update progress bar
    r <- raster(unique(df$fullname)[i])
    r[r == 2510] <- 1000 # make pole hole 100% ice cover
    r[r > 1000] <- NA
    r <- r/10
    c <- rasterToContour(r, levels = contour) %>% st_as_sf()
    
    df$dist[df$fullname == unique(df$fullname)[i]] <- df %>% filter(fullname == unique(fullname)[i]) %>%
      st_as_sf(coords = c("X", "Y")) %>% st_set_crs(projection(rdummy)) %>% st_distance(c,  by_element = T) %>% as.vector()
  }
  
  return(df)
}




