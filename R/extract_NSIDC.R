#############################################
### Extract  NSIDC sea ice data to points ###
#############################################

# Given a location dataframe and a path to NSIDC geotiffs
# Import NSIDC daily data and extract sea ice concentration to specified locations

extract_NSIDC <- function(df, path){

  ## build a file-db
  files <- tibble(fullname = list.files(normalizePath(path), full.names = TRUE, pattern = "tif$"), 
                  date = as.POSIXct(strptime(basename(fullname), "N_%Y%m%d", tz = "UTC")))
  
  ## map our data to this file-db (we might have multiple points per file)
  df$fullname <- files$fullname[findInterval(df$date, files$date)]
  
  ## set up progress bar based on number of files that will be loaded
  pb <- progress_estimated(length(unique(df$fullname)))
  
  (rdummy <- raster(df$fullname[1]))
  
  ## project the query points
  df[c("X", "Y")] <- as_tibble(rgdal::project(as.matrix(df[c("lon", "lat")]), projection(rdummy)))
  
  ## now, extract per file 
  df <- purrr::map_df(split(df, df$fullname)[unique(df$fullname)], 
                      function(.x) {
                        pb$tick()$print()
                        .x["concentration"] <- raster::extract(raster(.x$fullname[1]), as.matrix(.x[c("X", "Y")]))
                        .x
                      })
  # if animal locations are further south than ice raster extent set concentration to 0
  df <- df %>% mutate(concentration = if_else(Y < raster::extent(r)[3], 0, concentration))

  # #0 is ocean; 2510 pole hole; 2530 coast line; 2540 land; 2550 missing
  # 0-1000 so divide by 10 to get percentage
  df$concentration[df$concentration > 1000] <- NA
  df$concentration <- df$concentration/10
  
  # check order of dataframe
  df <- df %>% arrange(id, date)
  df <- df %>% dplyr::select(-c("fullname", "X", "Y"))
  
  return(df)
}