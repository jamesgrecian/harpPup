####################################
### Download  NSIDC sea ice data ###
####################################

# Given a vector of dates and a path to write files to, download NSIDC daily data and save geotiffs
# Should be possible to use purrr::map with a progress bar...?!

get_NSIDC <- function(dates, path){
  
  dates <- unique(date(dates)) # strip the unique dates from the datetime vector
  
  # Strip date components and construct ftp path to file
  yr <- year(dates)
  mo <- month(dates)
  mo <- ifelse(mo < 10, paste("0", mo, sep=""), mo)
  mon <- month(dates, label = T, abbr = T)
  dy = day(dates)
  dy <- ifelse(dy < 10, paste("0", dy, sep=""), dy)
  fn <- paste0("ftp://anonymous:wjg5@sidads.colorado.edu/DATASETS/NOAA/G02135/north/daily/geotiff/",
               yr, "/", mo, "_", mon, "/N_", yr, mo, dy, "_concentration_v3.0.tif")
  
  dir.create(path, showWarnings = F)
  
  ## get all the files, though avoid re-downloading
  download_it <- function(x) {
    file <- file.path(path, basename(x))
    if (!file.exists(file)) {
      curl::curl_download(x, file)
      return(TRUE)
    } 
    FALSE
  } 
  
  pb1 <- txtProgressBar(min = 1, max = length(dates), style = 3) 
  
  #purrr::map_lgl(fn, download_it)
  for (i in 1:length(dates)){
    setTxtProgressBar(pb1, i) # update progress bar
    download_it(fn[i])
  }
  
}

