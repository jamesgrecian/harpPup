########################################################
### Calculate 6 hour dive summary times for foieGras ###
########################################################

# Given a dataframe of location data from a SMRU SRDL
# Calculate the 6 hour intervals so that foieGras interpolated location match dive summary data

SRDL_time <- function(locs){
  
  time <- locs %>%
    group_by(id) %>%
    summarise(min = floor_date(min(locs$date), "6 hours"),
              max = ceiling_date(max(locs$date), "6 hours")) %>%
    mutate(date = map2(min, max, seq, by = "6 hours")) %>% # Create a list column with dates
    unnest(cols = c(date)) %>%
    select(id, date)
  
  return(time)
}
