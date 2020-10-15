library(rcarbon)
library(parallel)
ncores = (detectCores() - 1)

rcarbon.spd <- function(c14.data, 
                        timeRange,
                        kernel,
                        oxcalnorm = FALSE, 
                        median = FALSE, 
                        spdnormalised =  FALSE
                        ){
  
  cal = rcarbon::calibrate(x = c14.data$C14AGE,
                           errors = c14.data$C14STD,
                           calCurves = 'intcal13', 
                           ncores = ncores, 
                           normalised = TRUE) #running calibration over 3 cores
  
  spd <- rcarbon::spd(cal,
                      timeRange = timeRange, 
                      spdnormalised = spdnormalised, 
                      runm = kernel)
  
  spd <- as.data.frame(spd[2])
  
  if(oxcalnorm) {
    # raise to max() == 1 like OxCal does!
    spd$grid.PrDens <- spd$grid.PrDens/max(spd$grid.PrDens, 
                                           na.rm = TRUE)
  }
  
  spd <- spd[spd$grid.PrDens != 0,] # remove 0 values
  spd$calBCAD <- 1950 - spd$grid.calBP

  start <- min(spd$grid.calBP)
  
  if(median){
    med <- list()
    for(k in 1:length(cal$grids)){
      m <- 1950 - median(cal$grids[[k]]$calBP, na.rm = T)
      med[k] <- m
    }
    median <- do.call(rbind, med)
    
    output <- list(dat = spd, 
                   median = median, 
                   start = start)
  }else{
    output <- list(dat = spd, 
                   start = start)
  }
}

basemap <- function(){
  
  white <- sf::st_read("input/whitesveg/Whites vegetation.shp") %>%
    st_set_crs(4326) %>%
    dplyr::filter(DESCRIPTIO %in% c("Anthropic landscapes",
                                    "Dry forest and thicket",
                                    "Swamp forest and mangrove",
                                    "Tropical lowland rainforest"))  

  # Vector layers ----
  rivers10 <- ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical", returnclass="sf")
  lakes10 <- ne_download(scale = 10, type = "lakes", category = "physical", returnclass="sf")
  coast10 <- ne_download(scale = 10, type = "coastline", category = "physical", returnclass="sf")
  land10 <- ne_download(scale = 10, type = "land", category = "physical", returnclass="sf")
  boundary_lines_land10 <- ne_download(scale = 10, type = "boundary_lines_land", category = "cultural", returnclass="sf")
  
  plt <- ggplot() + 
    # base map ----
    #geom_raster(data = rfs.rainforst, aes(y = y, x = x), fill = '#c8c8c8') + 
    #geom_raster(data = rfs.water, aes(y = y, x = x), fill = '#808080') + 
    geom_sf(data = white, fill = "grey", color = NA) + 
    geom_sf(data = coast10, size = .5, color = '#808080') + 
    geom_sf(data = rivers10, size = .5, color = '#808080') + 
    geom_sf(data = lakes10, fill = '#808080', color = NA) + 
    geom_sf(data = boundary_lines_land10, size = .1, color = 'black')# + 
    #coord_sf(xlim = c(7.5, 29), 
    #         ylim = c(-9.2, 6))

  return(plt)
}