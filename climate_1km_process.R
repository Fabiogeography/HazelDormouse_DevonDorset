##############################################################
##### PROCESS CLIMATE DATA AT 1KM  ###########################
##### Monthly data from: http://data.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.2.1/1km
##### Written by: Regan Early ################################
##### Written on: 15th Feb 2021 ##############################
##### Modified on:   #########################################
##############################################################

##### NOTE on extent #####
## climate might not be quite aligned with the 1km OSGB grid which I use for a lot of other analyses
# grid <- readOGR(dsn="E:\\GIS_DATA\\UK\\OSGB_GRIDS", layer="1km_grid_region")
# proj4string(grid) <- OSGB.proj
# 
# ## To fix:
# nfi.c <- stack(paste0(wd.env,"\\",nfi.layers, "_change.tif")) ## Change in habitat cover between 2011 and 2017 inclusive
# test <- crop(clim, c(0, extent(nfi.c)[2], 0, extent(nfi.c)[4]))
# extent(test) <- extent(nfi.c)
# test2 <- stack(nfi.c, test)
#####

library(raster)
wd.env <- "E:\\GIS_DATA\\CLIMATE\\UK\\HadUK-Grid_1km\\"
OSGB.proj <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'

months <- c("jan","feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

##### Rainfall variables (sum) #####
for (y in 1986:2019) {

b <- stack(paste0(wd.env, "rainfall_monthly\\rainfall_hadukgrid_uk_1km_mon_",y,"01-",y,"12.nc"))
names(b) <- months 
proj4string(b) <- OSGB.proj

spr <- sum(b[[c("mar","apr","may")]]) ## March - May
names(spr) <- paste0("rainfall_spr_",y)
writeRaster(spr, paste0(wd.env, "rainfall_monthly\\rainfall_spr_",y,".tif"))

smr <- sum(b[[c("jun","jul","aug")]]) ## June - August
names(smr) <- paste0("rainfall_smr_",y)
writeRaster(smr, paste0(wd.env, "rainfall_monthly\\rainfall_smr_",y,".tif"))

aut <- sum(b[[c("sep","oct","nov")]]) ## September - November 
names(aut) <- paste0("rainfall_aut_",y)
writeRaster(aut, paste0(wd.env, "rainfall_monthly\\rainfall_aut_",y,".tif"))

}

##### Sunshine hours variables (mean) #####
for (y in 1986:2019) {
  
  b <- stack(paste0(wd.env, "sun_monthly\\sun_hadukgrid_uk_1km_mon_",y,"01-",y,"12.nc"))
  names(b) <- months 
  proj4string(b) <- OSGB.proj
  
  spr <- mean(b[[c("mar","apr","may")]]) ## March - May
  names(spr) <- paste0("sun_spr_",y)
  writeRaster(spr, paste0(wd.env, "sun_monthly\\sun_spr_",y,".tif"))
  
  smr <- mean(b[[c("jun","jul","aug")]]) ## June - August
  names(smr) <- paste0("sun_smr_",y)
  writeRaster(smr, paste0(wd.env, "sun_monthly\\sun_smr_",y,".tif"))
  
  aut <- mean(b[[c("sep","oct","nov")]]) ## September - November 
  names(aut) <- paste0("sun_aut_",y)
  writeRaster(aut, paste0(wd.env, "sun_monthly\\sun_aut_",y,".tif"))
  
}


##### Minimum temperature variables (mean) #####
for (y in 1986:2019) {
  
  b <- stack(paste0(wd.env, "tasmin_monthly\\tasmin_hadukgrid_uk_1km_mon_",y,"01-",y,"12.nc"))
  names(b) <- months 
  proj4string(b) <- OSGB.proj
  
  spr <- mean(b[[c("mar","apr","may")]]) ## March - May
  names(spr) <- paste0("tasmin_spr_",y)
  writeRaster(spr, paste0(wd.env, "tasmin_monthly\\tasmin_spr_",y,".tif"))

  smr <- mean(b[[c("jun","jul","aug")]]) ## June - August
  names(smr) <- paste0("tasmin_smr_",y)
  writeRaster(smr, paste0(wd.env, "tasmin_monthly\\tasmin_smr_",y,".tif"))

  aut <- mean(b[[c("sep","oct","nov")]]) ## September - November
  names(aut) <- paste0("tasmin_aut_",y)
  writeRaster(aut, paste0(wd.env, "tasmin_monthly\\tasmin_aut_",y,".tif"))
  
  win <- mean(b[[c("dec","jan","feb")]]) ## December - February
  names(win) <- paste0("tasmin_win_",y)
  writeRaster(win, paste0(wd.env, "tasmin_monthly\\tasmin_win_",y,".tif"))
  
}

##### Maximum temperature variables (mean) #####
for (y in 1986:2019) {
  
  b <- stack(paste0(wd.env, "tasmax_monthly\\tasmax_hadukgrid_uk_1km_mon_",y,"01-",y,"12.nc"))
  names(b) <- months 
  proj4string(b) <- OSGB.proj
  
  spr <- mean(b[[c("mar","apr","may")]]) ## March - May
  names(spr) <- paste0("tasmax_spr_",y)
  writeRaster(spr, paste0(wd.env, "tasmax_monthly\\tasmax_spr_",y,".tif"))
  
  smr <- mean(b[[c("jun","jul","aug")]]) ## June - August
  names(smr) <- paste0("tasmax_smr_",y)
  writeRaster(smr, paste0(wd.env, "tasmax_monthly\\tasmax_smr_",y,".tif"))
  
  aut <- mean(b[[c("sep","oct","nov")]]) ## September - November 
  names(aut) <- paste0("tasmax_aut_",y)
  writeRaster(aut, paste0(wd.env, "tasmax_monthly\\tasmax_aut_",y,".tif"))
  
}

##### Winter temperature range #####
for (y in 1986:2019) {
  
  b.max <- stack(paste0(wd.env, "tasmax_monthly\\tasmax_hadukgrid_uk_1km_mon_",y,"01-",y,"12.nc"))
  names(b.max) <- months 
  proj4string(b.max) <- OSGB.proj
  win.max <- max(b.max[[c("dec","jan","feb")]]) ## Maximum mean daily temperature in December - February 

  b.min <- stack(paste0(wd.env, "tasmin_monthly\\tasmin_hadukgrid_uk_1km_mon_",y,"01-",y,"12.nc"))
  names(b.min) <- months 
  proj4string(b.min) <- OSGB.proj
  win.min <- min(b.min[[c("dec","jan","feb")]]) ## Minimum mean daily temperature in December - February 

  win.rng <- win.max - win.min
  
  names(win.rng) <- paste0("tasrng_win_",y)
  writeRaster(win.rng, paste0(wd.env, "\\tasrng_win_",y,".tif"))
  
}

##### Calculate 5 year means #####

for (i in 1990:2019) {

  ### Rainfall
  s <- mean(stack(paste0(wd.env, "rainfall_monthly\\rainfall_spr_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "rainfall_monthly\\rainfall_spr_5yrmn_",i,".tif"))

  s <- mean(stack(paste0(wd.env, "rainfall_monthly\\rainfall_smr_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "rainfall_monthly\\rainfall_smr_5yrmn_",i,".tif"))

  s <- mean(stack(paste0(wd.env, "rainfall_monthly\\rainfall_aut_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "rainfall_monthly\\rainfall_aut_5yrmn_",i,".tif"))

  ### Sunshine hours
  s <- mean(stack(paste0(wd.env, "sun_monthly\\sun_spr_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "sun_monthly\\sun_spr_5yrmn_",i,".tif"))

  s <- mean(stack(paste0(wd.env, "sun_monthly\\sun_smr_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "sun_monthly\\sun_smr_5yrmn_",i,".tif"))

  s <- mean(stack(paste0(wd.env, "sun_monthly\\sun_aut_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "sun_monthly\\sun_aut_5yrmn_",i,".tif"))

  ### Minimum temperature variables (mean)
  s <- mean(stack(paste0(wd.env, "tasmin_monthly\\tasmin_spr_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "tasmin_monthly\\tasmin_spr_5yrmn_",i,".tif"))

  s <- mean(stack(paste0(wd.env, "tasmin_monthly\\tasmin_smr_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "tasmin_monthly\\tasmin_smr_5yrmn_",i,".tif"))

  s <- mean(stack(paste0(wd.env, "tasmin_monthly\\tasmin_aut_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "tasmin_monthly\\tasmin_aut_5yrmn_",i,".tif"))

  ### Maximum temperature variables (mean)
  s <- mean(stack(paste0(wd.env, "tasmax_monthly\\tasmax_spr_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "tasmax_monthly\\tasmax_spr_5yrmn_",i,".tif"))

  s <- mean(stack(paste0(wd.env, "tasmax_monthly\\tasmax_smr_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "tasmax_monthly\\tasmax_smr_5yrmn_",i,".tif"))

  s <- mean(stack(paste0(wd.env, "tasmax_monthly\\tasmax_aut_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "tasmax_monthly\\tasmax_aut_5yrmn_",i,".tif"))
  
  ### Winter temperature range 
  s <- mean(stack(paste0(wd.env, "tasrng_monthly\\tasrng_win_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "tasrng_monthly\\tasrng_win_5yrmn_",i,".tif"))
  
  s <- mean(stack(paste0(wd.env, "tasmin_monthly\\tasmin_win_", i-c(0:4),".tif")))
  writeRaster(s, paste0(wd.env, "tasmin_monthly\\tasmin_win_5yrmn_",i,".tif"))
  
}
