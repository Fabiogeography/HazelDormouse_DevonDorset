##############################################################
##### PROCESS ENVIRONMENTAL DATA - add hedgerow data #############################
##### Written by: Regan Early ################################
##### Written on: 30th November 2021 #############################
##############################################################

# .libPaths("C:/Rpackages") 
# .libPaths("D:\\SOFTWARE\\R-4.1.1\\library")
library(raster)#, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(rgdal)#, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(pROC)#, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(tidyverse)#, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library") ## glimpse and useful functions
library(corrplot)#, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(Hmisc)#, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library") ## rcorr
library(ggplot2)#, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(fasterize)#, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")

wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
wd.distn <- "E:/NON_PROJECT/DORMOUSE_NE/DISTRIBUTION"
wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/EXPLORE"
wd.clim <- "E:/NON_PROJECT/DORMOUSE_WALES/GIS/100m" # "E:/GIS_DATA/CLIMATE/UK/HadUK-Grid_1km"
OSGB.proj <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'

##### Environmental data created by code1_process_env_data #####
env.r <- stack(paste0(wd.env, "/env_100m7"))

## crop and mask to Devon Dorset area
dd <- readOGR(dsn="E:/NON_PROJECT/DORMOUSE_NE/GIS", layer="devondorset_buffer5km", stringsAsFactors = F)
env.dd <- crop(env.r, extent(dd))
env.dd <- raster::mask(env.dd, dd) ## mask is also in Hmisc
writeRaster(env.dd, paste0(wd.env, "/env_100m7_dd"))
env.dd <- stack(paste0(wd.env, "/env_100m7_dd"))

##### Hedgerow data #####
thaw <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/Devon_THaW_11_11_21.tif")

## Trim to same extent as env data
test <- crop(env.r, thaw)
thaw <- crop(thaw, test) ## env.r

## Classify categories 2-4 and 6 as trees / hedges (see phase 2 notes for working on this)
m <- matrix(c(0,1.5,0, 1.5,4.5,1, 4.5,5.5,0,  5.5,6.5,1), nrow=4, byrow=T)
thaw2 <- reclassify(thaw, m)
# table(getValues(thaw2, row=5000)) ## check numbers add up
# table(getValues(thaw, row=5000)) ## check numbers add up
# writeRaster(thaw2, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/trees_hedge_1m.tif")
thaw2 <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/trees_hedge_1m.tif")

## Calculate the sum of 1m tree / hedge cells in each 100m grid-cell
thaw100 <- aggregate(thaw2, fact=100, fun=sum)/1000 ## onvert to percentage
writeRaster(thaw100, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/trees_hedge_100m.tif")
thaw100 <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/trees_hedge_100m.tif")
thaw100 <- projectRaster(thaw100, crs=OSGB.proj)

## Exclude the NFI and LCM forest areas (created in ArcMap map3_phase2) from the 100m aggregated THAW data. 
forest <- readOGR(dsn="E:/NON_PROJECT/DORMOUSE_NE/GIS", layer="NFI_LCM_broadl_conif_mix", stringsAsFactors = F) ## 2015 data
hedge100 <- raster::mask(thaw100, forest, inverse=T) ## mask is also in Hmisc
### This approach seems to remove raster grid-cells that are completely contained by the polygon (see ArcMAP)
### Therefore grid-cells that overlap forest polygons substantially will have high NFI/LCM forest cover AND high hedge cover.
### Could lead to weird predictions at forest edges?
### Perhaps instead substract forest cover at 100m resolution from hedge?
## Make a raster with coverage of LCM and NFI forest in 2015 
forest.r <- raster::rasterize(forest, thaw100, getCover=T) ## also a function in terra. ## Can't use fasterize as doesn't have function for getting cover.
writeRaster(forest.r, "E:/NON_PROJECT/DORMOUSE_NE/GIS/NFI_LCM_broadl_conif_mix.tif")

hedge100 <- thaw100 - forest.r
writeRaster(hedge100, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/hedge_100m.tif")

## Exclude the NFI and LCM forest areas (created in ArcMap map3_phase2) from the 1m THAW data. 
hedge <- raster::mask(thaw, forest, inverse=T) ## mask is also in Hmisc
writeRaster(hedge, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/hedge.tif")




# hedge100 <- aggregate(hedge, fact=100)
# writeRaster(hedge100, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/hedge_100m.tif")


