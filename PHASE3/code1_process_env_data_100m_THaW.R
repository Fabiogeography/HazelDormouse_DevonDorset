##############################################################
##### PROCESS ENVIRONMENTAL DATA - add hedgerow data #############################
##### Written by: Regan Early ################################
##### Written on: 30th November 2021 #############################
##### Modified on: 8th Feb 2022 #############################
##############################################################

.libPaths("D:\\SOFTWARE\\R-4.1.2\\library")
library(raster)
library(rgdal)
library(pROC)
library(tidyverse)# ## glimpse and useful functions
library(corrplot)
library(Hmisc) ## rcorr
library(ggplot2)
library(fasterize)

wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
wd.distn <- "E:/NON_PROJECT/DORMOUSE_NE/DISTRIBUTION"
wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/EXPLORE"
wd.clim <- "E:/NON_PROJECT/DORMOUSE_WALES/GIS/100m" # "E:/GIS_DATA/CLIMATE/UK/HadUK-Grid_1km"
OSGB.proj <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'

##### Region data #####
## Make 100m mask raster
t <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/Devon_THaW_11_11_21.tif")
# t100 <- aggregate(t, fact=100, fun=sum) ## convert to 100m resolution. Lengthy. Doesn't match grid100m extent.
# m <- matrix(c(0,10,1), ncol=3, byrow=T)
# t100 <- reclassify(t100, m) ## Make all values 1
# t100 <- resample(t100, grid100m) ## ensure aligns with base grid
# writeRaster(t100, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/thaw_region.tif")
t100 <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/thaw_region.tif")

t <- crop(t, t100) ## Crop 1m thaw data to extent of t100, which means that aggregate functions below give a 100m raster that aligns with the UK grid. Lengthy.
writeRaster(t, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/Devon_THaW_11_11_21crop.tif")

##### Hedgerow data #####
## Classify categories 2-4 and 6 as trees / hedges (see phase 2 notes for working on this)
m <- matrix(c(0,1.5,0, 1.5,4.5,1, 4.5,5.5,0,  5.5,6.5,1), nrow=4, byrow=T)
treehedge <- reclassify(t, m)
# table(getValues(thaw2, row=5000)) ## check numbers add up
# table(getValues(thaw, row=5000)) ## check numbers add up
writeRaster(treehedge, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/trees_hedge_1m.tif")
# treehedge <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/trees_hedge_1m.tif")

## Calculate the sum of 1m tree / hedge cells in each 100m grid-cell.
# thaw3 <- crop(thaw2, t100) ## Crop thaw2 to extent of t100, which will hopefully mean that aggregate gives a 100m raster that aligns with the UK grid 
treehedge100 <- aggregate(treehedge, fact=100, fun=sum)/10000 ## convert to proportion. This function is much faster in ArcMap.
writeRaster(treehedge100, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/trees_hedge_100m.tif")
# treehedge100 <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/trees_hedge_100m.tif")

## Exclude the NFI and LCM forest areas (created in ArcMap map3_phase2) from the 100m aggregated THAW data. 
forest <- readOGR(dsn="E:/NON_PROJECT/DORMOUSE_NE/GIS", layer="NFI_LCM_broadl_conif_mix", stringsAsFactors = F) ## 2015 data
## Not run:
# hedge100 <- raster::mask(treehedge100, forest, inverse=T) ## mask is also in Hmisc
## This approach seems to remove raster grid-cells that are completely contained by the polygon (see ArcMAP)
## Therefore grid-cells that overlap forest polygons substantially will have high NFI/LCM forest cover AND high hedge cover.
## Could lead to weird predictions at forest edges.
## Instead subtract forest cover at 100m resolution from hedge.

### Make a raster with coverage of LCM and NFI forest in 2015 
## Not run:
# forest.r <- raster::rasterize(forest, treehedge100, getCover=T) ## also a function in terra. ## Can't use fasterize as doesn't have function for getting cover.
# writeRaster(forest.r, "E:/NON_PROJECT/DORMOUSE_NE/GIS/NFI_LCM_broadl_conif_mix.tif")
## The above approach fails. Made in ArcMap instead, so load up here:
# forest.r100 <- raster("E:/NON_PROJECT/DORMOUSE_NE/GIS/nfilcm_bcm10b") ## Converted NFI_LCM_broadl_conif_mix.shp to 1m raster and used aggregate to sum number of 1m treehedge cells in each 200m cell. This function is much faster in ArcMap. But ArcMap raster calculator fails, so calculate proportion here., 
# forest.r100 <- forest.r100/10000 ## convert count to proportion
# forest.r100 <- projectRaster(forest.r100, crs=OSGB.proj) ## Doesn't change projection,m just ensures defined the same as other rasters 
# writeRaster(forest.r100, "E:/NON_PROJECT/DORMOUSE_NE/GIS/nfilcm_bcm100.tif")
forest.r100 <- raster("E:/NON_PROJECT/DORMOUSE_NE/GIS/nfilcm_bcm100.tif")

hedge100 <- treehedge100 - forest.r100
writeRaster(hedge100, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/hedge_100m.tif")

##### Scrub data #####
## Classify categories 1 and 5 as scrub below 1.3m (see phase 2 notes for working on this)
m <- matrix(c(0.5,1.5,1,  1.5,4.5,0,  4.5,5.5,1,  5.5,6.5,0), nrow=4, byrow=T)
scrub <- reclassify(t, m)
# table(getValues(t, row=5000)) ## check numbers add up
# table(getValues(scrub, row=5000)) ## check numbers add up
writeRaster(scrub, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/scrub_1m.tif")
# scrub <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/scrub_1m.tif")

## Modify the sum of 1m scrub cells in each 100m grid-cell. THe aggregate was done in ArcMap as that is much faster.
scrub100 <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/scrub_100mX.tif")
scrub100 <- scrub100/10000 ## convert sum of 1m cells to a proportion
m <- matrix(c(NA,NA,0), nrow=1, byrow=T) ## Convert NAs to 0s.
scrub100 <- reclassify(scrub100, m)
scrub100 <- scrub100*t100 ## Remove any cells that should actually be NA
plot(scrub100)
writeRaster(scrub100, file="E:/GIS_DATA/LandCover/UK/Devon_THaW/scrub_100m.tif")


