##### MAKE PRESENCE AND ABSENCE SHAPEFILES ########################
##### Written by: Regan Early ################################
##### Written on: 7th February 2022 ##############################
##### Modified on: 7th February 2022  #########################
##############################################################

.libPaths("C:\\SOFTWARE\\R-4.1.2\\library")
library(rgdal)
library(raster)
library(sp) ## over
library(FNN) ## get.knnx

wd.distn <- "E:\\NON_PROJECT\\DORMOUSE_NE\\DISTRIBUTION"
wd.env <- "E:\\NON_PROJECT\\DORMOUSE_NE\\GIS\\100m"

OSGB.proj <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
wgs.proj <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

### Thaw region
t100 <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/thaw_region.tif")

##### PRESENCES #####
p100.1990 <- readOGR(dsn=wd.distn, layer="p100m1990_dd") ## Already gridded and filtered to 1km, but needs to be trimmed to THaW area

## Obtain 100m grid
grid100m <- raster("E:\\GIS_DATA\\UK\\OSGB_GRIDS\\OSGB_Grid_100m_clim.tif") ## This was made in ArcMap by converting the original raster to a points shapefile, trimming to the 1km_grid_region_clim polygon shapefile, and converting back to raster. Process in R didn't produce unique cell ids for each 100m grid-cell

## Trim presences to THaW area
p100.1990$t <- raster::extract(t100, p100.1990)
p100t <- p100.1990[!is.na(p100.1990$t),]
p100t$t <- NULL

## Save
writeOGR(p100t, dsn=wd.distn, layer="p100m1990_t", driver="ESRI Shapefile")

##### ABSENCES #####
# a <- readOGR(dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="XYsurv_ab_24Jan2022")
a <- read.csv(paste0(wd.distn,"\\SURVEY_ABSENCES\\surv_ab_24Jan2022.csv"))
a$x <- a$easting ## make new coordinate columns so that easting and northing can remain in @data
a$y <- a$northing
coordinates(a) <- ~x+y
proj4string(a) <- OSGB.proj
writeOGR(a, dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="XYsurv_ab_24Jan2022", driver="ESRI Shapefile")

### Trim down columns in presence data. See metadata file for column definitions. SiteName and Source may be useful down the line. 
a <- a[,c("site","easting", "northing", "year", "ptes", "sei")]
colnames(a@data) <- c("site", "easting","northing", "year", "source", "sei")

### Trim to THaW region
a100 <- a
a100$m100 <- raster::extract(t100, a100)
summary(a100$m100) ## 246 fall outside THaW region
a100 <- a100[!is.na(a100$m100),] ## remove points that fall outside env data. 335 remaining records

### Filter to point point per 100m grid-cell. Keep the row with the most recent date (some cells have > 1 record on the same date, and one of these is chosen randomly)
## Don't trim to 1990 since all points are from 2002-2021.
a100$m100 <- raster::extract(grid100m, a100) ## assign unique value to points in each 1km grid-cell
a100 <- a100[!is.na(a100$m100),] ## remove points that fall outside env data. 1 record removed
a100 <- a100[order(a100$year, decreasing=TRUE),] 
a100 <- a100[!duplicated(a100$m100),] ## 309 in seperate 1km grid-cells records. Nice.

a100$occ <- 0
writeOGR(a100, dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="a100m_t", driver="ESRI Shapefile")


##### Make occupancy file (in PHAS2 this was created in code 3 - explore) #####
p <- readOGR(dsn=wd.distn, layer="p100m1990_t")
a <- readOGR(dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="a100m_t")

proj4string(p) <- proj4string(a) <- OSGB.proj ## Already in OSGB so doesn't transform, just ensures projection is written identically to other objects

## Change name of a column so it matches columns in p
a$IDorig <- a$site
a$site <- NULL

## Add extra columns to p
x <- names(a)[!(names(a) %in% names(p))] ## Add columns in p to a
p@data[,x] <- NA

occ <- rbind(p, a)

##### Remove points that are dubious or could cause problems for modelling #####
### Remove presences and absences that fall in the same 100m cell. Note that these presences are not from NDMP sites, so wouldn't be picked up in the next step.
dups <- which(duplicated(occ$m100)) ## Identify duplicates. Only finds the second row, not the first
dups <- occ$m100[dups] ## duplicated m100 (grid-cell) values
occ <- occ[which(!(occ$m100 %in% dups)),] ## Remove the rows where presence and absence points fall in the same grid-cell

### Remove absences within 1km of an NDMP presence
ndmp <- occ[occ$source=="NDMP",]
survAbs <- occ[occ$occ==0,]

## get the distances
g <- get.knnx(coordinates(ndmp), coordinates(survAbs),k=1) ## Inputs must be in planar coordinates
summary(g$nn.dist) ## The distances between each absence point and the nearest NDMP presence point.
tooClose <- occ$m100[g$nn.dist<=1000]

## Remake occ
occ <- occ[!(occ$m100 %in% tooClose & occ$occ==0),] ## Removes 8 points. This is fewer than the 12 mentioned in the notes, but this is because here I'm removing 100m grid-cells, whereas in the notes I was referring to points, and there may be multiple points within each grid-cell.


writeOGR(occ, dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="occ100mt", driver="ESRI Shapefile")
