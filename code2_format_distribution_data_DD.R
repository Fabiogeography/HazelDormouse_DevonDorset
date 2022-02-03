##### MAKE PRESENCE AND ABSENCE SHAPEFILES ########################
##### Written by: Regan Early ################################
##### Written on: 7th Jan 2021 ##############################
##### Modified on: 1st November 2021  #########################
##############################################################

library(rgdal, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(raster, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(sp, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library") ## over

wd.distn <- "E:\\NON_PROJECT\\DORMOUSE_NE\\DISTRIBUTION"
wd.env <- "E:\\NON_PROJECT\\DORMOUSE_NE\\GIS\\100m"

OSGB.proj <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
wgs.proj <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

## Devon Dorset area
dd <- readOGR(dsn="E:/NON_PROJECT/DORMOUSE_NE/GIS", layer="devondorset_buffer5km", stringsAsFactors = F)

##### PRESENCES #####
##### PTES records #####
### 1988-2019. See metadata doc in same folder. 
ptes <- readOGR(dsn=paste0(wd.distn, "\\PTES"), layer="XYNDD_2011_2019_DevonDorset")

## Use all records with Precision 100, 10, or 1m as when only 'good' records were used this was 17170 of 17365 records. I removed 3 records with 0 precision as I wasn't sure what this meant.
ptes <- ptes[ptes$Precision %in% c(100,10,1),]

## Use all records with RecordTypeReliability 'Good' as this is 17232 of all 17365 records and have either been collected by trained surveyors or verified by PTES.
ptes <- ptes[ptes$RecordTy_1=="Good",]

## Remove invalid coordinates
invalid <- ptes@coords[,"coords.x1"] > 660000 ## One coordinates are invalid
ptes <- ptes[invalid==F,]

## extract year from StartDate (checked that this is always the same year as EndDate)
ptes$year <- as.numeric(paste0("20", substr(ptes$StartDate, 8,9)))

## Make all locations that aren't recorded as NDMP site or not 'No'
ptes$NDMPSite[is.na(ptes$NDMPSite)] <- "No"

## Trim down columns in presence data. See metadata file for column definitions. SiteName and Source may be useful down the line. 
ptes <- ptes[,c("RecordKey","Easting", "Northing", "year", "Source")]
colnames(ptes@data) <- c("IDorig", "easting","northing", "year", "source")
ptes$source[ptes$source %in% c("LRD", "NDD")] <- "LRD+NDD"

# env.p <- cbind(ptes[c("Source", "NDMPSite", "RecordType")], raster::extract(env.dd, ptes))
# env.p <- distinct(env.p@data) ## 792 records. Occ dataset indicates only 623 100m grid-cells occupied, so some of the occurrences will be duplicates. 

##### Grid and filter at 100m #####
## Obtain 100m grid and trim to match the environmental data more closely
grid100m <- raster("E:\\GIS_DATA\\UK\\OSGB_GRIDS\\OSGB_Grid_100m_clim.tif") ## This wsa made in ArcMap by converting the original raster to a points shapefile, trimming to the 1km_grid_region_clim polygon shapefile, and converting back to raster. Process in R didn't produce unique cell ids for each 100m grid-cell

## crop and mask to Devon Dorset area
grid100m <- crop(grid100m, extent(dd))
grid100m <- raster::mask(grid100m, dd) ## mask is also in Hmisc

p100 <- ptes
p100$m100 <- raster::extract(grid100m, p100)
summary(p100$m100) ## 24 points fall outside the environmental data. Too few to worry about. 
p100 <- na.omit(p100) ## remove points that fall outside env data

## How many grid-cells have both NDMP and NRD+LDD records?
test <- na.omit(p100@data[,c("source", "m100")]) ## 3072 rows
test <- distinct(test) ## 624 rows
length(table(test$m100)) ## 620, so 4 grid-cells have records from both NDMP and LRD+NDD sources. Too few to worry about. 
tail(sort(table(test$m100))) ## Identify the grid-cells 

### Filter to point point per 100m grid-cell. Keep the row with the most recent date (some cells have > 1 record on the same date, and one of these is chosen randomly)
p100 <- p100[order(p100$year, decreasing=TRUE),] ## 3096 records
p100 <- p100[!duplicated(p100$m100),] ## 621 records

p100$occ <- 1
writeOGR(p100, dsn=wd.distn, layer="p100m_dd", driver="ESRI Shapefile")
# p100 <- readOGR(dsn=wd.distn, layer="p100m_dd")

### Filter time periods
hist(p100$year)
hist(p100$year[p100$year>1990])
hist(p100$year[p100$year>1990], breaks=c(1990:2021))

p100.1990 <- p100[p100$year %in% c(1990:2020),]
writeOGR(p100.1990, dsn=wd.distn, layer="p100m1990_dd", driver="ESRI Shapefile")

##### ABSENCES #####
##### See dormouse Wales project for investigation of various absence options.
n <- nrow(p100.1990) ## Number of absences to select
# e <- extent(dd) ## geographic extent in which to place the absences
# e@ymax <- 396000 ## trim it a bit, as there are just six presences more northerly than the cutoff
# test <- readOGR(dsn="E:\\GIS_DATA\\UK\\OSGB_GRIDS", layer="1km_grid_region_WalesEngland")


### Randomly selected to match 100m records from 1990 onwards. 

a100 <- as.data.frame(sampleRandom(grid100m, n*2, xy=T)) ## select some random grid-cells, more than needed
match <- a100$OSGB_Grid_100m %in% p100.1990$m100 ## grid-cells where presences are recorded. OSGB_Grid_100m
a100 <- a100[match==F,]
a100 <- a100[sample(c(1:nrow(a100)), n),] ## select correct number of grid-cells
coordinates(a100) <- ~x+y ## Obtain the centroids of the 100m cells where absences are placed. coords
a100@data$occ <- 0
proj4string(a100) <- OSGB.proj
writeOGR(a100, dsn=wd.distn, layer="a100m1990_dd", driver="ESRI Shapefile")

### 100m absences within 5km buffer, randomly selected to match 100m records from 1990 onwards.
## Note probably not much use as covers most of dd area

buf <- buffer(p100, width=5000) ## Units are in m so 5,000 = 5km
buf <- spTransform(buf, OSGB.proj)
grid.buf <- mask(grid100m, buf)
a100b <- as.data.frame(sampleRandom(grid.buf, n*2, xy=T, ext=e)) ## select some random grid-cells, more than needed
match <- a100b$OSGB_Grid_100m %in% p100.1990$m100 ## grid-cells where presences are recorded. OSGB_Grid_100m
a100b <- a100b[match==F,]
a100b <- a100b[sample(c(1:nrow(a100)), n),] ## select correct number of grid-cells
coordinates(a100b) <- ~x+y ## Obtain the centroids of the 100m cells where absences are placed. coords
a100b@data$occ <- 0
proj4string(a100b) <- OSGB.proj
writeOGR(a100b, dsn=wd.distn, layer="a100m1990buf5km_dd", driver="ESRI Shapefile")

# ### 100m absences within grid-cells with > 5% broadleaved forest grid-cells, randomly selected to match 100m records from 1990 onwards.
# ## Note that in Wales investigation the occupied cells contained less broadleaved woodland than all cells! I judged this approach not very useful. 
# nfi <- raster(paste0(wd.env,"\\broadl_m.tif"))
# r <- matrix(c(-1,0.05,NA,  0.05,1,1), byrow=T, nrow=2)
# nfi <- reclassify(nfi, r) ## identify cells with 5% or more broadleaved coverage 
# nfi[p100.1990] <- NA ## remove cells that contain dormouse records from 1990 onwards
# 
# a100 <- as.data.frame(sampleRandom(nfi, n, xy=T, ext=e)) ## select some random grid-cells, more than needed
# # a100 <- a100[sample(c(1:nrow(a100)), n),] ## select correct number of grid-cells
# coordinates(a100) <- ~x+y ## Obtain the centroids of the 100m cells where absences are placed. coords
# a100@data$occ <- 0
# proj4string(a100) <- OSGB.proj
# writeOGR(a100, dsn=wd.distn, layer="a100m1990broadl5pc", driver="ESRI Shapefile")



