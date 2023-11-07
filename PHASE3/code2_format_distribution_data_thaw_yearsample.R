##### MAKE PRESENCE AND ABSENCE SHAPEFILES ########################
##### In this variant of the code the presences are sampled to match the years of the survey absences #####
##### Written by: Regan Early ################################
##### Written on: 7th February 2022 ##############################
##### Modified on: 21st February 2022  #########################
##############################################################

.libPaths("C:\\SOFTWARE\\R-4.1.2\\library")
library(rgdal)
library(raster)
library(sp) ## over
library(FNN) ## get.knnx
library(dplyr) ## distinct
library(corrplot)

wd.distn <- "E:\\NON_PROJECT\\DORMOUSE_NE\\DISTRIBUTION"
wd.env <- "E:\\NON_PROJECT\\DORMOUSE_NE\\GIS\\100m"

OSGB.proj <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
wgs.proj <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

grid100m <- raster("E:\\GIS_DATA\\UK\\OSGB_GRIDS\\OSGB_Grid_100m_clim.tif") ## This wsa made in ArcMap by converting the original raster to a points shapefile, trimming to the 1km_grid_region_clim polygon shapefile, and converting back to raster. Process in R didn't produce unique cell ids for each 100m grid-cell

### Thaw region
t100 <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/thaw_region.tif")

##### ABSENCES #####
# a <- read.csv(paste0(wd.distn,"\\SURVEY_ABSENCES\\surv_ab_24Jan2022.csv"))
# a$x <- a$easting ## make new coordinate columns so that easting and northing can remain in @data
# a$y <- a$northing
# coordinates(a) <- ~x+y
# proj4string(a) <- OSGB.proj
# writeOGR(a, dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="XYsurv_ab_24Jan2022", driver="ESRI Shapefile")
# a <- readOGR(dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="XYsurv_ab_24Jan2022")

### Trim down columns in presence data. See metadata file for column definitions. SiteName and Source may be useful down the line.
# a <- a[,c("site","easting", "northing", "year", "ptes", "sei")]
# colnames(a@data) <- c("site", "easting","northing", "year", "source", "sei")

### Trim to THaW region
# a100 <- a
# a100$m100 <- raster::extract(t100, a100)
# summary(a100$m100) ## 246 fall outside THaW region
# a100 <- a100[!is.na(a100$m100),] ## remove points that fall outside env data. 335 remaining records

### Filter to point point per 100m grid-cell. Keep the row with the most recent date (some cells have > 1 record on the same date, and one of these is chosen randomly)
## Don't trim to 1990 since all points are from 2002-2021.
a100$m100 <- raster::extract(grid100m, a100) ## assign unique value to points in each 1km grid-cell
a100 <- a100[!is.na(a100$m100),] ## remove points that fall outside env data. 1 record removed
a100 <- a100[order(a100$year, decreasing=TRUE),] ## 335 points
a100 <- a100[!duplicated(a100$m100),] ## 309 in seperate 1km grid-cells records. Nice.
# 
# a100$occ <- 0
# writeOGR(a100, dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="a100m_t", driver="ESRI Shapefile")
# a100 <- readOGR(dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="a100m_t")

##### Make occupancy file (in PHASE2 this was created in code 3 - explore) #####
a <- readOGR(dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="a100m_t")

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

##### Grid and filter at 100m #####
## Obtain 100m grid and trim to match the environmental data more closely

## crop and mask to Devon Dorset area
dd <- readOGR(dsn="E:/NON_PROJECT/DORMOUSE_NE/GIS", layer="devondorset_buffer5km", stringsAsFactors = F) ## Devon Dorset area
grid100m <- crop(grid100m, extent(dd))
grid100m <- raster::mask(grid100m, dd) ## mask is also in Hmisc

p100 <- ptes
p100$m100 <- raster::extract(grid100m, p100)
summary(p100$m100) ## 24 points fall outside the environmental data. Too few to worry about. 
p100 <- p100[!is.na(p100$m100),] ## remove points that fall outside env data

## How many grid-cells have both NDMP and NRD+LDD records?
test <- na.omit(p100@data[,c("source", "m100")]) ## 3072 rows
test <- distinct(test) ## 624 rows
length(table(test$m100)) ## 620, so 4 grid-cells have records from both NDMP and LRD+NDD sources. Too few to worry about. 
tail(sort(table(test$m100))) ## Identify the grid-cells 

## Trim to desired time period (matches presence data)
p <- p100[p100$year %in% c(2011:2019),] ## Legacy code. p100$year is already within this range
a <- a[a$year %in% c(2011:2019),]  ## Note there are no absences in 2011. There are 26 absences in 2020

### Filter to point point per 100m grid-cell. Sample years according to distribution within the absence data. 
set.seed(1) # make the order reproducible
p.sampled <- lapply(split(p, p$m100),
              function(subdf){
                if(nrow(subdf) >1) {   
                  x <- sample(a$year, size=nrow(subdf), replace=T) ## create a distribution of absence years the length needed to sample the presence years
                  y <- sample(subdf$year, size=1, prob=x)
                } else { y <- subdf$year}
                
                out <- cbind(subdf[1, c("IDorig", "easting", "northing", "source", "m100")], y)
                names(out) <- c("IDorig", "easting", "northing", "source", "m100", "year")
                return(out)
              }
)

p.sampled <- do.call('rbind', p.sampled) ## convert to data frame

### Compare to using a random year. set.seed(1) # make the order reproducible
rows <- sample(nrow(p100))
p100 <- p100[rows,] ## 3096 records
random <- p100[!duplicated(p100$m100),]

cor(random$year, p.sampled$year, method="spearman") ## low correlation! -0.06249314

par(mfrow=c(2,2))
hist(p$year, breaks=c(2011:2019), main="All presence records")
hist(p.sampled$year, breaks=c(2011:2019), main="Presence records sampled according to temporal distribution of absences") # More unifirm distribution
hist(a$year, breaks=c(2011:2019), main="All absences")
hist(random$year, breaks=c(2011:2019), main="Presence records sampled randomly")

##### For report #####
png(paste0(wd.out, "/THaW_occ_years.png"), width=800, height=300)

par(mfrow=c(1,3))
hist(a$year, breaks=c(2011:2019), main="All absences", xlab="Year")
hist(p$year, breaks=c(2011:2019), main="All presences", xlab="Year")
hist(occ$year[occ$source %in% c("LRD+NDD", "NDMP")], breaks=c(2011:2019),
     main="Presence records sampled according\nto temporal distribution of absences",
     xlab="year")
dev.off()

##### Continue #####
par(mfrow=c(1,1))
test <- merge(p.sampled@data, random@data, by.x="m100", by.y="m100", all=T)
sunflowerplot(test$year.x, test$year.y, xlab="Year sampled", ylab="Year random")

## Output sampled data
p.sampled$occ <- 1
writeOGR(p.sampled, dsn=wd.distn, layer="p100m1990_dd_sampled", driver="ESRI Shapefile")

## Trim presences to THaW area
p.sampled$t <- raster::extract(t100, p.sampled)
p100t <- p.sampled[!is.na(p.sampled$t),]
p100t$t <- NULL

## Save
writeOGR(p100t, dsn=wd.distn, layer="p100m1990_t_sampled", driver="ESRI Shapefile")

##### Make occupancy file #####
p <- readOGR(dsn=wd.distn, layer="p100m1990_t_sampled")
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

writeOGR(occ, dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="occ100mt_sampled", driver="ESRI Shapefile")
