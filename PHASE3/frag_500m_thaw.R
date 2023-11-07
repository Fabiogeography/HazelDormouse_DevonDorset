##############################################################
##### CALCULATE FRAGMENTATION WITHIN A 500M RADIUS AROUND TEH FOCAL 100M CELL#############################
##### Written by: Shari Mang and Hannah Wauchope ################################
##### Modified by: Regan Early ###################
##### Modified on: 15 Feb 2022 #############################
##############################################################

.libPaths("C:\\SOFTWARE\\R-4.1.2\\library")
library(raster)
library(rgdal)
library(landscapemetrics)
library(landscapetools)
library(pbapply) ## adds progress bar to apply
library(data.table) ## rbindlist - makes one data table of many
library(tidyverse) ## tibble functions
library(rgeos)
library(sp) ## coordinates

wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
OSGB.proj <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'

### Get the 100m grid on which to calculate fragmentation. Created by original frag_500m code.
m100.dd <- raster(paste0(wd.env, "/OSGB_Grid_100m_clim_dd.tif"))
t100 <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/thaw_region.tif")
m100.t <- m100.dd*t100

## Vegetation from THaW data
treehedge <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/trees_hedge_1m.tif")
scrub <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/scrub_1m.tif")

## Get 100m cell numbers inside the ThaW region
xy <- xyFromCell(m100.t, 1:ncell(m100.t))
xyz <- na.omit(cbind(xy, raster::extract(m100.t, xy)))
cells <- cellFromXY(m100.t, xyz[,c("x","y")])
rm(xy)
# name used later

### clumpiness calculations

clumpy.fn <- function (cell) {
  # Make a mask for the target 100m cell and the surrounding cells within a 500m radius
  xy <- as.data.frame(xyFromCell(m100.t, cell)) ## coordinates of the focal 100mcell
  coordinates(xy) <- ~x+y ## make spatial object
  buf500 <- gBuffer(xy, width=500) ## calculate the buffer
  
  # Extract the 1m cells that are within the buffer around the focal 100m cell
  f1 <- crop(treehedge, extent(buf500))
  f1 <- raster::mask(f1, buf500)

  if(cellStats(f1, sum) == 0) { ## There is no trees/hedge in the 500m buffer
    clumpy <- connect <- pacfrac <- -3
    ta <- 0
  } 
  if(cellStats(f1, sum) == 1) { ## Only one pixel in the 500m buffer is trees/hedge
    clumpy <- connect <- pacfrac <- -2
    ta <- 1
  } 
  if(cellStats(f1, sum) == 10000) { ## All values in the 500m buffer are trees/hedge
    clumpy <- connect <- pacfrac <- 3 ## (not 2, because pacfrac can be up to 2)
    ta <- raster::area(buf500) ## m2
  } 
  if(!(cellStats(f1, sum) %in% c(0, 1, 10000))) { ## Clumpiness can be calculated
    clumpy <- lsm_c_clumpy(f1)
    clumpy <- clumpy[clumpy$class==1,]$value ## The clumpiness of the forest only
    connect <- lsm_l_contag(f1)$value ## name is 'connectance' in list_lsm
    pacfrac <- lsm_l_pafrac(f1)$value ## perimter area fractal dimension. values will be between 1 and 2
    # prox_cv <- ... enn?
    m <- matrix(c(-1,0.5,NA), byrow=T, nrow=1)
    test <- reclassify(f1, m)
    ta <- lsm_l_ta(test)$value ## total area of focal habitat type

  }
  out <- cbind(cell, coordinates(xy), clumpy, connect, pacfrac, ta)
  # f1.clumpy$prop <- cellStats(f1, sum)/10000 ## this is the proportion of teh entire buffer, and the buffer may include non-landscape areas
  print(out)
  write.table(out, file=paste0(wd.out, "/fragstats_500m_treehedge_t.csv"), sep=",", append=TRUE, quote=FALSE,
                        col.names=FALSE, row.names=FALSE)
  

  return(out)
}

### Apply over all cells in 100m raster
fragstats <- pblapply(cells[1:57620], clumpy.fn) ## 1:ncell(m100)
which(cells==156482)/length(cells) ## Calculate how far through the list of cells we are
# cell <- cells[which(cells==152111)+1] ## The cell id of the value that failed
# plot(clumpy$prop, clumpy$value, ylim=c(0.99, 1))

# write.csv(fragstats, paste0(wd.out, "/fragstats_500m_treehedge_t.csv"), row.names = FALSE)

f <- read.csv(paste0(wd.out, "/fragstats_500m_treehedge_t.csv"), header=F)
colnames(f) <- c("cell", "x", "y", "clumpy", "connect", "pacfrac", "ta")

### Assign values to the raster
clumpy500 <- connect500 <- pacfrac500 <- ta500 <- m100.t

clumpy500[cells] <- f$clumpy
connect500[cells] <- f$connect
pacfrac500[cells] <- f$pacfrac
ta500[cells] <- f$ta

### Plot results
jpeg(paste0(wd.out,"/fragstats500_treehedge.jpg"), width=600, height=800)
par(mfrow=c(2,2))
rcl <- matrix(c(-515045,-0.1,NA,  1.01,Inf,NA), nrow=2, byrow=T) ## Don't know why there are values as low as -515044
clumpy500.plot <- reclassify(clumpy500, rcl) ## Approx 50 cells with values above 1 and less than 4.5. Didn't expect valid values >1.
plot(clumpy500.plot, main="Clumpiness")

plot(connect500, main="Connectance")

rcl <- matrix(c(-3.5,0.99999,NA,  2.0,2.0,NA), nrow=2, byrow=T) ## Values shouls beetween 1 and 2.
pacfrac500.plot <- reclassify(pacfrac500, rcl)
plot(pacfrac500.plot, main="Perimeter area fractal dimension")

plot(ta500, main="Total area")

dev.off()

### Write rasters
writeRaster(clumpy500.plot, paste0(wd.out,"/clumpy500_treehedge.tif"))
writeRaster(connect500, paste0(wd.out,"/connect500_treehedge.tif"))
writeRaster(pacfrac500.plot, paste0(wd.out,"/pacfrac500_treehedge.tif"))
writeRaster(ta500, paste0(wd.out,"/ta500_treehedge.tif"))

### Correlate rasters
s <- stack(clumpy500.plot, connect500, pacfrac500.plot, ta500)
names(s) <- c("clumpy500", "connect500", "pacfrac500", "ta500")
layerStats(s, 'pearson', na.rm=T)$'pearson correlation coefficient' ## Greatest correlation is between clumpy and connect, at -0.515

### Original code has some nice code for making bivariate plots of area and clumpiness ###
