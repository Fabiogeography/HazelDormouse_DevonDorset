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
# cell <- cells[100]

clumpy.fn <- function (cell) {
  # Calculate the extent of the target 100m cell
  xy <- as.data.frame(xyFromCell(m100.t, cell)) ## coordinates of the focal 100mcell
  xmin <- xy$x-50
  xmax <- xy$x+50
  ymin <- xy$y-50
  ymax <- xy$y+50
  e <- extent(xmin, xmax, ymin, ymax)
  
  # Extract the 1m cells that are within the 100m cell
  f1 <- crop(treehedge, e)
  # f1 <- raster::mask(f1, buf500)

  if(cellStats(f1, sum) == 0) { ## There are no trees/hedge in the 500m buffer
    clumpy <- connect <- pacfrac <- -3
    ta <- 0
  } 
  if(cellStats(f1, sum) == 1) { ## Only one pixel in the 500m buffer is trees/hedge
    clumpy <- connect <- pacfrac <- -2
    ta <- 1
  } 
  if(cellStats(f1, sum) == 10000) { ## All values in the 500m buffer are trees/hedge
    clumpy <- connect <- pacfrac <- 2
    ta <- raster::area(buf500) ## m2
  } 
  if(!(cellStats(f1, sum) %in% c(0, 1, 10000))) { ## Clumpiness can be calculated
    clumpy <- lsm_c_clumpy(f1)
    clumpy <- clumpy[clumpy$class==1,]$value ## The clumpiness of the forest only
    connect <- lsm_l_contag(f1)$value ## name is 'connectance' in list_lsm
    pacfrac <- lsm_l_pafrac(f1)$value
    # prox_cv <- ... enn?
    m <- matrix(c(-1,0.5,NA), byrow=T, nrow=1)
    test <- reclassify(f1, m)
    ta <- lsm_l_ta(test)$value ## total area of focal habitat type

  }
  out <- cbind(cell, coordinates(xy), clumpy, connect, pacfrac, ta)
  print(out)
  write.table(out, file=paste0(wd.out, "/fragstats_100m_treehedge_t.csv"), sep=",", append=TRUE, quote=FALSE,
                        col.names=FALSE, row.names=FALSE)
  
  return(out)
}

### Apply over all cells in 100m raster
fragstats <- pblapply(cells, clumpy.fn) ## 
which(cells==156482)/length(cells) ## Calculate how far through the list of cells we are
# cell <- cells[which(cells==152111)+1] ## The cell id of the value that failed
# plot(clumpy$prop, clumpy$value, ylim=c(0.99, 1))

# write.csv(fragstats, paste0(wd.out, "/fragstats_500m_treehedge_t.csv"), row.names = FALSE)

f <- read.csv(paste0(wd.out, "/fragstats_100m_treehedge_t.csv"), header=F)
colnames(f) <- c("cell", "x", "y", "clumpy", "connect", "pacfrac", "ta")

### Assign values to the raster
clumpy100 <- connect100 <- pacfrac100 <- ta100 <- m100.t

clumpy100[cells] <- f$clumpy
connect100[cells] <- f$connect
pacfrac100[cells] <- f$pacfrac
ta100[cells] <- f$ta

### Plot results
jpeg(paste0(wd.out,"/fragstats100_treehedge.jpg"), width=600, height=800)
par(mfrow=c(2,2))
IS THIS RIGHT?????
rcl <- matrix(c(-4057,-0.1,0,  1.01,Inf,1), nrow=2, byrow=T) ## Don't know why there are values as low as -4056 or above 1, but checking in Arcmap confirms that these are cells with no or complete coverage of treehedge, so reclassing as 0 and 1 respectively is apprpriate
clumpy100.plot <- reclassify(clumpy100, rcl) ## Approx 50 cells with values above 1 and less than 4.5. Didn't expect valid values >1.
plot(clumpy100.plot, main="Clumpiness")

plot(connect100, main="Connectance")

rcl <- matrix(c(-3.5,0.99999,NA,  2.0,2.0,NA), nrow=2, byrow=T) ## Values should between 1 and 2.
pacfrac100.plot <- reclassify(pacfrac100, rcl)
plot(pacfrac100.plot, main="Perimeter area fractal dimension")

plot(ta100, main="Total area")

dev.off()

### Write rasters
# writeRaster(clumpy100, paste0(wd.out,"/clumpy100_treehedge_errors.tif"))
writeRaster(clumpy100.plot, paste0(wd.out,"/clumpy100_treehedge.tif"))
writeRaster(connect100, paste0(wd.out,"/connect100_treehedge.tif"))
writeRaster(pacfrac100.plot, paste0(wd.out,"/pacfrac100_treehedge.tif"))
writeRaster(ta100, paste0(wd.out,"/ta100_treehedge.tif"))

### Correlate rasters
s <- stack(clumpy100.plot, connect100, pacfrac100.plot, ta100)
names(s) <- c("clumpy100", "connect100", "pacfrac100", "ta100")
layerStats(s, 'pearson', na.rm=T)$'pearson correlation coefficient' ## Greatest correlation is between clumpy and pacfrac, at -0.566

### Original code has some nice code for making bivariate plots of area and clumpiness ###

