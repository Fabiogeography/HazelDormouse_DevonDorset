##############################################################
##### UPDATE SDM PREDICTIONS BY CHANGING BROADLEAF COVERAGE ########################
##### Written by: Regan Early ################################
##### Written on: 16th December 2021 ##############################
##### Modified on: #########################
##############################################################

.libPaths("D:/SOFTWARE/R-4.1.1/library")
library(raster)
library(rgdal)

wd <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/100m/PHASE2_FOR_ZELDA"

### Environmental data - the standardised variables
ndmp.env <- stack(paste0(wd, "/env_100m7_dd_fullclim_ndmp_std"))
nddlrd.env <- stack(paste0(wd, "/env_100m7_dd_fullclim_nddlrd_std")) 

### Updated broadleaf raster (NFI + LCM broadleaved coverage)
## Note that a legacy issue means that the coverage of forest scales from 0-2 rather than 0-1. It doesn't make any difference to the SDM calculations, but has to be considered when manually altering the coverage of broadleaf habitat. E.g. if you want a coverage of 0.8 you should specify a value of 1.6 
combo_broadl_m <- raster(paste0(wd, "/combo_broadl_m.tif")) ## This is the original broadleaf raster
## Experiment with updating, e.g. reduce coverage in all cells by 15 percentage points
combo_broadl_m <- calc(combo_broadl_m, fun=function(x){ max(0, (x-0.3)) })

### Convert raw environmental data to the standardised versions
destdize.ndmp <- read.csv(paste0(wd, "/dat_means_sds_30Nov2021_dd_ndmp.csv"))
destdize.nddlrd <- read.csv(paste0(wd, "/dat_means_sds_30Nov2021_dd_nddlrd.csv")) 

ndmp.env[["combo_broadl_m"]] <- (combo_broadl_m - destdize.ndmp[1,"combo_broadl_m"]) / destdize.ndmp[2,"combo_broadl_m"] ## Subtract mean of original data and divide by standard deviation
nddlrd.env[["combo_broadl_m"]] <- (combo_broadl_m - destdize.nddlrd[1,"combo_broadl_m"]) / destdize.nddlrd[2,"combo_broadl_m"] ## Subtract mean of original data and divide by standard deviation

### Original SDMs
load(file=paste0(wd,"/mods_final_ndmp")) ## The GLM for NDMP data
ndmp.m <- final; rm(final)

load(file=paste0(wd,"/mods_final_nddlrd")) ## The GLM for NDD+LRD data
nddlrd.m <- final; rm(final)

### Make the prediction for the averaged model
ndmp.p <- predict(ndmp.env, ndmp.m, re.form=NA, type="response", full=T, se.fit=F)
writeRaster(ndmp.p, file=paste0(wd,"/update_preds_ndmp.tiff"), overwrite=T)

nddlrd.p <- predict(nddlrd.env, nddlrd.m, re.form=NA, type="response", full=T, se.fit=F)
writeRaster(nddlrd.p, file=paste0(wd,"/update_preds_nddlrd.tiff"), overwrite=T)

##### Threshold and combine rasters #####
## Threshold the ndmp model according to the ss method
rcl <- matrix(c(0,0.2637337,0,   0.2637337,1,1), nrow=2, byrow=T) 
ndmp.pss <- reclassify(ndmp.p, rcl)

## Combine the thresholded ndmp raster and the raw nddlrd raster
combo <- max(ndmp.pss, nddlrd.p)

## Threshold the combined model according to the ss method, calculated when combining SDM predictions (code 7)
rcl <- matrix(c(0,0.4846393,0,   0.4846393,1,1), nrow=2, byrow=T) 
combo.pss <- reclassify(combo, rcl)
writeRaster(combo.pss, paste0(wd, "/update_predss_nddlrd_ndmp.tiff"))


