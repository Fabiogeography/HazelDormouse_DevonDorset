##############################################################
##### CORRELATE PREDICTIONS ########################
##### Written by: Regan Early ################################
##### Written on: 29th March 2022 ##############################
##### Modified on: 29th March 2022  #########################
##############################################################

.libPaths("D:/SOFTWARE/R-4.1.2/library")
library(raster)
library(caret) ## for confusion matrix and kappa

wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/100m/"

dd_ndmp_pa <- raster(paste0(wd.out,"/final_preds_ndmp_THaW.tiff")) ## Devon and Dorset model made with NDMP presences and pseudo-absences, trimmed to THaW region
dd_nddlrd_pa <- raster(paste0(wd.out,"/final_preds_nddlrd_THaW.tiff")) ## Devon and Dorset model made with NDD+LRD presences and pseudo-absences, trimmed to THaW region
t_nddlrd_ndmp_sa <- raster(paste0(wd.out,"/final_preds_t_sampled.tiff")) ## THaW model made with NDMP and NDDLRD presences and survey absences
t_nddlrd_sa <- raster(paste0(wd.out,"/final_preds_t_sampled_nddlrd.tiff")) ## THaW model made withNDDLRD presences and survey absences
t_nddlrd_pa <- raster(paste0(wd.out,"/final_preds_t_sampled_nddlrd_pabs.tiff")) ## THaW model made with NDDLRD presences and pseudo-absences
t_nddlrd_ndmp_nfilcm <- raster(paste0(wd.out,"/final_preds_t_sampled_NFI+LCM.tiff")) ## THaW model made with NDMP and NDDLRD presences and survey absences but without treehedge data

t_nddlrd_ndmp_sa <- crop(t_nddlrd_ndmp_sa, dd_ndmp_pa)
t_nddlrd_sa <- crop(t_nddlrd_sa, dd_ndmp_pa)
t_nddlrd_pa <- crop(t_nddlrd_pa, dd_ndmp_pa)
t_nddlrd_ndmp_nfilcm <- crop(t_nddlrd_ndmp_nfilcm, dd_ndmp_pa)

##### Investigate relationship between NDMP and NDD+LRD predictions #####
s <- stack(dd_ndmp_pa, dd_nddlrd_pa, t_nddlrd_ndmp_sa, t_nddlrd_sa, t_nddlrd_pa, t_nddlrd_ndmp_nfilcm)
# pairs(s) ## Guessing this uses Pearson's correlation - not stated in documentation - which I've found elseswhere to be much higher than Speaman's
s.dat <- na.omit(as.data.frame(values(s)))

## Correlation
cor(s.dat, method="spearman") ##
# cor(s.dat, method="pearson") ## 


##### Categorical #####
c.dd_ndmp_pa <- raster(paste0(wd.out,"/final_predss_ndmp_THaW.tiff")) ## Devon and Dorset model made with NDMP presences and pseudo-absences, trimmed to THaW region
c.dd_nddlrd_pa <- raster(paste0(wd.out,"/final_predss_nddlrd_THaW.tiff")) ## Devon and Dorset model made with NDD+LRD presences and pseudo-absences, trimmed to THaW region
c.t_nddlrd_ndmp_sa <- raster(paste0(wd.out,"/final_predss_t_sampled.tiff")) ## THaW model made with NDMP and NDDLRD presences and survey absences
c.t_nddlrd_sa <- raster(paste0(wd.out,"/final_predss_t_sampled_nddlrd.tiff")) ## THaW model made withNDDLRD presences and survey absences
c.t_nddlrd_pa <- raster(paste0(wd.out,"/final_predss_t_sampled_nddlrd_pabs.tiff")) ## THaW model made with NDDLRD presences and pseudo-absences
c.t_nddlrd_ndmp_nfilcm <- raster(paste0(wd.out,"/final_predss_t_sampled_NFI+LCM.tiff")) ## THaW model made with NDMP and NDDLRD presences and survey absences but without treehedge data

## Crop to same extent
c.t_nddlrd_ndmp_sa <- crop(c.t_nddlrd_ndmp_sa, c.dd_ndmp_pa)
c.t_nddlrd_sa <- crop(c.t_nddlrd_sa, dd_ndmp_pa)
c.t_nddlrd_pa <- crop(c.t_nddlrd_pa, dd_ndmp_pa)
c.t_nddlrd_ndmp_nfilcm <- crop(c.t_nddlrd_ndmp_nfilcm, dd_ndmp_pa)

### Extract values
xy <- xyFromCell(c.dd_ndmp_pa, 1:ncell(c.dd_ndmp_pa))
z.dd_ndmp_pa <- extract(c.dd_ndmp_pa, xy) ## as.factor(na.omit())
z.dd_nddlrd_pa <-extract(c.dd_nddlrd_pa, xy)
z.t_nddlrd_ndmp_sa <- extract(c.t_nddlrd_ndmp_sa, xy)
z.t_nddlrd_sa <- extract(c.t_nddlrd_sa, xy)
z.t_nddlrd_pa <- extract(c.t_nddlrd_pa, xy)
z.t_nddlrd_ndmp_nfilcm <- extract(c.t_nddlrd_ndmp_nfilcm, xy)

## Ensure NAs are removed in all vectors
na1 <- z.t_nddlrd_pa
na1[!is.na(na1)] <- 1
z.dd_ndmp_pa <- as.factor(na.omit(z.dd_ndmp_pa * na1))
z.dd_nddlrd_pa <- as.factor(na.omit(z.dd_nddlrd_pa * na1))
z.t_nddlrd_ndmp_sa <- as.factor(na.omit(z.t_nddlrd_ndmp_sa * na1))
z.t_nddlrd_sa <- as.factor(na.omit(z.t_nddlrd_sa * na1))
z.t_nddlrd_pa <- as.factor(na.omit(z.t_nddlrd_pa * na1))
z.t_nddlrd_ndmp_nfilcm <- as.factor(na.omit(z.t_nddlrd_ndmp_nfilcm * na1))

##

c <- confusionMatrix(z.dd_ndmp_pa, z.dd_nddlrd_pa) 
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.dd_ndmp_pa, z.t_nddlrd_ndmp_sa) 
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.dd_ndmp_pa, z.t_nddlrd_sa)
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.dd_ndmp_pa, z.t_nddlrd_pa) 
round(c$table / sum(c$table), 2)

##
c <- confusionMatrix(z.dd_nddlrd_pa, z.t_nddlrd_ndmp_sa) 
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.dd_nddlrd_pa, z.t_nddlrd_sa) 
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.dd_nddlrd_pa, z.t_nddlrd_pa) 
round(c$table / sum(c$table), 2)

##

c <- confusionMatrix(z.t_nddlrd_ndmp_sa, z.t_nddlrd_sa) 
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.t_nddlrd_ndmp_sa, z.t_nddlrd_pa) 
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.t_nddlrd_sa, z.t_nddlrd_pa)
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.dd_ndmp_pa, z.t_nddlrd_ndmp_nfilcm) 
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.dd_nddlrd_pa, z.t_nddlrd_ndmp_nfilcm) 
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.t_nddlrd_ndmp_sa, z.t_nddlrd_ndmp_nfilcm) 
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.t_nddlrd_sa, z.t_nddlrd_ndmp_nfilcm)
round(c$table / sum(c$table), 2)

c <- confusionMatrix(z.t_nddlrd_pa, z.t_nddlrd_ndmp_nfilcm) 
round(c$table / sum(c$table), 2)
