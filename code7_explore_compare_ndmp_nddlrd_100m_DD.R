##############################################################
##### MAKE AND EVALUATE FINAL MODEL ########################
##### Written by: Regan Early ################################
##### Written on: 9th December 2021 ##############################
##### Modified on: 9th December 2021  #########################
##############################################################

# .libPaths("C:\\Rpackages")
.libPaths("D:/SOFTWARE/R-4.1.1/library")
# library(dplyr)#, lib.loc="D:/SOFTWARE/R-4.1.1/library") ## bind_rows, sample_frac
library(raster)#, lib.loc="D:/SOFTWARE/R-4.1.1/library")
library(rgdal)#, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(pROC)#, lib.loc="D:/SOFTWARE/R-4.1.1/library")
library(modEvA)
library(tidyverse)
# library(lmodel2) ## Model II regression
# library(energy) ## mvnorm.etest

wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/100m"
wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
wd.distn <- "E:/NON_PROJECT/DORMOUSE_NE/DISTRIBUTION"

occ <- readOGR(dsn=wd.distn, layer="occ100m1990_dd") ## load, so always working with the same dataset
occ$source[is.na(occ$source)] <- "pseudo-absence"

p.ndmp <- raster(paste0(wd.out, "/final_preds_ndmp.tif"))
p.nddlrd <- raster(paste0(wd.out,"/final_preds_nddlrd_16Dec2021.tif"))

##### Combine continuous predictions from NDMP and NDD+LRD models - not used due to different prevalences #####
# p.max <- max(p.ndmp, p.nddlrd)               
# writeRaster(p.max, file=paste0(wd.out,"/final_pred_max.tiff"))
# 
# p.mean <- mean(p.ndmp, p.nddlrd)               
# writeRaster(p.mean, file=paste0(wd.out,"/final_pred_mean.tiff")) ## Doesn't look like the right approach

# ## How well does the combined model predict all presences?
# p.pred <- raster::extract(p.max, occ)
# r <- roc(occ$occ, p.pred)
# r # 0.8102 Not bad

# ## Is  this is better than the ndmp data alone?
# p.pred <- raster::extract(p.ndmp, occ)
# r <- roc(occ$occ, p.pred)
# r # 0.6445 - yes

##### Combine models by calculating continuous favourability and selecting the highest value in each grid-cell #####
### Calculate favourability
occ.ndmp <- occ[occ$source!="LRD+NDD",]
prev.ndmp <- table(occ.ndmp$source)[1] / table(occ.ndmp$source)[2] ## Ratio of presence ; absence in the ndmp model
f.ndmp <- exp(p.ndmp)/(prev.ndmp+exp(p.ndmp))
writeRaster(f.ndmp, file=paste0(wd.out,"/final_fav_ndmp.tiff"))

occ.nddlrd <- occ[occ$source!="NDMP",]
prev.nddlrd <- table(occ.nddlrd$source)[1] / table(occ.nddlrd$source)[2] ## Ratio of presence ; absence in the nddlrd model
f.nddlrd <- exp(p.nddlrd)/(prev.nddlrd+exp(p.nddlrd))

f.nddlrd <- crop(f.nddlrd, f.ndmp)
writeRaster(f.nddlrd, file=paste0(wd.out,"/final_fav_nddlrd.tiff"))

f.max <- max(f.ndmp, f.nddlrd)
writeRaster(f.max, file=paste0(wd.out,"/final_fav_max.tiff")) ## This is identical to the NDMP model

f.mean <- mean(f.ndmp, f.nddlrd)
writeRaster(f.mean, file=paste0(wd.out,"/final_fav_mean.tiff")) ## 

## How well does the combined model predict all presences?
f.pred <- raster::extract(f.max, occ)
r <- roc(occ$occ, f.pred)
r # 0.6445 :(

## Is this is better than the ndmp data alone?
f.pred <- raster::extract(f.ndmp, occ)
r <- roc(occ$occ, f.pred)
r # 0.6445 - no, identical, obviously

## How well does the mean prediction do?
f.pred <- raster::extract(f.mean, occ)
r <- roc(occ$occ, f.pred)
r # 0.6865 - A bit better than the max / NDMP prediction

### Apply thresholds
r <- roc(occ$occ, f.pred)
r # 0.6865
df <- data.frame(spec=r$specificities, sens=r$sensitivities, thresh=r$thresholds)

df[which.min(abs(0.95-df$sens)), ] ## The threshold that yields 95% sensitivity: 0.7539412
df[which.min(abs(df$spec-df$sens)), ] ## The threshold that minimises the difference between sensitivity and specificity: 0.7651836
df[which.max(df$spec+df$sens), ]## minimises difference between sensitivity and specificity 0.7643491

##### Combine models using ss threshold on each model individually, then combining #####

### What is the sens / spec / TSS at the chosen thresholds for the two models? Want them to be equal. 
## ss threshold, NDDLRD
occ.nddlrd <- occ
occ.nddlrd <- occ.nddlrd[occ.nddlrd$source!="NDMP",] ## remove points with non-NDD+LRD presences
p.nddlrd.df <- raster::extract(p.nddlrd, occ.nddlrd)
p.nddlrd.df.eval.ss <- threshMeasures(obs=occ.nddlrd$occ, pred=p.nddlrd.df, thresh=0.4673669) ## Sensitivity 0.61142857, specificity 0.61066236, sTSS 0.61104547. Values before 16th Dec were 0.86095238, Specificity 0.86106624, sTSS 0.86100931

## ss threshold, NDMP
occ.ndmp <- occ
occ.ndmp <- occ.ndmp[occ.ndmp$source!="LRD+NDD",] ## remove points with non-NDD+LRD presences
p.ndmp.df <- raster::extract(p.ndmp, occ.ndmp)
p.ndmp.df.eval.ss <- threshMeasures(obs=occ.ndmp$occ, pred=p.ndmp.df, thresh=0.2637337) ## Sensitivity 0.81250000, Specificity 0.81744750, sTSS 0.81497375

## 95% threshold
p.nddlrd.df.eval.95 <- threshMeasures(obs=occ.nddlrd$occ, pred=p.nddlrd.df, thresh=0.3752346) ## Sensitivity 0.95047619 (0.95047619), Specificity 0.1922455 (0.81421648), sTSS 0.619481499 (0.88234633)
p.ndmp.df.eval.95 <- threshMeasures(obs=occ.ndmp$occ, pred=p.ndmp.df, thresh=0.05963497) ## Sensitivity 0.94791667, Specificity 0.43941842, sTSS 0.69366754

### based on the above, the ss threshold results in fairly similar sens / spec / TSS for the two models. I think it's safe to go ahead and use these thresholded values to combine the two models
p.nddlrd.ndmp.ss <- max(p.nddlrdss, p.ndmpss)
p.nddlrd.ndmp.df <- raster::extract(p.nddlrd.ndmp.ss, occ)
p.nddlrd.ndmp.df.eval.ss <- threshMeasures(obs=occ$occ, pred=p.nddlrd.ndmp.df, thresh=0.5) ## Sensitivity 0.91626409, Specificity 0.75605816, sTSS 0.83616112

writeRaster(p.nddlrd.ndmp.ss, file=paste0(wd.out, "/final_predss_nddlrd_ndmp.tif"))

p.ndmp95 <- raster(paste0(wd.out, "/final_pred95_ndmp.tif"))
p.nddlrd95 <- raster(paste0(wd.out, "/final_pred95_nddlrd.tif"))
p.ndmpss <- raster(paste0(wd.out, "/final_predss_ndmp.tif"))
p.nddlrdss <- raster(paste0(wd.out, "/final_predss_nddlrd.tif"))

##### Combine models by retaining all points above the ss threshold in the NDMP model, and then add cells from NDDLRD model that improve predictions - FINAL METHOD #####
p.ndmpss <- raster(paste0(wd.out, "/final_predss_ndmp.tif"))

p.nddlrd <- crop(p.nddlrd, p.ndmpss)
test <- max(p.ndmpss, p.nddlrd) ## combine the rasters initially

## Calculate the ROC for the nddlrd model and the remaining distribution data
p.pred.combo <- raster::extract(test, occ)
r <- roc(occ$occ, p.pred.combo)
r # 0.6552 ## How well the nddlrd model + max ndmp predicts the presences

df <- data.frame(spec=r$specificities, sens=r$sensitivities, thresh=r$thresholds)

df[which.min(abs(0.95-df$sens)), ] ## The threshold that yields 95% sensitivity: 0.3801677
df[which.min(abs(df$spec-df$sens)), ] ## The threshold that minimises the difference between sensitivity and specificity: 0.4846393
df[which.max(df$spec+df$sens), ]## minimises difference between sensitivity and specificity 0.4626313

## Threshold the nddlrd model according to the ss method
rcl <- matrix(c(0,0.4846393,0,   0.4846393,1,1), nrow=2, byrow=T) 
pss.combo <- reclassify(test, rcl)
writeRaster(pss.combo, paste0(wd.out, "/final_predss_nddlrd_ndmp.tiff"))

##### Investigate relationship between NDMP and NDD+LRD predictions #####
s <- stack(p.ndmp, p.nddlrd)
s.dat <- na.omit(as.data.frame(values(s)))

## Correlation
cor(s.dat, method="spearman") ## 0.038
cor(s.dat, method="pearson") ## 0.367 - much higher than spearman value, indicating that the non-independence introduced by spatial autocorrelation affects the results 

## Contingency table
t <- stack(p.ndmpss, p.nddlrdss)
t.dat <- na.omit(as.data.frame(values(t)))
tt <- table(t.dat) 
tt[2,1]/(tt[2,1] + tt[2,2]) ## 54% of ndmp suitable sites are not suitable in nddlrd
tt[1,2]/(tt[1,2] + tt[2,2]) ## 49% of nddlrd suitable sites are not suitable in ndmp

## TSS of the thresholded final prediction on all records
p.combo.df.eval.ss <- threshMeasures(obs=occ$occ, pred=p.pred.combo, thresh=0.4846393) ## Sensitivity 0.613526570, specificity 0.613893376, sTSS 0.613709973.

## TSS of the thresholded NDMP prediction on all records
p.pred <- raster::extract(p.ndmpss, occ)
p.ndmp.df.eval.ss <- threshMeasures(obs=occ$occ, pred=p.pred, thresh=0.5) ## Sensitivity 0.53945250, specificity 0.66720517, sTSS 0.60332883

##### Investigate outliers in regression between NDMP and NDD+LRD predictions - not continued because predictions are so dissimilar it wouldn't be worthwhile #####
## LM
m1 <- lm(s.dat$final_preds_nddlrd ~ s.dat$final_preds_ndmp)
par(mfrow=c(2,2))
plot(m1) ## Doesn't look great

resid.m1 <- raster(p.ndmp) * NA
values(resid.m1)[-m1$na.action] <- m1$residuals

### Model II regression
## Check for bivariate normality
par(mfrow=c(2,2))
hist(s.dat[,1])
hist(log(s.dat[,1]))
hist(sqrt(s.dat[,1])) ## best but not great

par(mfrow=c(2,2))
hist(s.dat[,2])
hist(log(s.dat[,2])) ## best
hist(sqrt(s.dat[,2]))


plot(s.dat) ## umm no
plot(s.dat[,1], log(s.dat[,2]))
plot(sqrt(s.dat[,1]), log(s.dat[,2]))

# mvnorm.etest(s.dat, R=2)

lmodel2(s.dat$final_preds_nddlrd ~ s.dat$final_preds_ndmp, nperm=10)
lmodel2(sqrt(s.dat$final_preds_nddlrd) ~ log(s.dat$final_preds_ndmp), nperm=10) ## most normal version of data I can find


resid_lm <- raster(temp_chl_s, 1) * NA
values(resid_lm)[-lm1$na.action] <- lm1$residuals

##### COMPARE LANDSCAPE CHARACTERISTICS FOR NDMP AND NDDLRD PREDICTED HABITAT M#####
env.r <- stack(paste0(wd.env, "/env_100m7_dd"))

# ndmp95.env <- p.ndmp95 * env.r
# names(ndmp95.env) <- names(env.r)

# nddlrd95.env <- p.nddlrd95 * env.r
# names(nddlrd95.env) <- names(env.r)

ndmpss.env <- p.ndmpss * env.r
names(ndmpss.env) <- names(env.r)

nddlrdss.env <- p.nddlrdss * env.r
names(nddlrdss.env) <- names(env.r)

### Visual comparison

ndmp <- as.data.frame(na.omit(values(ndmpss.env)))
ndmp$source <- "ndmp"
nddlrd <- as.data.frame(na.omit(values(nddlrdss.env)))
nddlrd$source <- "nddlrd"

dat <- rbind(ndmp, nddlrd)
dat[dat==0] <- NA ## Otherwise 0 values swamp the data. 

vars <- c("coppice_m","copp_stds_m","yngtree_m","lodens_m","shrub_m","felled_m",
          "coppice_c","copp_stds_c","yngtree_c","lodens_c","shrub_c","felled_c",
          "combo_broadl_m","combo_conif_m","combo_broadl_c","combo_conif_c")

var.names <- c("NFI coppice", "NFI coppice stds", "NFI yngtree", "NFI lodens", "NFI shrub", "NFI felled",
               "NFI coppice change", "NFI coppice stds change", "NFI yngtree change", "NFi lodens change", "NFI shrub change", "NFI felled change",
               "LCM NFI broadl+mix", "LCM NFI conif+mix", "LCM NFI broadl+mix change", "LCM NFI conif+mix change"
)


dat.1 <- dat[,c("source", vars)]
colnames(dat.1) <- c("source", var.names)

jpeg(paste0(wd.out,"/hist_nddlrd_vs_ndmp_env.jpg"), width=1000, height=600)
dat.1 %>%
  gather(key = "var", value = "value", factor_key=T, -source) %>% ## factor_key means order of columns is preserved
  ggplot(aes(x=value, colour=source, fill=source)) +
  geom_histogram(alpha=0.6) +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
dev.off()



vars <- c("anc_rep_wood", "anc_sn_wood","OS_Terrain_100_NS","OS_Terrain_100_WE","OS_Terrain_100_slope_pct",
          "rainfall_spr_5yrmn","rainfall_smr_5yrmn","rainfall_aut_5yrmn",
          "sun_spr_5yrmn","sun_smr_5yrmn","sun_aut_5yrmn",
          "tasmax_spr_5yrmn","tasmax_smr_5yrmn","tasmax_aut_5yrmn",
          "tasmin_win_5yrmn","tasrng_win_5yrmn")

var.names <- c("Ancient replanted wood", "Ancient semi-nat wood", "Aspect NS", "Aspect WE", "Slope %",
                "Spring rainfall", "Summer rainfall", "Autumn rainfall",
                "Sun in spring", "Sun in summer", "Sun in autumn",
                "Max spring temperature", "Max summer temperature", "Max autumn temperature",
                "Min winter temperature", "Winter temperature range")

dat.2 <- dat[,c("source", vars)]
colnames(dat.2) <- c("source", var.names)

jpeg(paste0(wd.out,"/hist_nddlrd_vs_ndmp_env_2.jpg"), width=1000, height=600)
dat.2 %>%
  gather(key = "var", value = "value", factor_key=T, -source) %>% ## factor_key means order of columns is preserved
  ggplot(aes(x=value, colour=source, fill=source)) +
  geom_histogram(alpha=0.6) +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
dev.off()









par(mfrow=c(4,4))
for (v in vars) {
  ndmp <- values(ndmpss.env)[,v]
  ndmp <- na.omit(ndmp[ndmp>0])
  
  nddlrd <- values(nddlrdss.env)[,v]
  nddlrd <- na.omit(nddlrd[nddlrd>0])
  
  hist(ndmp, xlab=v, main="", col=rgb(0,0,1,0.2))
  hist(nddlrd, xlab=v, main="", col=rgb(1,0,0,0.2), add=T)
  
  legend('topright', c('NDMP', 'NDDLRD'),
         fill=c(rgb(0,0,1,0.2), rgb(1,0,0,0.2)))
  
}

par(mfrow=c(4,4))
for (v in vars) {
  h <- values(nddlrdss.env)[,v]
  h <- na.omit(h[h>0])
  hist(h, xlab=v, main="")
}

