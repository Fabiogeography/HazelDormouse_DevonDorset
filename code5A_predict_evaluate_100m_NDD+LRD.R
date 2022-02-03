##############################################################
##### PREDICT AND EVALUATE ########################
##### Written by: Regan Early ################################
##### Written on: 27th March 2021 ##############################
##### Modified on: 7th December 2021  #########################
##############################################################

# .libPaths("C:\\Rpackages")
.libPaths("D:/SOFTWARE/R-4.1.1/library")
library(dplyr) ## bind_rows, sample_frac
library(raster)
library(rgdal)
library(pROC)
library(MuMIn) # dredge() stdize() model.sel()

wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/100m/"
wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
wd.clim <- "E:/NON_PROJECT/DORMOUSE_WALES/GIS/100m" # "E:/GIS_DATA/CLIMATE/UK/HadUK-Grid_1km"
wd.distn <- "E:/NON_PROJECT/DORMOUSE_NE/DISTRIBUTION"
# wd.out <- "~/re259/UoE_U_Drive/DORMOUSE_NE/SDM/100m"
# wd.env <- "~/re259/UoE_U_Drive/DORMOUSE_NE/SDM/100m"

load(file=paste0(wd.out,"/mods_final_nddlrd"))

destdize <- read.csv(paste0(wd.out, "/dat_means_sds_30Nov2021_dd_nddlrd.csv"))

# w <- readOGR(dsn="E:/GIS_DATA/UK/bdline_essh_gb/Data/Wales", layer="Wales_mainland")

###########################################################################################
##### Test robustness of model with the same variables calibrated with different data #####
###########################################################################################
### Predict the original validation data
dat.v <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_v_nddlrd.csv"))
dat.m <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_m_nddlrd.csv"))

vars <- as.vector(colnames(final$coefficients))
vars <- vars[!grepl("Intercept|I|:", vars)] ## remove intercept, quadratic, and interaction terms from the variables

dat <- read.csv(paste0(wd.env,"/dat_30Nov2021_dd.csv"))
dat <- as.data.frame(cbind(dat[,1:7], dat[,vars])) ## Checked and this is the same data as used on the cluster - this is not the problem.

p <- predict(final, dat.v, re.form=NA, type="response") ## Predict model with validation data

roc_obj <- roc(dat.v$occ, p)
auc(roc_obj) ## 0.6163. Aargh!

### Test AUC on multiple validation datasets. Make model with same variables with new calibration and validation data
vars.mlqdi <- colnames(final$coefficients)
vars.mlqdi <- vars.mlqdi[vars.mlqdi!="(Intercept)"]

auc.fn <- function(x) {
  dat.cx <- sample_frac(dat.m, size=0.7)
  dat.vx <- dat[!(dat.m$id %in% dat.cx$id),]
  m <- glm(as.formula(paste0("occ ~ ", vars.mlqdi)), family=binomial, dat.cx) ## Make model with calibration data
  px <- predict(m, dat.vx, re.form=NA, type="response") ## Predict model with validation data
  
  roc_obj <- roc(dat.vx$occ, px)
  print(round(auc(roc_obj),3))
}

eval <- unlist(lapply(1:3, auc.fn))
mean(eval) ## 0.6343
sd(eval) ## 0.0093

###########################################################
##### Predict suitability for all Devon and Dorset ########
###########################################################

##### ENVIRONMENT DATA #####
env.r <- stack(paste0(wd.env, "/env_100m7_dd_16thDec2021"))
env.r <- env.r[[vars]]

### Need to add climate for all grid-cells, not just ones with distribution data
clim.vars <- c("tasmin_win")

i <- 2019 ## The 5 year climate mean is calculated for the most recent period
s <- stack(paste0(wd.clim, "/",clim.vars,"_5yrmn_",i,".tif"))

s <- crop(s, extent(env.r))

nn <- names(s)[grepl("_5yrmn", names(s))]
n <- nn[1]

for (n in nn) {
  n2 <- substr(n, start=1, stop=nchar(n)-5) ## trim year from end
  env.r[[n2]] <- s[[n]]
} ## These values are higher than the ones in teh raster I originally wrote on 7th Dec, even though env_100m7_dd was written on 1st Dec

writeRaster(env.r, file=paste0(wd.env, "/env_100m7_dd_fullclim_nddlrd_16thDec2021"))  ## This produces different values for *all* variables, even tasmin, than the file I have saved (though it might just be a sampling issue for some of themn)

### Put environmental data through same standardisation process as the original environmental data
# vars.mlqdi.main <- vars.mlqdi[!grepl("I", vars.mlqdi)]
# vars.mlqdi.main <- vars.mlqdi.main[!grepl(":", vars.mlqdi.main)]
## Enter manually in the interest of time - need to remove '_5yrmn' from names
# vars.mlqdi.main <- c("broadl_all_m","lcm_broadl_m","lcm_conif_m","OS_Terrain_100_NS", "OS_Terrain_100_slope_pct",
                     # "rainfall_aut_5yrmn","sun_spr_5yrmn", "tasmin_win_5yrmn","lcm_broadl_c","broadl_c" )

# n <- vars[1]
for (n in vars) {
  std <- destdize[,n] ## values of model input environment data to standardise raster values
  env.r[[n]] <- (env.r[[n]] - std[1]) / std[2] ## subtract original mean and divide by original standard deviation
}

raster::writeRaster(env.r, paste0(wd.env,"/env_100m7_dd_fullclim_nddlrd_std_16thDec2021")) ## Some variables have been added
# env.r.orig <- stack(paste0(wd.env, "/env_100m7_dd_fullclim_nddlrd_stdX"))  ## This produces different values for combo_broadl_m and tasmin than the fine I have saved

### Make the prediction for the averaged model
p <- predict(env.r, final, re.form=NA, type="response", full=T, se.fit=T)
writeRaster(p, file=paste0(wd.out,"/final_preds_nddlrd_16Dec2021.tiff"))
# p <- raster(paste0(wd.out,"/final_preds_nddlrd.tif"))

### Standard error of the prediction
env.dat <- as.data.frame(extract(env.r, 1:ncell(env.r)))
p.nddlrd.se <- stats::predict(final, newdata=env.dat, re.form=NA, type="response", full=T, se.fit=T)
test <- env.r[[1]]
values(test) <- NA
values(test) <- p.nddlrd.se$se.fit
writeRaster(test, file=paste0(wd.out,"/final_preds_nddlrd_se_16Dec2021.tiff"))

##### Calculate thresholds #####

## Inspect the sensitivity and specificity
occ <- readOGR(dsn=wd.distn, layer="occ100m1990_dd") ## load, so always working with the same dataset

occ.nddlrd <- occ
occ.nddlrd$source[occ.nddlrd$occ==0] <- "pseudo-absence"
occ.nddlrd <- occ.nddlrd[occ.nddlrd$source!="NDMP",] ## remove points with non-NDD+LRD presences

p.pred <- raster::extract(p, occ.nddlrd)
r <- roc(occ.nddlrd$occ, p.pred)
r # 0.6584 (as opposed to 0.8887, the new value seems more likely to be correct)

df <- data.frame(spec=r$specificities, sens=r$sensitivities, thresh=r$thresholds)

df[which.min(abs(0.95-df$sens)), ] ## The threshold that yields 95% sensitivity: 0.3752346 (was 0.06211754)
df[which.min(abs(df$spec-df$sens)), ] ## The threshold that minimises the difference between sensitivity and specificity: 0.4673669 (was 0.1227153)
df[which.max(df$spec+df$sens), ]## minimises difference between sensitivity and specificity 0.443513 (was 0.06589685)

### Reclassify prediction raster according to two thresholds (from pa_model3_eval_ply)
## 95% sensitivity threshold
rcl <- matrix(c(0,0.3752346,0,   0.3752346,1,1), nrow=2, byrow=T) 
p95 <- reclassify(p, rcl)
writeRaster(p95, file=paste0(wd.out,"/final_pred95_nddlrd_16Dec2021.tiff"))

## Threshold that minimises difference between sensitivity and specificity
rcl <- matrix(c(0,0.4673669,0,   0.4673669,1,1), nrow=2, byrow=T) 
pss <- reclassify(p, rcl)
writeRaster(pss, file=paste0(wd.out,"/final_predss_nddlrd_16Dec2021.tiff"))

# ## Manual threshold - to remove the large area of over-prediction in the SW, which appears to be caused by climate.
# ## Increasing the threshold from the ss threshold causes very little loss of 'suitable' habitat elsewhere, particularly where dormouse is present
# rcl <- matrix(c(0,0.68,0,   0.68,1,1), nrow=2, byrow=T) 
# pm <- reclassify(p, rcl)
# writeRaster(pm, file=paste0(wd.out,"/final_predm.tiff"))
# pm <- raster(paste0(wd.out,"/final_predm.tif"))

########################################
##### Calculate R-squared for all models in best subset #####
########################################
dat.c <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_c_nddlrd.csv"))

load(paste0(wd.out,"/linear_quadratic_interactions_dredged_nddlrd_combo1"))
load(paste0(wd.out,"/linear_quadratic_interactions_dredged_nddlrd_combo2"))

mlqid.1 <- get.models(mlqid.1, subset=delta<2) ## mlqid.1 is an object containing the best models. 
mlqid.2 <- get.models(mlqid.2, subset=delta<2) ## mlqid.2 is an object containing the best models. 

rsq.fn <- function(x) {
  (x$null.deviance - x$deviance) / x$null.deviance
}

rsq1 <- unlist(lapply(mlqid.1, rsq.fn))
rsq2 <- unlist(lapply(mlqid.2, rsq.fn))

mean(c(rsq1, rsq2)) ## 0.05086
sd(c(rsq1, rsq2)) ## 0.00305

