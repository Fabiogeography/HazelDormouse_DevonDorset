##############################################################
##### PREDICT AND EVALUATE ########################
##### Written by: Regan Early ################################
##### Written on: 27th March 2021 ##############################
##### Modified on: 8th December 2021  #########################
##############################################################

# .libPaths("C:\\Rpackages")
.libPaths("D:/SOFTWARE/R-4.1.1/library")
library(dplyr)#, lib.loc="D:/SOFTWARE/R-4.1.1/library") ## bind_rows, sample_frac
library(raster)#, lib.loc="D:/SOFTWARE/R-4.1.1/library")
library(rgdal)#, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(pROC)#, lib.loc="D:/SOFTWARE/R-4.1.1/library")
library(MuMIn)#, lib.loc="D:/SOFTWARE/R-4.1.1/library") # dredge() stdize() model.sel()

wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/100m"
wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
wd.clim <- "E:/NON_PROJECT/DORMOUSE_WALES/GIS/100m" # "E:/GIS_DATA/CLIMATE/UK/HadUK-Grid_1km"
wd.distn <- "E:/NON_PROJECT/DORMOUSE_NE/DISTRIBUTION"
# wd.out <- "~/re259/UoE_U_Drive/DORMOUSE_NE/SDM/100m"
# wd.env <- "~/re259/UoE_U_Drive/DORMOUSE_NE/SDM/100m"

load(file=paste0(wd.out,"/mods_final_ndmp"))

destdize <- read.csv(paste0(wd.out, "/dat_means_sds_30Nov2021_dd_ndmp.csv"))

# w <- readOGR(dsn="E:/GIS_DATA/UK/bdline_essh_gb/Data/Wales", layer="Wales_mainland")

###########################################################################################
##### Test robustness of model with the same variables calibrated with different data #####
###########################################################################################
### Predict the original validation data
dat.v <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_v_ndmp.csv"))
dat.m <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_m_ndmp.csv"))

vars <- as.vector(colnames(final$coefficients))
vars <- vars[!grepl("Intercept|I|:", vars)] ## remove intercept, quadratic, and interaction terms from the variables

dat <- read.csv(paste0(wd.env,"/dat_30Nov2021_dd.csv"))
dat <- as.data.frame(cbind(dat[,1:7], dat[,vars]))

p <- predict(final, dat.v, re.form=NA, type="response") ## Predict model with validation data

roc_obj <- roc(dat.v$occ, p) ## AUC for NDMP validation records
auc(roc_obj) ## 0.8719. Nice!

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
mean(eval) ## 0.511 ## awful, but not sure that this is a good way to test the model
sd(eval) ## 0.01513

### Predict the NDD+LRD validation data
dat.v.nddlrd <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_v_nddlrd.csv"))

p.nddlrd <- predict(final, dat.v.nddlrd, re.form=NA, type="response") ## Predict model with validation data

roc_obj <- roc(dat.v.nddlrd$occ, p.nddlrd) ## AUC for NDD+LRD validation records
auc(roc_obj) ## 0.5117 So NDD+LRD sites are in totally different habitat to the NDMP sites

###########################################################
##### Predict suitability for all Devon and Dorset ########
###########################################################

##### ENVIRONMENT DATA #####
env.r <- stack(paste0(wd.env, "/env_100m7_dd"))
env.r <- env.r[[vars]]

### Need to add climate for all grid-cells, not just ones with distribution data
clim.vars <- c("rainfall_smr","tasmax_aut","tasmin_win","tasrng_win")

i <- 2019 ## The 5 year climate mean is calculated for the most recent period
s <- stack(paste0(wd.clim, "/",clim.vars,"_5yrmn_",i,".tif"))

s <- crop(s, extent(env.r))

nn <- names(s)[grepl("_5yrmn", names(s))]
n <- nn[1]

for (n in nn) {
  n2 <- substr(n, start=1, stop=nchar(n)-5) ## trim year from end
  env.r[[n2]] <- s[[n]]
}

writeRaster(env.r, file=paste0(wd.env, "/env_100m7_dd_fullclim_ndmp"))

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
raster::writeRaster(env.r, paste0(wd.env,"/env_100m7_dd_fullclim_ndmp_std")) ## Some variables have been added

env.r <- stack(paste0(wd.env,"/env_100m7_dd_fullclim_ndmp_std"))
  
### Make the prediction for the averaged model
p.ndmp <- predict(env.r, final, re.form=NA, type="response", full=T, se.fit=T)
writeRaster(p.ndmp, file=paste0(wd.out,"/final_preds_ndmp.tiff"))

### Standard error of the prediction
env.dat <- as.data.frame(extract(env.r, 1:ncell(env.r)))
p.ndmp.se <- stats::predict(final, newdata=env.dat, re.form=NA, type="response", full=T, se.fit=T)
test <- env.r[[1]]
values(test) <- NA
values(test) <- p.ndmp.se$se.fit
writeRaster(test, file=paste0(wd.out,"/final_preds_ndmp_se.tiff"))

##### Calculate thresholds #####

## Inspect the sensitivity and specificity
occ <- readOGR(dsn=wd.distn, layer="occ100m1990_dd") ## load, so always working with the same dataset

occ.ndmp <- occ
occ.ndmp$source[occ.ndmp$occ==0] <- "pseudo-absence"
occ.ndmp <- occ.ndmp[occ.ndmp$source!="LRD+NDD",] ## remove points with non-NDMP presences

p.pred <- raster::extract(p.ndmp, occ.ndmp)
r <- roc(occ.ndmp$occ, p.pred)
r # 0.8878

df <- data.frame(spec=r$specificities, sens=r$sensitivities, thresh=r$thresholds)

df[which.min(abs(0.95-df$sens)), ] ## The threshold that yields 95% sensitivity: 0.05963497
df[which.min(abs(df$spec-df$sens)), ] ## The threshold that minimises the difference between sensitivity and specificity: 0.2637337
df[which.max(df$spec+df$sens), ]## minimises difference between sensitivity and specificity 0.2633958

### Reclassify prediction raster according to two thresholds (from pa_model3_eval_ply)
## 95% sensitivity threshold
rcl <- matrix(c(0,0.05963497,0,   0.05963497,1,1), nrow=2, byrow=T) 
p95 <- reclassify(p.ndmp, rcl)
writeRaster(p95, file=paste0(wd.out,"/final_pred95_ndmp.tiff"))

## Threshold that minimises difference between sensitivity and specificity
rcl <- matrix(c(0,0.2637337,0,   0.2637337,1,1), nrow=2, byrow=T) 
pss <- reclassify(p.ndmp, rcl)
writeRaster(pss, file=paste0(wd.out,"/final_predss_ndmp.tiff"))

### Calculate sensitivity of predictions of LRD+NDD data at the different thresholds
occ.nddlrd <- occ
occ.nddlrd$source[occ.nddlrd$occ==0] <- "pseudo-absence"
occ.nddlrd <- occ.nddlrd[occ.nddlrd$source=="LRD+NDD",] ## remove points with absences and non-NDMP presences

p.pred.nddlrd.95 <- raster::extract(p95, occ.nddlrd)
sum(p.pred.nddlrd.95) ## 362 pf 525 presences are captured by this raster (69%, as opposed to 95% of the ndmp presences)

p.pred.nddlrd.ss <- raster::extract(pss, occ.nddlrd)
sum(p.pred.nddlrd.ss) ## 129 pf 525 presences are captured by this raster (25%)


########################################
##### Calculate R-squared for all models in best subset #####
########################################
dat.c <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_c_ndmp.csv"))

# load(paste0(wd.out,"/linear_quadratic_interactions_dredged_ndmp_combo1"))
# load(paste0(wd.out,"/linear_quadratic_interactions_dredged_ndmp_combo2"))
load(paste0(wd.out,"/linear_quadratic_interactions_dredged_ndmp_combo3"))

# mlqid.1 <- get.models(mlqid.1, subset=delta<2) ## mlqid.1 is an object containing the best models. 
# mlqid.2 <- get.models(mlqid.2, subset=delta<2) ## mlqid.2 is an object containing the best models. 
mlqid.3 <- get.models(mlqid.3, subset=delta<2) ## mlqid.3 is an object containing the best models. 

rsq.fn <- function(x) {
  (x$null.deviance - x$deviance) / x$null.deviance
}

# rsq1 <- unlist(lapply(mlqid.1, rsq.fn))
# rsq2 <- unlist(lapply(mlqid.2, rsq.fn))
rsq3 <- unlist(lapply(mlqid.3, rsq.fn))

mean(rsq3) ## 0.473862
sd(rsq3) ## 0.00404

