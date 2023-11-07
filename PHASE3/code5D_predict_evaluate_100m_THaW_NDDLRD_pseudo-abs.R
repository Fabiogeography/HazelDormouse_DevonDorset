##############################################################
##### PREDICT AND EVALUATE ########################
##### NDDLRD pseudo-abs #####
##### Written by: Regan Early ################################
##### Written on: 27th March 2021 ##############################
##### Modified on: 22nd March 2022  #########################
##############################################################

# .libPaths("C:\\Rpackages")
.libPaths("D:/SOFTWARE/R-4.1.2/library")
library(dplyr) ## bind_rows, sample_frac
library(raster)
library(rgdal)
library(pROC)
library(MuMIn) # dredge() stdize() model.sel()
library(dismo) ## MESS with rasters

wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/100m/"
wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
wd.clim <- "E:/NON_PROJECT/DORMOUSE_WALES/GIS/100m" # "E:/GIS_DATA/CLIMATE/UK/HadUK-Grid_1km"
wd.distn <- "E:/NON_PROJECT/DORMOUSE_NE/DISTRIBUTION"
# wd.out <- "~/re259/UoE_U_Drive/DORMOUSE_NE/SDM/100m"
# wd.env <- "~/re259/UoE_U_Drive/DORMOUSE_NE/SDM/100m"

load(file=paste0(wd.out,"/mods_final_t_sampled_nddlrd_pabs_frag"))

destdize <- read.csv(paste0(wd.out, "/dat_means_sds_23Mar2022_t_frag_sampled_prepped.csv"))

# w <- readOGR(dsn="E:/GIS_DATA/UK/bdline_essh_gb/Data/Wales", layer="Wales_mainland")

###########################################################################################
##### Test robustness of model with the same variables calibrated with different data #####
###########################################################################################
### Predict the original validation data
dat.v <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_frag_v_sampled_nddlrd_pabs.csv"))
dat.m <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_frag_m_sampled_nddlrd_pabs.csv"))

vars <- as.vector(colnames(final$coefficients))
vars <- vars[!grepl("Intercept|I|:", vars)] ## remove intercept, quadratic, and interaction terms from the variables list

dat <- read.csv(paste0(wd.env,"/dat_23Mar2022_t_frag_sampled_prepped.csv"), as.is=T) 
# dat <- read.csv(paste0(wd.env,"/dat_30Nov2021_dd.csv"))
dat <- as.data.frame(cbind(dat[,1:8], dat[,vars])) 

p <- predict(final, dat.v, re.form=NA, type="response") ## Predict model with validation data

roc_obj <- roc(dat.v$occ, p, plot=T)
auc(roc_obj) ## 0.7224

###########################################################
##### Predict suitability for all THaW region ########
###########################################################

##### ENVIRONMENT DATA #####
env.r <- stack(paste0(wd.env, "/env_100m7_t_sampled_frag")) 
# env.r <- stack(paste0(wd.env, "/env_100m7_dd"))
env.r <- env.r[[vars]]

### Need to add climate for all grid-cells, not just ones with distribution data
clim.vars <- c("rainfall_spr","sun_spr","tasmax_spr")

i <- 2019 ## The 5 year climate mean is calculated for the most recent period
s <- stack(paste0(wd.clim, "/",clim.vars,"_5yrmn_",i,".tif"))

s <- crop(s, extent(env.r))

nn <- names(s)[grepl("_5yrmn", names(s))]
n <- nn[1]

for (n in nn) {
  n2 <- substr(n, start=1, stop=nchar(n)-5) ## trim year from end
  env.r[[n2]] <- s[[n]]
} ## These values are higher than the ones in teh raster I originally wrote on 7th Dec, even though env_100m7_dd was written on 1st Dec

writeRaster(env.r, file=paste0(wd.env, "/env_100m7_t_sampled_nddlrd_pabs_fullclim"))  ## )

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

raster::writeRaster(env.r, paste0(wd.env,"/env_100m7_t_sampled_nddlrd_pabs_fullclim_std")) ## Some variables have been added

### Make the prediction for the averaged model
p <- predict(env.r, final, re.form=NA, type="response", full=T, se.fit=T)
writeRaster(p, file=paste0(wd.out,"/final_preds_t_sampled_nddlrd_pabs.tiff"))
# p <- raster(paste0(wd.out,"final_preds_t_sampled.tiff"))

### Standard error of the prediction
env.dat <- as.data.frame(extract(env.r, 1:ncell(env.r)))
p.se <- stats::predict(final, newdata=env.dat, re.form=NA, type="response", full=T, se.fit=T)
test <- env.r[[1]]
values(test) <- NA
values(test) <- p.se$se.fit
writeRaster(test, file=paste0(wd.out,"/final_preds_t_sampled_nddlrd_pabs_se.tiff"))

# ##### MESS #####
# occ <- readOGR(dsn=wd.distn, layer="occ100m1990_dd") ## load, so always working with the same dataset
occ <- readOGR(dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="occ100mt_sampled_dd") ## load, so always working with the same dataset. Created in code3
# env.pa.mat <- extract(env.r, occ)
# 
# mess.out <- mess(x=env.r, v=env.pa.mat, full=TRUE)
# names(mess.out) <- c(vars, "rmess")
# plot(mess.out[["rmess"]])
# writeRaster(mess.out$rmess, "rmess") ## the MESS surface for all variables
# 
# for(i in vars) {
#   writeRaster(mess.out[[i]], paste0(pathway, i, "_mess.tif")) ## the MESS surfaces for all individual variables
# }


##### Calculate thresholds #####

## Inspect the sensitivity and specificity
occ.nddlrd <- occ
# occ.nddlrd$source[occ.nddlrd$occ==0] <- "pseudo-absence"
occ.nddlrd <- occ.nddlrd[!(occ.nddlrd$source %in% c("NDMP", "Non-PTES", "PTES")),] ## remove points with non-NDD+LRD presences and with survey absences

p.pred <- raster::extract(p, occ.nddlrd)
r <- roc(occ.nddlrd$occ, p.pred)
r #  0.824

df <- data.frame(spec=r$specificities, sens=r$sensitivities, thresh=r$thresholds)

df[which.min(abs(0.95-df$sens)), ] ## The threshold that yields 95% sensitivity: 0.2710008
df[which.min(abs(df$spec-df$sens)), ] ## The threshold that minimises the difference between sensitivity and specificity: 0.66068
df[which.max(df$spec+df$sens), ]## minimises difference between sensitivity and specificity 0.6438594

### Reclassify prediction raster according to two thresholds (from pa_model3_eval_ply)
## 95% sensitivity threshold
rcl <- matrix(c(0,0.2710008,0,   0.2710008,1,1), nrow=2, byrow=T) 
p95 <- reclassify(p, rcl)
writeRaster(p95, file=paste0(wd.out,"/final_pred95_t_sampled_nddlrd_pabs.tiff"))
p.pred95 <- raster::extract(p95, occ.nddlrd)
table(p.pred95)
113/(113+384) ## overpred

## Threshold that minimises difference between sensitivity and specificity
rcl <- matrix(c(0,0.6438594,0,   0.6438594,1,1), nrow=2, byrow=T) 
pss <- reclassify(p, rcl)
writeRaster(pss, file=paste0(wd.out,"/final_predss_t_sampled_nddlrd_pabs.tiff"))

########################################
##### Calculate R-squared for all models in best subset #####
########################################
dat.c <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_frag_c_sampled_nddlrd_pabs.csv"))

load(paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs1"))
load(paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs2"))
load(paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs3"))

mlqid.1 <- get.models(mlqid.1, subset=delta<2) ## mlqid.1 is an object containing the best models. 
mlqid.2 <- get.models(mlqid.2, subset=delta<2) ## mlqid.2 is an object containing the best models. 
mlqid.3 <- get.models(mlqid.3, subset=delta<2) ## mlqid.2 is an object containing the best models. 

rsq.fn <- function(x) {
  (x$null.deviance - x$deviance) / x$null.deviance
}

rsq1 <- unlist(lapply(mlqid.1, rsq.fn))
rsq2 <- unlist(lapply(mlqid.2, rsq.fn))
rsq3 <- unlist(lapply(mlqid.3, rsq.fn))

mean(c(rsq1, rsq2, rsq3)) ## 0.2323364
sd(c(rsq1, rsq2, rsq3)) ## 0.01137785

