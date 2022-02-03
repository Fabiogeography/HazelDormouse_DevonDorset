##############################################################
##### SDM ########################
##### Written by: Regan Early ################################
##### Written on: 22nd March 2021 ##############################
##### Modified on: 30th November 2021  #########################
##############################################################

# library(ggplot2)#, lib.loc="D:/SOFTWARE/R-4.1.1/library") ## Needed to make interact_plot run and needs to be version for R 3.3.0 or later
# library(interactions) ## interact_plot
# library(Hmisc) ## rcorr ## load after ggplot2 so it doesn't inadvertently load the wrong version of the package
library(dplyr, lib.loc="D:/SOFTWARE/R-4.1.1/library") ## bind_rows, sample_frac
# library(corrplot)
# library(raster)
# library(rgdal)
# library(spatialEco)
# library(sp)
# library(pROC)
# library(tidyverse) ## glimpse and useful functions
# library(tidyr) ## gather
library(MuMIn, lib.loc="D:/SOFTWARE/R-4.1.1/library") # dredge() stdize() model.sel()
library(car, lib.loc="D:/SOFTWARE/R-4.1.1/library") ## Variance Inflation Factor analysis
# library(snow) ## for running on multiple cores of a cluster. Used for dredge function (best when used on a server rather than a desktop).

wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/100m"
wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"

##### MAKE MULTIVARIATE MODELS #####
### Variables to be attempted - see word doc
vars <- c("broadl_all_m", "broadl_c", "conif_all_m", "conif_c", 
          "lcm_broadl_m", "lcm_broadl_c", "lcm_conif_m", "lcm_conif_c", "anc_wood",
          "rainfall_smr_5yrmn",	"sun_aut_5yrmn", "tasmin_win_5yrmn", "tasmax_aut_5yrmn", "tasrng_win_5yrmn",
          "OS_Terrain_100_NS", "OS_Terrain_100_WE", "OS_Terrain_100_slope_pct")
dat <- read.csv(paste0(wd.env,"/dat_30Nov2021_dd.csv"))
dat <- as.data.frame(cbind(dat[,1:7], dat[,vars]))
# dat$easting[is.na(dat$easting)] <- dat$northing[is.na(dat$northing)] <- dat$m100[is.na(dat$m100)] <- -9999

##### Scale explanatory variables as large numerical scope of variables can cause estimation problems (https://r-sig-mixed-models.r-project.narkive.com/fXcaHABA/r-sig-me-cholmod-warning-with-glmer). #####
## First save the values used to stdize, so as to graph and map later
means <- apply(na.omit(dat[,vars]), 2, mean)
sds <- apply(na.omit(dat[,vars]), 2, sd)
destdize <- rbind(means, sds)
# write.csv(destdize, paste0(wd.out, "/dat_means_sds_30Nov2021_dd.csv"), row.names=F)
# destdize <- read.csv(paste0(wd.out, "/dat_means_sds_30Nov2021_dd.csv"))

## Use method stdize from MuMin, which subtracts the mean and then divides by the standard deviation
dat.m <- na.omit(dat[,c("occ", vars)]) ## NAs have to be removed for dredge function
dat.m[,vars ] <- stdize(dat.m[,vars])
dat.m <- as.data.frame(cbind(dat.m, id=c(1:nrow(dat.m))))
# write.csv(dat.m, paste0(wd.out, "/dat_std_30Nov2021_dd_m.csv"), row.names=F)
# dat.m <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_m.csv"))

### Calibration and validation data to make the model
dat.c <- sample_frac(cbind(dat.m), size=0.7)
# write.csv(dat.c, paste0(wd.out, "/dat_std_30Nov2021_dd_c.csv"), row.names=F)
# dat.c <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_c.csv"))

dat.v <- dat.m[!(dat.m$id %in% dat.c$id),]
# write.csv(dat.v, paste0(wd.out, "/dat_std_30Nov2021_dd_v.csv"), row.names=F)
# dat.v <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_v.csv"))

##### Model with linear terms only #####
vars.m <- paste(vars, collapse=" + ")
ml <- glm(as.formula(paste0("occ ~ ", vars.m)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data
## A model to be dredged must not have na.action as na.omit or na.exclude, or dredge will fail.

## Check for variance inflation - an effect of collinear explanatory variables
vif(ml, singular.ok = TRUE) ## all ok. broadl and conif borderline

## Dredge model with linear terms. Previous version only permitted one of the two highly correlated terms
mld.all <- dredge(ml, beta="sd", evaluate=T, trace=2)
save(mld.all, file=paste0(wd.out, "/linear_dredged_all"))
# load(paste0(wd.out, "/linear_dredged_all"))
write.csv(mld, file=paste0(wd.out,"/linear_dredged_all.csv"), row.names=F)
mld.all[mld.all$delta<=2,] ## 20 models fall into the best model subset. 
# m.best <- get.models(md, 1)[[1]]
# mld.all <- get.models(mld.all, subset=delta<2) ## mld is now an object containing the ten best models. 

mld.sub2 <- dredge(ml, beta="sd", evaluate=T, trace=2) ## , subset = !(lcm_conif_m && conif_all_m) & !(lcm_broadl_m && broadl_all_m) 
                  # & !(lcm_conif_c && conif_c) & !(lcm_broadl_c && broadl_c) ) ## Not permitting variables to be in the same model actually led to a much broader list of best mofdels, with poorer AIC than when variables could be considered together. 
## Even if subset not enforced, the best model didn't contain the mean values of nfi and lcm broadleaved and coniferous forest together. NFI broadleaved forest not selected
save(mld.sub, file=paste0(wd.out, "/linear_dredged_sub"))
# load(paste0(wd.out, "/linear_dredged_sub"))
write.csv(mld.sub, file=paste0(wd.out,"/linear_dredged_sub.csv"), row.names=F)
mld.sub2[mld.sub2$delta<=2,] ## 10 models fall into the best model subset. The best model has a delta AIC of -0.05, but the second -0.86. So the first two are indistinguishable, one has nfi broadleaf mean, ant the other doesn't . Higher ranked models less likely to contain nfi forest variables (almost all contain lcm)
# m.best <- get.models(md, 1)[[1]]
# mld.sub <- get.models(mld.sub, subset=delta<2) ## mld is now an object containing the ten best models. 

##### Model with linear and quadratic terms #####
### Get the linear terms that were selected in any of the best model subset and make into quadratic terms
# vars.mld <- unique(unlist(lapply(mld , function(x) {rownames(summary(x)$coefficients)})))
# ## Actually just do this for the best model, in the interests of time (dredging the larger model takes ages)
# ## The other five include a single extra variable each with a substantially smaller coefficient than the other variables,
# ## which are retained in all of the best model subset.
# ## Additionally the coefficient estimats of the always-retained models are very little altered by the inclusion of the single extra variable 
# vars.mld <- rownames(summary(mld[[1]])$coefficients)
# vars.mld <- vars.mld[vars.mld!="(Intercept)"]

vars.mq <- paste0("I(", vars, "^2)") ## quadratic terms

vars.mlq <- paste(c(vars, vars.mq), collapse=" + ") ## all linear and quadratic terms. was vars.mld

### Dredge the model with all linear terms retained in the best model subset, and all of their quadratic terms
mlq <- glm(as.formula(paste0("occ ~ ", vars.mlq)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data
vif(mlq) ## ok

## Retain quadratic term only if linear term is added
## Done manually, but should be done automatically
msubset <- expression(dc("broadl_all_m", "I(broadl_all_m^2)") &
                        dc("broadl_c", "I(broadl_c^2)") &
                        dc("conif_all_m", "I(conif_all_m^2)") &
                        dc("conif_c", "I(conif_c^2)") &
                        dc("lcm_broadl_m", "I(lcm_broadl_m^2)") &
                        dc("lcm_broadl_c", "I(lcm_broadl_c^2)") &
                        dc("lcm_conif_m", "I(lcm_conif_m^2)") &
                        dc("lcm_conif_c", "I(lcm_conif_c^2)") &
                        dc("rainfall_aut_5yrmn", "I(rainfall_aut_5yrmn^2)") &
                        dc("tasmin_win_5yrmn", "I(tasmin_win_5yrmn^2)") &
                        dc("tasmax_aut_5yrmn", "I(tasmax_aut_5yrmn^2)") &
                        dc("tasrng_win_5yrmn", "I(tasrng_win_5yrmn^2)") &
                        dc("sun_aut_5yrmn", "I(sun_aut_5yrmn^2)") &
                        dc("OS_Terrain_100_NS", "I(OS_Terrain_100_NS^2)") &
                        dc("OS_Terrain_100_WE", "I(OS_Terrain_100_WE^2)") &
                        dc("OS_Terrain_100_slope_pct", "I(OS_Terrain_100_slope_pct^2)"))
  
# ## Set up the cluster. Number of cores limited by computing resources - best on server
# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 10), type = clusterType))
# clusterExport(clust, "dat.c") ## Send the data to each node of the cluster
# clusterCall(clust, function() library(MuMIn)) ## Call the function to load the library on each node of the cluster

mlqd <- pdredge(mlq, beta="sd", evaluate=T, trace=2, subset = msubset)#, cluster=clust) 

save(mlqd, file=paste0(wd.out, "/linear_quadratic_dredged"))
write.csv(mlqd, file=paste0(wd.out,"/linear_quadratic_dredged.csv"), row.names=F)

## Best subset of models with delta AICc <=2
mlqd[mlqd$delta<=2,] ## 11 models in the best model subset. Two quadratic forms reliably included: coniferous and tasmax_smr

mlqd <- get.models(mlqd, subset=delta<2) ## mlqd is now an object containing the three best models. 

## Check for variance inflation - an effect of collinear explanatory variables
for(i in 1:length(mlqd)) {print(vif(mlqd[[i]]), singular.ok = TRUE)} ## tasmax_aut has high VIF, but seems to be with its own quadratic form. shrub is just over the limit of 5.

# Inspect models (ask whether fixed effect relationships meaningful)
for(i in 1:length(mlqd)) {print(summary(mlqd[[i]]))}

##### Model with interactions between linear terms #####
## Create list of interactions to test
vars.mlqd <- unique(unlist(lapply(mlqd , function(x) {rownames(summary(x)$coefficients)})))
vars.mlqd <- vars.mlqd[vars.mlqd!="(Intercept)"]

vars.mi <- combn(subset(vars.mlqd, !grepl("2", vars.mlqd)),2)
vars.mi <- apply(vars.mi, MARGIN=2, FUN=function(x) {paste(x, collapse="*")})

### Group interactions into 3 subsets, otherwise there are too many fixed terms for model selection to work
vars.mi.1 <- paste0("(", paste(vars.mi[1:(length(vars.mi)/3)], collapse=") + ("),")")
start <- 1+(length(vars.mi)/3); stop <- (length(vars.mi)/3)*2
vars.mi.2 <- paste0("(", paste(vars.mi[start : stop], collapse=") + ("),")")
start <- 1+(length(vars.mi)/2); stop <- length(vars.mi)
vars.mi.3 <- paste0("(", paste(vars.mi[start : stop], collapse=") + ("),")")

mlqi.full1 <- glm(as.formula(paste0("occ ~ ", paste(c(vars.mlqd, vars.mi.1), collapse=" + "))), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data
mlqi.full2 <- glm(as.formula(paste0("occ ~ ", paste(c(vars.mlqd, vars.mi.2), collapse=" + "))), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data
mlqi.full3 <- glm(as.formula(paste0("occ ~ ", paste(c(vars.mlqd, vars.mi.3), collapse=" + "))), family=binomial, dat.c, na.action=na.fail) ##  

# Compare models using AICc and save selection table
# mlqid.1 <- pdredge(mlqi.full1, beta="sd", evaluate=T, trace=2, subset = msubset)#, cluster=clust) ## 
# save(mlqid.1, file=paste0(wd.out, "/linear_quadratic_interactions_dredged1"))
# write.csv(mlqid.1, file=paste0(wd.out,"/linear_quadratic_interactions_dredged1.csv"), row.names=F)

mlqid.2 <- pdredge(mlqi.full2, beta="sd", evaluate=T, trace=2, subset = msubset)#, cluster=clust)
save(mlqid.2, file=paste0(wd.out, "/linear_quadratic_interactions_dredged2"))
write.csv(mlqid.2, file=paste0(wd.out,"/linear_quadratic_interactions_dredged2.csv"), row.names=F)

mlqid.3 <- pdredge(mlqi.full3, beta="sd", evaluate=T, trace=2, subset = msubset)#, cluster=clust) ## Fails. Likely due to starting models being unfeasible
save(mlqid.3, file=paste0(wd.out, "/linear_quadratic_interactions_dredged3"))
write.csv(mlqid.3, file=paste0(wd.out,"/linear_quadratic_interactions_dredged3.csv"), row.names=F)

## Best subset of models with delta AICc <=2. Notes about interactions included are not correct
mlqid.1[mlqid.1$delta<=2,] ## 42 models in the best model subset. Best AICc= 1115.4 
 ## brd_all_m:tsmn_win_5yr (-ve), brd_all_m:cnf_all_m (-ve), broadl_all:tasmin_win (-ve), conif_all_m:rainfall_aut (+ve), and conif_all_m:tasmax_aut (-ve) are in all of the best models.
mlqid.2[mlqid.2$delta<=2,] ## 24 models in the best model subset. Best AICc = 1113.5 
 ## cnf_all_m:rnf_aut_5yr (+ve), cnf_all_m:tsmx_aut_5yr (-ve), cpp_m:fll_c (-ve) and cpp_m:fll_m (+ve) are in all of the best models. 
mlqid.3[mlqid.3$delta<=2,] ## 33 models in the best model subset. Best AIC=1114.5

## Combine all three model selection objects
mlqid <- merge(mlqid.1[mlqid.1$delta<=2,], mlqid.2[mlqid.2$delta<=2,], mlqid.3[mlqid.3$delta<=2,])
mlqid <- mlqid[-22,] ## There are duplicate models
mlqid[mlqid$delta<=2,]

save(mlqid, file=paste0(wd.out, "/linear_quadratic_interactions_dredged_best"))
write.csv(mlqid, file=paste0(wd.out,"/linear_quadratic_interactions_dredged_best.csv"), row.names=F)

#################################################
###### Produce & save final model ###############
#################################################
load(file=paste0(wd.out, "/linear_quadratic_interactions_dredged_best"))

final <- model.avg(mlqid, subset = delta<2, fit=T) ## Average the models with delta AIC < value and ensure the models are actually fit

### Output the coefficients in a table
final.coefs <- final$coefficients["full",] ## Average parameter estimates over all models (value will be 0 in models wehre the parameter wasn't included)
vars <- names(final.coefs)
coefs <- as.vector(final.coefs)
final.coefs <- data.frame(Variables=vars, MeanCoefficients=coefs)
rownames(final.coefs) <- vars

## Write outputs
save(final, file=paste0(wd.out,"/mods_final"))
write.csv(final.coefs, file=paste0(wd.out,"/mods_final.csv"), row.names=F)

