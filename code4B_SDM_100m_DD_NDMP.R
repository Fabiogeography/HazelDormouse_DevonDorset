##############################################################
##### SDM ########################
##### Written by: Regan Early ################################
##### Written on: 22nd March 2021 ##############################
##### Modified on: 7th December 2021  #########################
##############################################################

library(dplyr, lib.loc="D:/SOFTWARE/R-4.1.1/library") ## bind_rows, sample_frac
library(MuMIn, lib.loc="D:/SOFTWARE/R-4.1.1/library") # dredge() stdize() model.sel()
library(car, lib.loc="D:/SOFTWARE/R-4.1.1/library") ## Variance Inflation Factor analysis
library(snow) ## for running on multiple cores of a cluster. Used for dredge function (best when used on a server rather than a desktop).

# wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/100m"
# wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
wd.out <- "~/re259/UoE_U_Drive/DORMOUSE_NE/SDM/100m"
wd.env <- "~/re259/UoE_U_Drive/DORMOUSE_NE/SDM/100m"

##### MAKE MULTIVARIATE MODELS #####
### Variables to be attempted - see word doc
vars <- c("combo_broadl_m", "combo_broadl_c", "combo_conif_m", "combo_conif_c", "anc_wood",
          "rainfall_smr_5yrmn",	"sun_aut_5yrmn", "tasmin_win_5yrmn", "tasmax_aut_5yrmn", "tasrng_win_5yrmn",
          "OS_Terrain_100_NS", "OS_Terrain_100_WE", "OS_Terrain_100_slope_pct")
dat <- read.csv(paste0(wd.env,"/dat_30Nov2021_dd.csv"), as.is=T) 
dat <- as.data.frame(cbind(dat[,1:7], dat[,vars]))
dat$source[dat$occ==0] <- "pseudo-absence"
dat <- dat[dat$source!="LRD+NDD",] ## remove points with non-NDMP presences

##### Scale explanatory variables as large numerical scope of variables can cause estimation problems (https://r-sig-mixed-models.r-project.narkive.com/fXcaHABA/r-sig-me-cholmod-warning-with-glmer). #####
## First save the values used to stdize, so as to graph and map later
means <- apply(na.omit(dat[,vars]), 2, mean)
sds <- apply(na.omit(dat[,vars]), 2, sd)
destdize <- rbind(means, sds)
# write.csv(destdize, paste0(wd.out, "/dat_means_sds_30Nov2021_dd_ndmp.csv"), row.names=F)
# destdize <- read.csv(paste0(wd.out, "/dat_means_sds_30Nov2021_dd_ndmp.csv"))

## Use method stdize from MuMin, which subtracts the mean and then divides by the standard deviation
dat.m <- na.omit(dat[,c("occ", vars)]) ## NAs have to be removed for dredge function
dat.m[,vars ] <- stdize(dat.m[,vars])
dat.m <- as.data.frame(cbind(dat.m, id=c(1:nrow(dat.m))))
# write.csv(dat.m, paste0(wd.out, "/dat_std_30Nov2021_dd_m_ndmp.csv"), row.names=F)
# dat.m <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_m_ndmp.csv"))

### Calibration and validation data to make the model
dat.c <- sample_frac(cbind(dat.m), size=0.7)
# write.csv(dat.c, paste0(wd.out, "/dat_std_30Nov2021_dd_c_ndmp.csv"), row.names=F)
# dat.c <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_c_ndmp.csv"))

dat.v <- dat.m[!(dat.m$id %in% dat.c$id),]
# write.csv(dat.v, paste0(wd.out, "/dat_std_30Nov2021_dd_v_ndmp.csv"), row.names=F)
# dat.v <- read.csv(paste0(wd.out, "/dat_std_30Nov2021_dd_v_ndmp.csv"))

##### Model with linear terms only #####
vars.m <- paste(vars, collapse=" + ")
ml <- glm(as.formula(paste0("occ ~ ", vars.m)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. 
## A model to be dredged must not have na.action as na.omit or na.exclude, or dredge will fail.

## Check for variance inflation - an effect of collinear explanatory variables
vif(ml, singular.ok = TRUE) ## all fine

## Set up the cluster. Number of cores limited by computing resources - best on server
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 10), type = clusterType))
clusterExport(clust, "dat.c") ## Send the data to each node of the cluster
clusterCall(clust, function() library(MuMIn)) ## Call the function to load the library on each node of the cluster

## Dredge model with linear terms. Previous version only permitted one of the two highly correlated terms
mld.all <- pdredge(ml, beta="sd", evaluate=T, trace=2, cluster=clust)
save(mld.all, file=paste0(wd.out, "/linear_dredged_all_ndmp_combo"))
# load(paste0(wd.out, "/linear_dredged_all"))
write.csv(mld.all, file=paste0(wd.out,"/linear_dredged_all_ndmp_combo.csv"), row.names=F)
mld.all[mld.all$delta<=2,] ## Four models fall into the best model subset. Broadleaf mean and change, conifer change, slope, rainfall, tasmax, tasmin, tasrange arte the consistent variables
# m.best <- get.models(md, 1)[[1]]
mld.all <- get.models(mld.all, subset=delta<2) ## mld is now an object containing the best models.

# mld.sub2 <- dredge(ml, beta="sd", evaluate=T, trace=2, subset = !(lcm_conif_m && conif_all_m) & !(lcm_broadl_m && broadl_all_m) 
#                   & !(lcm_conif_c && conif_c) & !(lcm_broadl_c && broadl_c))
# save(mld.sub, file=paste0(wd.out, "/linear_dredged_sub_ndmp"))
# load(paste0(wd.out, "/linear_dredged_sub"))
# write.csv(mld.sub, file=paste0(wd.out,"/linear_dredged_sub_ndmp.csv"), row.names=F)
# mld.sub2[mld.sub2$delta<=2,] ## 10 models fall into the best model subset. The best model has a delta AIC of -0.05, but the second -0.86. So the first two are indistinguishable, one has nfi broadleaf mean, ant the other doesn't . Higher ranked models less likely to contain nfi forest variables (almost all contain lcm)
# m.best <- get.models(md, 1)[[1]]
# mld.sub <- get.models(mld.sub, subset=delta<2) ## mld is now an object containing the ten best models. 

##### Model with linear and quadratic terms #####
### Get the linear terms that were selected in any of the best model subset and make into quadratic terms
vars.mld <- unique(unlist(lapply(mld.all, function(x) {rownames(summary(x)$coefficients)})))
# ## Actually just do this for the best model, in the interests of time (dredging the larger model takes ages)
# ## The other five include a single extra variable each with a substantially smaller coefficient than the other variables,
# ## which are retained in all of the best model subset.
# ## Additionally the coefficient estimats of the always-retained models are very little altered by the inclusion of the single extra variable 
# vars.mld <- rownames(summary(mld[[1]])$coefficients)
vars.mld <- vars.mld[vars.mld!="(Intercept)"]

vars.mq <- paste0("I(", vars.mld, "^2)") ## quadratic terms
vars.mlq <- paste(c(vars, vars.mq), collapse=" + ") ## all linear and quadratic terms selected in the linear best model subset. 11 variables
## Do manually so only have 31 variables. Tasmin and WE aspect arehardly represented in best model subset for linear terms.
# vars.mlq <- "broadl_all_m + broadl_c + conif_all_m + conif_c + lcm_broadl_m + lcm_broadl_c + lcm_conif_m + lcm_conif_c + anc_wood + rainfall_smr_5yrmn + sun_aut_5yrmn + tasmax_aut_5yrmn + tasrng_win_5yrmn + OS_Terrain_100_NS + OS_Terrain_100_slope_pct + I(anc_wood^2) + I(lcm_broadl_m^2) + I(OS_Terrain_100_NS^2) + I(rainfall_smr_5yrmn^2) + I(sun_aut_5yrmn^2) + I(tasmax_aut_5yrmn^2) + I(OS_Terrain_100_slope_pct^2) + I(tasrng_win_5yrmn^2) + I(broadl_c^2) + I(lcm_conif_c^2) + I(conif_c^2) + I(lcm_broadl_c^2) + I(conif_all_m^2) + I(broadl_all_m^2)"
## OR
# vars.mlq <- paste0("poly(", vars.mld, ", 2)", collapse=" + ") ## alternative coding of quadratics, removing collinearity

### Dredge the model with all linear terms retained in the best model subset, and all of their quadratic terms
mlq <- glm(as.formula(paste0("occ ~ ", vars.mlq)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. Fitted probabilities 0 or 1 occurred.
vif(mlq) ## ok

## Retain quadratic term only if linear term is added
## Done manually, but should be done automatically
msubset <- expression(dc("combo_broadl_m", "I(combo_broadl_m^2)") &
                        dc("combo_conif_m", "I(combo_conif_m^2)") &
                        dc("combo_broadl_c", "I(combo_broadl_c^2)") &
                        dc("combo_conif_c", "I(combo_conif_c^2)") &
dc("rainfall_smr_5yrmn", "I(rainfall_smr_5yrmn^2)") &
dc("tasmin_win_5yrmn", "I(tasmin_win_5yrmn^2)") &
dc("tasmax_aut_5yrmn", "I(tasmax_aut_5yrmn^2)") &
dc("tasrng_win_5yrmn", "I(tasrng_win_5yrmn^2)") &
dc("sun_aut_5yrmn", "I(sun_aut_5yrmn^2)") &
dc("OS_Terrain_100_NS", "I(OS_Terrain_100_NS^2)") &
dc("OS_Terrain_100_WE", "I(OS_Terrain_100_WE^2)") &
dc("OS_Terrain_100_slope_pct", "I(OS_Terrain_100_slope_pct^2)"))
                        
mlqd <- pdredge(mlq, beta="sd", evaluate=T, trace=2, subset=msubset, cluster=clust)

save(mlqd, file=paste0(wd.out, "/linear_quadratic_dredged_ndmp_combo"))
write.csv(mlqd, file=paste0(wd.out,"/linear_quadratic_dredged_ndmp_combo.csv"), row.names=F)

## Best subset of models with delta AICc <=2
mlqd[mlqd$delta<=2,] ## 68 models in the best model subset. 24 variables retained in at least one model (all of them)
## Ancient woodland was kiocked out in teh linear stage 
## No models have a strikingly lower AIC than others 
## Need to look for interactions. Only test for interactions amongst variables in every model?
## Linear effects selected in every model: combo_broad_c combo_broad_m, combo_conif_c, slope, rainfall_smr, tasmax_autumn, tasmin_win, tasrange_win. WE selected in majority.

# mlqd <- get.models(mlqd, subset=delta<2) ## mlqd is now an object containing the three best models.

## Check for variance inflation - an effect of collinear explanatory variables
for(i in 1:length(mlqd)) {print(vif(mlqd[[i]]), singular.ok = TRUE)} ## ancient woodland has high VIF (~14), but seems to be with its own quadratic form, as only happens in models with quadratic.

# Inspect models (ask whether fixed effect relationships meaningful)
for(i in 1:length(mlqd)) {print(summary(mlqd[[i]]))}

##### Model with interactions between linear terms #####
## Create list of interactions to test
vars.mlqd2 <- unique(unlist(lapply(mlqd, function(x) {rownames(summary(x)$coefficients)})))
vars.mlqd2 <- vars.mlqd2[!grepl("Intercept|I|:", vars.mlqd2)] ## remove intercept, quadratic, and interaction terms from the variables

# vars.mi <- as.data.frame(combn(subset(vars.mlqd2, !grepl("2", vars.mlqd2)),2)) ## This is the old method, selects every linear effect in the best models.
# vars.mi <- apply(vars.mi, MARGIN=2, FUN=function(x) {paste(x, collapse="*")}) ## 

vars.mi <- combn(c("combo_broadl_m", "combo_broadl_c", "combo_conif_c", "OS_Terrain_100_slope_pct", "rainfall_smr_5yrmn", "tasmax_aut_5yrmn", "tasmin_win_5yrmn", "tasrng_win_5yrmn"), 2) ## New method. Keeps only linear terms selecged in every model.
vars.mi <- apply(vars.mi, MARGIN=2, FUN=function(x) {paste(x, collapse="*")})

## Group interactions into 2 subsets, otherwise there are too many fixed terms for model selection to work
vars.mi.1 <- paste(c(vars.mlqd2, vars.mi[1:9]), collapse= " + ")
vars.mi.2 <- paste(c(vars.mlqd2, vars.mi[10:18]), collapse= " + ")
vars.mi.3 <- paste(c(vars.mlqd2, vars.mi[19:28]), collapse= " + ")

## Old method
# vars.mi.1 <- paste0("(", paste(vars.mi[1:(length(vars.mi)/2)], collapse=") + ("),")")
# start <- 1+(length(vars.mi)/3); stop <- (length(vars.mi)/3)*2
# vars.mi.2 <- paste0("(", paste(vars.mi[start : stop], collapse=") + ("),")")
# start <- 1+(length(vars.mi)/2); stop <- length(vars.mi)
# vars.mi.3 <- paste0("(", paste(vars.mi[start : stop], collapse=") + ("),")")

mlqi.full1 <- glm(as.formula(paste0("occ ~ ", paste(c(vars.mlqd2, vars.mi.1), collapse=" + "))), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. 28 variables
mlqi.full2 <- glm(as.formula(paste0("occ ~ ", paste(c(vars.mlqd2, vars.mi.2), collapse=" + "))), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data
mlqi.full3 <- glm(as.formula(paste0("occ ~ ", paste(c(vars.mlqd2, vars.mi.3), collapse=" + "))), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data

mlqi.full <- glm(as.formula(paste0("occ ~ ", paste(c(vars.mlqd2, vars.mi.1, vars.mi.2, vars.mi.3), collapse=" + "))), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. 28 variables

# Compare models using AICc and save selection table
mlqid.1 <- pdredge(mlqi.full1, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust) ##
save(mlqid.1, file=paste0(wd.out, "/linear_quadratic_interactions_dredged_ndmp_combo1"))
write.csv(mlqid.1, file=paste0(wd.out,"/linear_quadratic_interactions_dredged_ndmp_combo1.csv"), row.names=F)

mlqid.2 <- pdredge(mlqi.full2, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
save(mlqid.2, file=paste0(wd.out, "/linear_quadratic_interactions_dredged_ndmp_combo2"))
write.csv(mlqid.2, file=paste0(wd.out,"/linear_quadratic_interactions_dredged_ndmp_combo2.csv"), row.names=F)

mlqid.3 <- pdredge(mlqi.full3, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust) ##
save(mlqid.3, file=paste0(wd.out, "/linear_quadratic_interactions_dredged_ndmp_combo3"))
write.csv(mlqid.3, file=paste0(wd.out,"/linear_quadratic_interactions_dredged_ndmp_combo3.csv"), row.names=F)

## Best subset of models with delta AICc <=2. 
d1 <- mlqid.1[mlqid.1$delta<=2,] ## 41 models in the best model subset. 20 starting variables, 20 in best subset.
d2 <- mlqid.2[mlqid.2$delta<=2,] ## 13 models in best subset. 20 starting variables, 20 in best subset.
d3 <- mlqid.3[mlqid.3$delta<=2,] ## 12 models in best subset. 21 starting variables, 21 in best subset.

d1$AICc
d2$AICc ## All lower than d1
d3$AICc ## All notably lower than d1 and d2.
save(mlqid.3, file=paste0(wd.out, "/linear_quadratic_interactions_dredged_ndmp_best"))
write.csv(mlqid.3, file=paste0(wd.out,"/linear_quadratic_interactions_dredged_ndmp_best.csv"), row.names=F)

## Get the terms in the best models
#################################################
###### Produce & save final model ###############
#################################################
load(file=paste0(wd.out, "/linear_quadratic_interactions_dredged_ndmp_best"))

final <- model.avg(mlqid.3, subset = delta<2, fit=T) ## Average the models with delta AIC < value and ensure the models are actually fit

### Output the coefficients in a table
final.coefs <- final$coefficients["full",] ## Average parameter estimates over all models (value will be 0 in models wehre the parameter wasn't included)
vars <- names(final.coefs)
coefs <- as.vector(final.coefs)
final.coefs <- data.frame(Variables=vars, MeanCoefficients=coefs)
rownames(final.coefs) <- vars

## Write outputs
save(final, file=paste0(wd.out,"/mods_final_ndmp"))
write.csv(final.coefs, file=paste0(wd.out,"/mods_final_ndmp.csv"), row.names=F)

