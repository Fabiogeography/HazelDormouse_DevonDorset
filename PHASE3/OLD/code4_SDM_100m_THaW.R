##############################################################
##### SDM ########################
##### Written by: Regan Early ################################
##### Written on: 22nd March 2021 ##############################
##### Modified on: 22 Feb 2022  #########################
##############################################################

.libPaths("C:\\SOFTWARE\\R-4.1.2\\library")
library(dplyr) ## bind_rows, sample_frac
library(MuMIn) # dredge() stdize() model.sel()
library(car) ## Variance Inflation Factor analysis
library(snow) ## for running on multiple cores of a cluster. Used for dredge function (best when used on a server rather than a desktop).

wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/100m"
wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
# wd.out <- "~/re259/UoE_U_Drive/DORMOUSE_NE/SDM/100m"
# wd.env <- "~/re259/UoE_U_Drive/DORMOUSE_NE/SDM/100m"

##### MAKE MULTIVARIATE MODELS #####
vars <- c('combo_broadl_m', 'combo_broadl_c', 'combo_conif_m', 'combo_conif_c', 'anc_wood', 'OS_Terrain_100_NS', 'OS_Terrain_100_WE',
          'OS_Terrain_100_slope_pct', 'treehedge', 'scrub', "tasrng_win_5yrmn", "tasmax_spr_5yrmn", "tasmin_win_5yrmn", "sun_spr_5yrmn", "rainfall_spr_5yrmn")

dat <- read.csv(paste0(wd.env,"/dat_21Feb2022_t_sampled_prepped.csv"), as.is=T) ## Currently using NDMP and NDD+LRD samples, as only 49 NDMP presences

##### Scale explanatory variables as large numerical scope of variables can cause estimation problems (https://r-sig-mixed-models.r-project.narkive.com/fXcaHABA/r-sig-me-cholmod-warning-with-glmer). #####
## First save the values used to stdize, so as to graph and map later
# means <- apply(na.omit(dat[,vars]), 2, mean)
# sds <- apply(na.omit(dat[,vars]), 2, sd)
# destdize <- rbind(means, sds)
# write.csv(destdize, paste0(wd.out, "/dat_means_sds_22Feb2022_t_sampled.csv"), row.names=F)
destdize <- read.csv(paste0(wd.out, "/dat_means_sds_22Feb2022_t_sampled.csv"))

## Use method stdize from MuMin, which subtracts the mean and then divides by the standard deviation
# dat.m <- na.omit(dat[,c("occ", vars, "source")]) ## NAs have to be removed for dredge function
# dat.m <- dat.m[dat.m$source!="pseudo-absence",] ## remove points with non-NDMP presences
# dat.m$source[dat.m$source %in% c("nonPTES-abs","PTES-abs")] <- "surv-abs"
# dat.m[, vars] <- stdize(dat.m[,vars])
# dat.m <- as.data.frame(cbind(dat.m, id=c(1:nrow(dat.m))))
# write.csv(dat.m, paste0(wd.out, "/dat_std_22Feb2022_t_m_sampled.csv"), row.names=F)
dat.m <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_m_sampled.csv"))

### Calibration and validation data to make the model
# dat.c <- sample_frac(cbind(dat.m), size=0.7)
# write.csv(dat.c, paste0(wd.out, "/dat_std_22Feb2022_t_c_sampled.csv"), row.names=F)
dat.c <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_c_sampled.csv"))

# dat.v <- dat.m[!(dat.m$id %in% dat.c$id),]
# write.csv(dat.v, paste0(wd.out, "/dat_std_22Feb2022_t_v_sampled.csv"), row.names=F)
dat.v <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_v_sampled.csv"))

##### Model with linear terms only #####
vars.m <- paste(vars, collapse=" + ")
ml <- glm(as.formula(paste0("occ ~ ", vars.m)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. 
## A model to be dredged must not have na.action as na.omit or na.exclude, or dredge will fail.

## Check for variance inflation - an effect of collinear explanatory variables
vif(ml, singular.ok = TRUE) ## broadl and treehedge a bit high (~6)

# ## Set up the cluster. Number of cores limited by computing resources - best on server
# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 10), type = clusterType))
# clusterExport(clust, "dat.c") ## Send the data to each node of the cluster
# clusterCall(clust, function() library(MuMIn)) ## Call the function to load the library on each node of the cluster

## Dredge model with linear terms. Previous version only permitted one of the two highly correlated terms. Actually runs in a few mins on desktop.
# mld.all <- pdredge(ml, beta="sd", evaluate=T, trace=2, cluster=clust)
mld.all <- dredge(ml, beta="sd", evaluate=T, trace=2)
# save(mld.all, file=paste0(wd.out, "/linear_dredged_all_t_sampled"))
# load(paste0(wd.out, "/linear_dredged_all_t_sampled"))
# write.csv(mld.all, file=paste0(wd.out,"/linear_dredged_all_t_sampled.csv"), row.names=F)
mld.all[mld.all$delta<=2,] ## 6 models fall into the best model subset. combo_broad_m, rainfall, scrub, sun, tasmax, tasmin, tasrange the consistent variables. NA, WE, and broadl mean flip in and out.
## Whi is this so different from the previous models with the outliying value for scrub???
# m.best <- get.models(md, 1)[[1]]
mld.all <- get.models(mld.all, subset=delta<2) ## mld is now an object containing the best models.

##### Model with linear and quadratic terms #####
### Get the linear terms that were selected in the top models and make into quadratic terms
vars.mld <- unique(unlist(lapply(mld.all, function(x) {rownames(summary(x)$coefficients)})))
vars.mld <- vars.mld[vars.mld!="(Intercept)"] 

## Slope, NS and WE aspect are each in just one of the 12 best models, combo_conif_c is in 2, so I'm not includng here
vars.mld <- vars.mld[vars.mld!="OS_Terrain_100_NS"] ## only in 1  model
vars.mld <- vars.mld[vars.mld!="OS_Terrain_100_WE"] ## only in 1  model
vars.mld <- vars.mld[vars.mld!="OS_Terrain_100_slope_pct"] ## only in 1  model
vars.mld <- vars.mld[vars.mld!="combo_conif_c"] ## only in 1 model

vars.mq <- paste0("I(", vars.mld, "^2)") ## quadratic terms
vars.mlq <- paste(c(vars, vars.mq), collapse=" + ") ## all linear and quadratic terms. was vars.mld

### Dredge the model with all linear terms retained in the best model subset, and all of their quadratic terms
mlq <- glm(as.formula(paste0("occ ~ ", vars.mlq)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. Fitted probabilities 0 or 1 occurred.
vif(mlq) ## combo_broadl_m (9.7) borderline. conif_m and anc_wood super high, but seems to be correlating strongly with quadratic forms

vars.mlq1 <- paste(c(vars, vars.mq[1:5]), collapse=" + ") ## all linear and first five quadratic terms
mlq1 <- glm(as.formula(paste0("occ ~ ", vars.mlq1)), family=binomial, dat.c, na.action=na.fail)
  
vars.mlq2 <- paste(c(vars, vars.mq[6:10]), collapse=" + ") ## all linear and first five quadratic terms
mlq2 <- glm(as.formula(paste0("occ ~ ", vars.mlq2)), family=binomial, dat.c, na.action=na.fail)

## Retain quadratic term only if linear term is added
msubset <- expression(dc("anc_wood", "I(anc_wood^2)") & ## linear effect is in all six models
                        dc("rainfall_spr_5yrmn", "I(rainfall_spr_5yrmn^2)") & ## linear effect is in all 12 models
                        dc("scrub", "I(scrub^2)") & ## linear effect is in all 12 models
                        dc("sun_spr_5yrmn", "I(sun_spr_5yrmn^2)") & ## linear effect is in all 12 models
                        dc("tasmax_spr_5yrmn", "I(tasmax_spr_5yrmn^2)") & ## linear effect is in all 12 models
                        dc("tasmin_win_5yrmn", "I(tasmin_win_5yrmn^2)") & ## linear effect is in all 12 models
                        dc("tasrng_win_5yrmn", "I(tasrng_win_5yrmn^2)") & ## linear effect is in 11 of 12 models
                        dc("treehedge", "I(treehedge^2)") & ## linear effect is in 10 of 12 models
                        dc("combo_conif_m", "I(combo_conif_m^2)") & ## linear effect is in 2 of 12 models
                        dc("combo_broadl_m", "I(combo_broadl_m^2)")) ## linear effect is in 3 of 12 models
                      
# mlqd <- pdredge(mlq, beta="sd", evaluate=T, trace=2, subset=msubset, cluster=clust)
mlqd1 <- dredge(mlq1, beta="sd", evaluate=T, trace=2, subset=msubset) ## super lengthy
mlqd2 <- dredge(mlq2, beta="sd", evaluate=T, trace=2, subset=msubset) ## super lengthy

save(mlqd, file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled"))
write.csv(mlqd, file=paste0(wd.out,"/linear_quadratic_dredged_t_sampled.csv"), row.names=F)

## Best subset of models with delta AICc <=2
mlqd[mlqd$delta<=2,] ## many models in the best model subset.
## sun, rainfall, treehedge, and combo broadleaf change quadratics rarely included. Aspect quadratic sometimes included. 
## combo conif quadratic always included
## linear effects of NS aspect and slope very rarely included

# mlqd <- get.models(mlqd, subset=delta<2) ## mlqd is now an object containing the three best models.

## Check for variance inflation - an effect of collinear explanatory variables
for(i in 1:length(mlqd)) {print(vif(mlqd[[i]]), singular.ok = TRUE)} ## ancient woodland has high VIF (~14), but seems to be with its own quadratic form, as only happens in models with quadratic.

# Inspect models (ask whether fixed effect relationships meaningful)
for(i in 1:length(mlqd)) {print(summary(mlqd[[i]]))}

##### Model with interactions between linear terms #####
## Create list of interactions to test
vars.mlqd2 <- unique(unlist(lapply(mlqd, function(x) {rownames(summary(x)$coefficients)})))
vars.mlqd2 <- vars.mlqd2[vars.mlqd2!="(Intercept)"]

# vars.mi <- combn(subset(vars.mlqd2, !grepl("2", vars.mlqd2)),2) ## selects every linear effect in the best models.
# vars.mi <- apply(vars.mi, MARGIN=2, FUN=function(x) {paste(x, collapse="*")})## From previous line
vars.mi <- combn(c("combo_broadl_c", "combo_conif_m", "OS_Terrain_100_WE", "rainfall_spr_5yrmn", "sun_spr_5yrmn", "treehedge"), 2) ## New method. Select only linear terms in every model.
vars.mi <- apply(vars.mi, MARGIN=2, FUN=function(x) {paste(x, collapse="*")})

## Group interactions into 2 subsets, otherwise there are too many fixed terms for model selection to work
vars.mi.1 <- paste(c(vars.mlqd2, vars.mi[1:8]), collapse= " + ")
vars.mi.2 <- paste(c(vars.mlqd2, vars.mi[9:15]), collapse= " + ")

## Old method
# vars.mi.1 <- paste0("(", paste(vars.mi[1:(length(vars.mi)/2)], collapse=") + ("),")")
# start <- 1+(length(vars.mi)/3); stop <- (length(vars.mi)/3)*2
# vars.mi.2 <- paste0("(", paste(vars.mi[start : stop], collapse=") + ("),")")
# start <- 1+(length(vars.mi)/2); stop <- length(vars.mi)
# vars.mi.3 <- paste0("(", paste(vars.mi[start : stop], collapse=") + ("),")")

mlqi.full1 <- glm(as.formula(paste0("occ ~ ", paste(c(vars.mlqd2, vars.mi.1), collapse=" + "))), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. 28 variables
mlqi.full2 <- glm(as.formula(paste0("occ ~ ", paste(c(vars.mlqd2, vars.mi.2), collapse=" + "))), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data

# Compare models using AICc and save selection table
mlqid.1 <- pdredge(mlqi.full1, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust) ##
save(mlqid.1, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled1"))
write.csv(mlqid.1, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled1.csv"), row.names=F)

mlqid.2 <- pdredge(mlqi.full2, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
save(mlqid.2, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled2"))
write.csv(mlqid.2, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled2.csv"), row.names=F)

## Best subset of models with delta AICc <=2. 
mlqid.1[mlqid.1$delta<=2,] ## 13 models in the best model subset. 12 variables retained. No interactions.
mlqid.2[mlqid.2$delta<=2,] ## Ditched ~ 15 variables. 14 models in best subset. 
 
## Get the terms in the best models
vars.mlqid1 <- names(model.sel(mlqid.1[mlqid.1$delta<=2,]))
vars.mlqid2 <- names(model.sel(mlqid.2[mlqid.2$delta<=2,]))

vars.mlqid <- c(vars.mlqid1, vars.mlqid2)
vars.mlqid <- unique(vars.mlqid[!(vars.mlqid %in% c("df", "logLik", "AICc", "delta", "weight", "(Intercept)"))]) ## 16 variables
vars.mlqid <- paste(vars.mlqid, collapse = " + ")

## Dredge the combined list of variables
mlqi.comb <- glm(as.formula(paste0("occ ~ ", vars.mlqid)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. Fitted probabilities 0 or 1 occurred.
vif(mlqid.comb) ## ok

mlqid.comb <- pdredge(mlqi.comb, beta="sd", evaluate=T, trace=2, subset=msubset, cluster=clust)

save(mlqid.comb, file=paste0(wd.out, "/linear_quadratic_interactions_dredged_nddlrd_best"))
write.csv(mlqid.comb, file=paste0(wd.out,"/linear_quadratic_interactions_dredged_nddlrd_best.csv"), row.names=F)

mlqid.comb[mlqid.comb$delta<=2,] ## 14 models in the best subset.


#################################################
###### Produce & save final model ###############
#################################################
load(file=paste0(wd.out, "/linear_quadratic_interactions_dredged_nddlrd_best"))

final <- model.avg(mlqid.comb, subset = delta<2, fit=T) ## Average the models with delta AIC < value and ensure the models are actually fit

### Output the coefficients in a table
final.coefs <- final$coefficients["full",] ## Average parameter estimates over all models (value will be 0 in models wehre the parameter wasn't included)
vars <- names(final.coefs)
coefs <- as.vector(final.coefs)
final.coefs <- data.frame(Variables=vars, MeanCoefficients=coefs)
rownames(final.coefs) <- vars

## Write outputs
save(final, file=paste0(wd.out,"/mods_final"))
write.csv(final.coefs, file=paste0(wd.out,"/mods_final.csv"), row.names=F)

