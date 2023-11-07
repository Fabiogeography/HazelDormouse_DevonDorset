##############################################################
##### SDM ########################
##### Models with NDDLRD presences, survey absences, THaW region #####
##### Written by: Regan Early ################################
##### Written on: 22nd March 2021 ##############################
##### Modified on: 30th March 2022  #########################
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
vars.plus <- c(vars, "connect500_treehedge","clumpy500_treehedge", "pacfrac500_treehedge")

dat <- read.csv(paste0(wd.env,"/dat_23Mar2022_t_frag_sampled_prepped.csv"), as.is=T) ## Currently using NDMP and NDD+LRD samples, as only 49 NDMP presences

##### Scale explanatory variables as large numerical scope of variables can cause estimation problems (https://r-sig-mixed-models.r-project.narkive.com/fXcaHABA/r-sig-me-cholmod-warning-with-glmer). #####
## First save the values used to stdize, so as to graph and map later
# means <- apply(na.omit(dat[,vars]), 2, mean)
# sds <- apply(na.omit(dat[,vars]), 2, sd)
# destdize <- rbind(means, sds)
# write.csv(destdize, paste0(wd.out, "/dat_means_sds_22Feb2022_t_sampled.csv"), row.names=F)
destdize <- read.csv(paste0(wd.out, "/dat_means_sds_23Mar2022_t_frag_sampled_prepped.csv"))

## Use method stdize from MuMin, which subtracts the mean and then divides by the standard deviation
# dat.m <- na.omit(dat[,c("occ", vars.plus, "source")]) ## NAs have to be removed for dredge function
# dat.m$source[dat.m$source %in% c("nonPTES-abs","PTES-abs")] <- "surv-abs"
# dat.m[, vars.plus] <- stdize(dat.m[,vars.plus])
# dat.m <- as.data.frame(cbind(dat.m, id=c(1:nrow(dat.m))))
# dat.m <- dat.m[dat.m$source!="NDMP",] ## remove 49 points with NDMP presences
# ## Trim the pseudo-absences to the same length as the survey absences
# n <- table(dat.m$source)["surv-abs"] ## number of survey absences
# set.seed(69)
# n <- c(sample(dat.m[dat.m$source=="pseudo-absence","id"], n), dat.m[dat.m$source=="LRD+NDD","id"]) ## cell ids to retain - random but repeatable sample of the pseudo-absences and all of the presences
# dat.m <- dat.m[dat.m$id %in% n,]
# write.csv(dat.m, paste0(wd.out, "/dat_std_22Feb2022_t_frag_m_sampled_nddlrd.csv"), row.names=F)
dat.m <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_frag_m_sampled_nddlrd.csv"))

### Calibration and validation data to make the model
# dat.c <- sample_frac(cbind(dat.m), size=0.7)
# write.csv(dat.c, paste0(wd.out, "/dat_std_22Feb2022_t_frag_c_sampled_nddlrd.csv"), row.names=F)
dat.c <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_frag_c_sampled_nddlrd.csv"))

# dat.v <- dat.m[!(dat.m$id %in% dat.c$id),]
# write.csv(dat.v, paste0(wd.out, "/dat_std_22Feb2022_t_frag_v_sampled_nddlrd.csv"), row.names=F)
dat.v <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_frag_v_sampled_nddlrd.csv"))

##### Model with linear terms only #####
vars.m <- paste(vars, collapse=" + ")
ml <- glm(as.formula(paste0("occ ~ ", vars.m)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. 
## A model to be dredged must not have na.action as na.omit or na.exclude, or dredge will fail.

## Check for variance inflation - an effect of collinear explanatory variables
vif(ml, singular.ok = TRUE) ## all fine

## Set up the cluster. Number of cores limited by computing resources - best on server
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 30), type = clusterType))
clusterExport(clust, "dat.c") ## Send the data to each node of the cluster
clusterCall(clust, function() library(MuMIn)) ## Call the function to load the library on each node of the cluster

## Dredge model with linear terms. Previous version only permitted one of the two highly correlated terms. Actually runs in a few mins on desktop.
# mld.all <- pdredge(ml, beta="sd", evaluate=T, trace=2, cluster=clust)
mld.all <- dredge(ml, beta="sd", evaluate=T, trace=2)
# save(mld.all, file=paste0(wd.out, "/linear_dredged_all_t_sampled_nddlrd"))
# load(paste0(wd.out, "/linear_dredged_all_t_sampled_nddlrd"))
# write.csv(mld.all, file=paste0(wd.out,"/linear_dredged_all_t_sampled_nddlrd.csv"), row.names=F)
mld.all[mld.all$delta<=2,] ## 23 models fall into the best model subset. anc_wood, slope, rainfall, sun, tasmax and treehedge are in all models. NS and WE are in the majority, tamin, tasrng, combo)cnifC and combo_conif_m re in several. Scrub is in 2 (positive).
# m.best <- get.models(md, 1)[[1]]
mld.all <- get.models(mld.all, subset=delta<2) ## mld is now an object containing the best models.

##### Model with linear and quadratic terms #####
### Get the linear terms that were selected in the top models and make into quadratic terms
vars.mld <- unique(unlist(lapply(mld.all, function(x) {rownames(summary(x)$coefficients)})))
vars.mld <- vars.mld[vars.mld!="(Intercept)"] ## 13 variables

## RemoveDon't try quadratic effects of variables rarely included
vars.mld <- vars.mld[vars.mld!="scrub"] ## only in 2  model
vars.mld <- vars.mld[vars.mld!="combo_conif_m"] ## only in 2 model

vars.mq <- paste0("I(", vars.mld, "^2)") ## quadratic terms
vars.mlq <- paste(c(vars, vars.mq), collapse=" + ") ## all linear and quadratic terms. was vars.mld

### Dredge the model with all linear terms retained in the best model subset, and all of their quadratic terms
mlq <- glm(as.formula(paste0("occ ~ ", vars.mlq)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data
vif(mlq) ## Some vry concerningly high values, but see,m to be only with own quadrativ forms: combo_conif_c + I(combo_conif_c^2), anc_wood + I(anc_wood^2), OS_Terrain_100_NS + I(OS_Terrain_100_NS^2), OS_Terrain_100_WE + I(OS_Terrain_100_WE^2), 

vars.mlq1 <- paste(c(vars, vars.mq[1:5]), collapse=" + ") ## all linear and first five quadratic terms
mlq1 <- glm(as.formula(paste0("occ ~ ", vars.mlq1)), family=binomial, dat.c, na.action=na.fail) ## glm.fit: fitted probabilities numerically 0 or 1 occurred 
  
vars.mlq2 <- paste(c(vars, vars.mq[6:11]), collapse=" + ") ## all linear and second six quadratic terms
mlq2 <- glm(as.formula(paste0("occ ~ ", vars.mlq2)), family=binomial, dat.c, na.action=na.fail)

## Retain quadratic term only if linear term is added
msubset <- expression(dc("anc_wood", "I(anc_wood^2)") & ## 
                        dc("rainfall_spr_5yrmn", "I(rainfall_spr_5yrmn^2)") & ## 
                        dc("sun_spr_5yrmn", "I(sun_spr_5yrmn^2)") & ## 
                        dc("tasmax_spr_5yrmn", "I(tasmax_spr_5yrmn^2)") & ## 
                        dc("tasmin_win_5yrmn", "I(tasmin_win_5yrmn^2)") & ## 
                        dc("tasrng_win_5yrmn", "I(tasrng_win_5yrmn^2)") & ##
                        dc("treehedge", "I(treehedge^2)") & ## 
                        dc("combo_conif_c", "I(combo_conif_c^2)") & 
                        dc("OS_Terrain_100_NS", "I(OS_Terrain_100_NS^2)") & 
                        dc("OS_Terrain_100_slope_pct", "I(OS_Terrain_100_slope_pct^2)") & 
                        dc("OS_Terrain_100_WE", "I(OS_Terrain_100_WE^2)"))
                      
mlqd1 <- pdredge(mlq1, beta="sd", evaluate=T, trace=2, subset=msubset, cluster=clust) ## lengthy
save(mlqd1, file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled1_nddlrd"))
write.csv(mlqd1, file=paste0(wd.out,"/linear_quadratic_dredged_t_sampled1_nddlrd.csv"), row.names=T)
mlqd1[mlqd1$delta<=2,] ## 14 models in the best model subset. Ancient woodland now positive but humped

mlqd2 <- pdredge(mlq2, beta="sd", evaluate=T, trace=2, subset=msubset, cluster=clust) ## lengthy
save(mlqd2, file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled2_nddlrd"))
write.csv(mlqd2, file=paste0(wd.out,"/linear_quadratic_dredged_t_sampled2_nddlrd.csv"), row.names=F)
mlqd2[mlqd2$delta<=2,] ## 27 models in best subset. 

### Combine models 1 and 2
load(file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled1_nddlrd"))
load(file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled2_nddlrd"))
mlqd <- merge(mlqd1, mlqd2) ## 14 models in the best model subset. All from mlqd1

save(mlqd, file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled12_nddlrd"))
write.csv(mlqd, file=paste0(wd.out,"/linear_quadratic_dredged_t_sampled12_nddlrd.csv"), row.names=F)
load(file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled12_nddlrd"))
mlqd[mlqd$delta<=2,]

mlqd <- get.models(mlqd, subset=delta<2) ## now an object containing the 14 best models.

unique(unlist(lapply(mlqd, function(x) {rownames(summary(x)$coefficients)}))) ## 17 variables. 14 main effects. 
c(vars, vars.mq) ## started from 26 variables. 15 main effects

## Check for variance inflation - an effect of collinear explanatory variables
for(i in 1:length(mlqd)) {print(vif(mlqd[[i]]), singular.ok = TRUE)} ## anc_wood and slope has high VIFs, but seems to be with their own quadratic form.

# Inspect models (ask whether fixed effect relationships meaningful)
for(i in 1:length(mlqd)) {print(summary(mlqd[[i]]))}

##### Model with interactions between linear terms #####
## Create list of interactions to test
vars.mlqd <- unique(unlist(lapply(mlqd, function(x) {rownames(summary(x)$coefficients)})))
vars.mlqd <- vars.mlqd[vars.mlqd!="(Intercept)"]

# vars.mi <- combn(subset(vars.mlqd, !grepl("2", vars.mlqd)),2) ## selects every linear effect in the best models. 91 combinations. Too many. 
vars.mi <- combn(c("anc_wood", "OS_Terrain_100_slope_pct", "OS_Terrain_100_WE", "rainfall_spr_5yrmn",	"sun_spr_5yrmn",	"tasmax_spr_5yrmn",	"tasmin_win_5yrmn",	"treehedge"), 2) ## Manual method. Select only linear terms in 10 or more of the 14 best models. 28 combinations
vars.mi <- apply(vars.mi, MARGIN=2, FUN=function(x) {paste(x, collapse="*")}) ## Turns into interactions. 28 combinations

## Group interactions into 2 subsets, otherwise there are too many fixed terms for model selection to work
vars.mi.1 <- paste(c(vars.mlqd, vars.mi[1:9]), collapse= " + ")
vars.mi.2 <- paste(c(vars.mlqd, vars.mi[10:18]), collapse= " + ")
vars.mi.3 <- paste(c(vars.mlqd, vars.mi[19:28]), collapse= " + ")

mlqi.full1 <- glm(as.formula(paste0("occ ~ ", vars.mi.1)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. glm.fit: fitted probabilities numerically 0 or 1 occurred 
mlqi.full2 <- glm(as.formula(paste0("occ ~ ", vars.mi.2)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data
mlqi.full3 <- glm(as.formula(paste0("occ ~ ", vars.mi.3)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. glm.fit: fitted probabilities numerically 0 or 1 occurred 

# Compare models using AICc and save selection table
mlqid.1 <- pdredge(mlqi.full1, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust) ##
save(mlqid.1, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled1_nddlrd"))
write.csv(mlqid.1, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled1_nddlrd.csv"), row.names=F)

mlqid.2 <- pdredge(mlqi.full2, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
save(mlqid.2, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled2_nddlrd"))
write.csv(mlqid.2, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled2_nddlrd.csv"), row.names=F)

mlqid.3 <- pdredge(mlqi.full3, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
save(mlqid.3, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled3_nddlrd"))
write.csv(mlqid.3, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled3_nddlrd.csv"), row.names=F)

## Best subset of models with delta AICc <=2. 
mlqid.1[mlqid.1$delta<=2,] ## 90 models in the best model subset. Many with GLM fits 0 or 1. Interactions inconsistently retained. Best AIC is 290.8.
mlqid.2[mlqid.2$delta<=2,] ## 25 models in the best model subset. 21 with GLM fits 0 or 1. rnf_spr_5yr:trh is only interaction consistently retained (all but one best model). Best AIC is 296.8.
mlqid.3[mlqid.3$delta<=2,] ## ...
mlqid <- merge(mlqid.1, mlqid.2)
mlqid <- merge(mlqid, mlqid.3)

test <- mlqid[mlqid$delta<=2,]
save(test, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled12_nddlrd_AIC2"))
write.csv(test, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled12_nddlrd_AIC2.csv"), row.names=F)

        
# ## Get the terms in the best models
# vars.mlqid1 <- names(model.sel(mlqid.1[mlqid.1$delta<=2,]))
# vars.mlqid2 <- names(model.sel(mlqid.2[mlqid.2$delta<=2,]))
# vars.mlqid3 <- names(model.sel(mlqid.3[mlqid.3$delta<=2,]))
# 
# vars.mlqid <- c(vars.mlqid1, vars.mlqid2)
# vars.mlqid <- unique(vars.mlqid[!(vars.mlqid %in% c("df", "logLik", "AICc", "delta", "weight", "(Intercept)"))]) ## 16 variables
# vars.mlqid <- paste(vars.mlqid, collapse = " + ")
# 
# ## Dredge the combined list of variables
# mlqi.comb <- glm(as.formula(paste0("occ ~ ", vars.mlqid)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. Fitted probabilities 0 or 1 occurred.
# vif(mlqid.comb) ## ok
# 
# mlqid.comb <- pdredge(mlqi.comb, beta="sd", evaluate=T, trace=2, subset=msubset, cluster=clust)
# 
# save(mlqid.comb, file=paste0(wd.out, "/linear_quadratic_interactions_dredged_nddlrd_best"))
# write.csv(mlqid.comb, file=paste0(wd.out,"/linear_quadratic_interactions_dredged_nddlrd_best.csv"), row.names=F)
# 
# mlqid.comb[mlqid.comb$delta<=2,] ## 14 models in the best subset.

##### Add fragmentation variables #####
frag.vars <- c("connect500_treehedge", "clumpy500_treehedge", "pacfrac500_treehedge") ## 

## Get the terms that **were** included in the best model subset
load(paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled12_nddlrd_AIC2"))
mlqid.comb <- test
mlqid.comb <- get.models(mlqid.comb, subset=delta<2) ## Many algorithms didn't converge or had 0 or 1 fitted probabilities
vars.mlqid.comb <- unique(unlist(lapply(mlqid.comb, function(x) {rownames(summary(x)$coefficients)}))) ## 27 vars
vars.mlqid.comb <- vars.mlqid.comb[vars.mlqid.comb!="(Intercept)"] 

## Get the maximal set of variables
vars.maximal <- paste(c(vars.mlqid.comb, frag.vars), collapse = " + ") ## 25 variables. frag.vars[2] is clumpy500_treehedge

## Dredge the maximal list of variables
mlqi.frag <- glm(as.formula(paste0("occ ~ ", vars.maximal)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data.  AIC of best model without frag is 290.8
AIC(mlqi.frag) ##  Model with all three frag vars AIC is 302.4. Higher than best model without frag Models with individual frag vars worse as well. 
vif(mlqi.frag) ## 

## Don't run as fragmentation doesn't improve models
# mlqid.frag <- dredge(mlqi.frag, beta="sd", evaluate=T, trace=2, subset=msubset) 
# save(mlqid.frag, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled_best_frag"))
# write.csv(mlqid.frag, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled_best_frag.csv"), row.names=F)

#################################################
###### Produce & save final model ###############
#################################################
load(file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled12_nddlrd_AIC2"))
mlqid.comb <- test
final <- model.avg(mlqid.comb, subset = delta<2, fit=T) ## Average the models with delta AIC < value and ensure the models are actually fit

### Output the coefficients in a table
final.coefs <- final$coefficients["full",] ## Average parameter estimates over all models (value will be 0 in models wehre the parameter wasn't included)
vars <- names(final.coefs)
coefs <- as.vector(final.coefs)
final.coefs <- data.frame(Variables=vars, MeanCoefficients=coefs)
rownames(final.coefs) <- vars

## Write outputs
save(final, file=paste0(wd.out,"/mods_final_surv-abs_nddlrd"))
write.csv(final.coefs, file=paste0(wd.out,"/mods_final_surv-abs_nddlrd.csv"), row.names=F)

