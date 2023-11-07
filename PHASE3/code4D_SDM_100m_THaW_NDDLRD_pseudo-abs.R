##############################################################
##### SDM ########################
##### Run SDM for THaW area with the pseudo-absences, to see comparison with survey-absences
##### Written by: Regan Early ################################
##### Written on: 22nd March 2021 ##############################
##### Modified on: 21st March 2022  #########################
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
# means <- apply(na.omit(dat[,vars.plus]), 2, mean)
# sds <- apply(na.omit(dat[,vars.plus]), 2, sd)
# destdize <- rbind(means, sds)
# write.csv(destdize, paste0(wd.out, "/dat_means_sds_23Mar2022_t_frag_sampled_prepped.csv"), row.names=F)
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
# write.csv(dat.m, paste0(wd.out, "/dat_std_22Feb2022_t_frag_m_sampled_nddlrd_pabs.csv"), row.names=F)
dat.m <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_frag_m_sampled_nddlrd_pabs.csv"))

### Calibration and validation data to make the model
# dat.c <- sample_frac(cbind(dat.m), size=0.7)
# write.csv(dat.c, paste0(wd.out, "/dat_std_22Feb2022_t_frag_c_sampled_nddlrd_pabs.csv"), row.names=F)
dat.c <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_frag_c_sampled_nddlrd_pabs.csv"))

# dat.v <- dat.m[!(dat.m$id %in% dat.c$id),]
# write.csv(dat.v, paste0(wd.out, "/dat_std_22Feb2022_t_frag_v_sampled_nddlrd_pabs.csv"), row.names=F)
dat.v <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_frag_v_sampled_nddlrd_pabs.csv"))

##### Model with linear terms only #####
vars.m <- paste(vars, collapse=" + ")
ml <- glm(as.formula(paste0("occ ~ ", vars.m)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. 
## A model to be dredged must not have na.action as na.omit or na.exclude, or dredge will fail.

## Check for variance inflation - an effect of collinear explanatory variables
vif(ml, singular.ok = TRUE) ## combo_broadl_m and treehedge are borderline high (6-7)

# ## Set up the cluster. Number of cores limited by computing resources - best on server
# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# clust <- try(makeCluster(getOption("cl.cores", 30), type = clusterType))
# clusterExport(clust, "dat.c") ## Send the data to each node of the cluster
# clusterCall(clust, function() library(MuMIn)) ## Call the function to load the library on each node of the cluster

## Dredge model with linear terms. Previous version only permitted one of the two highly correlated terms. Actually runs in a few mins on desktop.
# mld.all <- pdredge(ml, beta="sd", evaluate=T, trace=2, cluster=clust)
mld.all <- dredge(ml, beta="sd", evaluate=T, trace=2)
# save(mld.all, file=paste0(wd.out, "/linear_dredged_all_t_sampled_nddlrd_pabs"))
# load(paste0(wd.out, "/linear_dredged_all_t_sampled_nddlrd_pabs"))
# write.csv(mld.all, file=paste0(wd.out,"/linear_dredged_all_t_sampled_nddlrd_pabs.csv"), row.names=F)
mld.all[mld.all$delta<=2,] ## 18 models fall into the best model subset. Scrub (positive), tasmax (positive), treehedge(positive) are in all models. WE aspect (positive) is in 12 models and Slope (negative) is in 14 models. Sun (negative) is in 8 models, anc_wood (negative) is in 6 models. combo_broadl_m, combo_conif_m, NS, rainfall, are in 1-4 models. 
# m.best <- get.models(md, 1)[[1]]
mld.all <- get.models(mld.all, subset=delta<2) ## mld is now an object containing the best models.

##### Model with linear and quadratic terms #####
### Get the linear terms that were selected in the top models and make into quadratic terms
vars.mld <- unique(unlist(lapply(mld.all, function(x) {rownames(summary(x)$coefficients)})))
vars.mld <- vars.mld[vars.mld!="(Intercept)"] ## 13 variables

## Remove. Don't try quadratic effects of variables rarely included
vars.mld <- vars.mld[vars.mld!="combo_broadl_m"] ## only in 1 model
vars.mld <- vars.mld[vars.mld!="combo_conif_m"] ## only in 1 model
vars.mld <- vars.mld[vars.mld!="OS_Terrain_100_NS"] ## only in 4 model
vars.mld <- vars.mld[vars.mld!="rainfall_spr_5yrmn"] ## only in 4 model

vars.mq <- paste0("I(", vars.mld, "^2)") ## quadratic terms
vars.mlq <- paste(c(vars, vars.mq), collapse=" + ") ## all linear and quadratic terms. was vars.mld

### Dredge the model with all linear terms retained in the best model subset, and all of their quadratic terms
mlq <- glm(as.formula(paste0("occ ~ ", vars.mlq)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. glm.fit: fitted probabilities numerically 0 or 1 occurred 
vif(mlq) ## anc_wood correlated with own quadratic term. combo_broadl_m and treehedge somewhat correlated (8)

# vars.mlq1 <- paste(c(vars, vars.mq[1:5]), collapse=" + ") ## all linear and first five quadratic terms
# mlq1 <- glm(as.formula(paste0("occ ~ ", vars.mlq1)), family=binomial, dat.c, na.action=na.fail) ## glm.fit: fitted probabilities numerically 0 or 1 occurred 
#   
# vars.mlq2 <- paste(c(vars, vars.mq[6:11]), collapse=" + ") ## all linear and second six quadratic terms
# mlq2 <- glm(as.formula(paste0("occ ~ ", vars.mlq2)), family=binomial, dat.c, na.action=na.fail)

## Retain quadratic term only if linear term is added
msubset <- expression(dc("anc_wood", "I(anc_wood^2)") & ## 
                        dc("sun_spr_5yrmn", "I(sun_spr_5yrmn^2)") & ## 
                        dc("tasmax_spr_5yrmn", "I(tasmax_spr_5yrmn^2)") & ## 
                        dc("treehedge", "I(treehedge^2)") & ## 
                        dc("OS_Terrain_100_slope_pct", "I(OS_Terrain_100_slope_pct^2)") & 
                        dc("scrub", "I(scrub^2)") & 
                        dc("OS_Terrain_100_WE", "I(OS_Terrain_100_WE^2)"))
                      
# mlqd1 <- dredge(mlq1, beta="sd", evaluate=T, trace=2, subset=msubset) ## lengthy
# save(mlqd1, file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled1_nddlrd"))
# write.csv(mlqd1, file=paste0(wd.out,"/linear_quadratic_dredged_t_sampled1_nddlrd.csv"), row.names=T)
# mlqd1[mlqd1$delta<=2,] ## 14 models in the best model subset. Ancient woodland now positive but humped
# 
# mlqd2 <- dredge(mlq2, beta="sd", evaluate=T, trace=2, subset=msubset) ## lengthy
# save(mlqd2, file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled2_nddlrd"))
# write.csv(mlqd2, file=paste0(wd.out,"/linear_quadratic_dredged_t_sampled2_nddlrdwarn.csv"), row.names=F)
# mlqd2[mlqd2$delta<=2,] ## Some overlap in AIC between mlqd1 and 2. 

# ### Combine models 1 and 2
# load(file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled1"))
# load(file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled2"))
# mlqd <- merge(mlqd1, mlqd2)

mlqd <- dredge(mlq, beta="sd", evaluate=T, trace=2, subset=msubset) ## lengthy
# save(mlqd, file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled_nddlrd_pabs"))
# load(paste0(wd.out, "/linear_quadratic_dredged_t_sampled_nddlrd_pabs"))
# write.csv(mlqd, file=paste0(wd.out,"/linear_quadratic_dredged_t_sampled_nddlrd_pabs.csv"), row.names=F)
mlqd[mlqd$delta<=2,] ## 25 models. Slope always negative and asymptotic. Scrub, treehedge always positive. tasmax alwayus positibe and sometimes accelerating. WE normally positive, anc_wood quad nevewr selected, rainfall, sun, normally negative

# save(mlqd, file=paste0(wd.out, "/linear_quadratic_dredged_t_sampled12"))
# write.csv(mlqd, file=paste0(wd.out,"/linear_quadratic_dredged_t_sampled12.csv"), row.names=F)

mlqd <- get.models(mlqd, subset=delta<2) ## now an object containing the best models.

unique(unlist(lapply(mlqd, function(x) {rownames(summary(x)$coefficients)}))) ## 14 variables. 10 main effects
c(vars, vars.mq) ## started from 22 variables. 15 main effects

## Check for variance inflation - an effect of collinear explanatory variables
for(i in 1:length(mlqd)) {print(vif(mlqd[[i]]), singular.ok = TRUE)} ## combo_conif_m has high VIF (~11), but seems to be with its own quadratic form, as only happens in models with quadratic.

# Inspect models (ask whether fixed effect relationships meaningful)
for(i in 1:length(mlqd)) {print(summary(mlqd[[i]]))}

##### Model with interactions between linear terms #####
## Create list of interactions to test
vars.mlqd <- unique(unlist(lapply(mlqd, function(x) {rownames(summary(x)$coefficients)})))
vars.mlqd <- vars.mlqd[vars.mlqd!="(Intercept)"]

# vars.mi <- combn(subset(vars.mlqd, !grepl("2", vars.mlqd)),2) ## selects every linear effect in the best models. 55 combimnations

vars.mi <- combn(c("anc_wood", "OS_Terrain_100_slope_pct", "OS_Terrain_100_WE", "rainfall_spr_5yrmn", "scrub", "sun_spr_5yrmn", "tasmax_spr_5yrmn", "treehedge"), 2) ## Manual method. Select only linear terms in the majority of the best models - all are in at least 22 of 25 best models. 21 combinations
vars.mi <- apply(vars.mi, MARGIN=2, FUN=function(x) {paste(x, collapse="*")}) ## Turns into interactions. 28 combinations

## Group interactions into 2 subsets, otherwise there are too many fixed terms for model selection to work
vars.mi.1 <- paste(c(vars.mlqd, vars.mi[1:9]), collapse= " + ")
vars.mi.2 <- paste(c(vars.mlqd, vars.mi[10:19]), collapse= " + ")
vars.mi.3 <- paste(c(vars.mlqd, vars.mi[20:28]), collapse= " + ")

mlqi.full1 <- glm(as.formula(paste0("occ ~ ", vars.mi.1)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. glm.fit: fitted probabilities numerically 0 or 1 occurred 
mlqi.full2 <- glm(as.formula(paste0("occ ~ ", vars.mi.2)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data
mlqi.full3 <- glm(as.formula(paste0("occ ~ ", vars.mi.3)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. glm.fit: fitted probabilities numerically 0 or 1 occurred

# Compare models using AICc and save selection table
# mlqid.1 <- pdredge(mlqi.full1, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust) ##
mlqid.1 <- dredge(mlqi.full1, beta="sd", evaluate=T, trace=2, subset = msubset) ## not too slow?
save(mlqid.1, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs1")) ## 49 models in best model subset
write.csv(mlqid.1, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs1.csv"), row.names=F)

# mlqid.2 <- pdredge(mlqi.full2, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
mlqid.2 <- dredge(mlqi.full2, beta="sd", evaluate=T, trace=2, subset = msubset)
save(mlqid.2, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs2"))
write.csv(mlqid.2, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs2.csv"), row.names=F)

# mlqid.3 <- pdredge(mlqi.full3, beta="sd", evaluate=T, trace=2, subset = msubset, cluster=clust)
mlqid.3 <- dredge(mlqi.full3, beta="sd", evaluate=T, trace=2, subset = msubset)
save(mlqid.3, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs3"))
write.csv(mlqid.3, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs3.csv"), row.names=F)

## Best subset of models with delta AICc <=2. 
mlqid.1[mlqid.1$delta<=2,] ## 
mlqid.2[mlqid.2$delta<=2,] ## 
mlqid.3[mlqid.3$delta<=2,] ## 
mlqid <- merge(mlqid.1, mlqid.2)
mlqid <- merge(mlqid, mlqid.3)
mlqid[mlqid$delta<=2,] ## All from mlqid.2

## Get the terms in the best models
vars.mlqid1 <- names(model.sel(mlqid.1[mlqid.1$delta<=2,]))
vars.mlqid2 <- names(model.sel(mlqid.2[mlqid.2$delta<=2,])) ## models in best model subset are only from this set
vars.mlqid3 <- names(model.sel(mlqid.3[mlqid.3$delta<=2,]))

## Get the terms that weren't included in the best model subset
# extra <- vars.mlqid1[!(vars.mlqid1 %in% vars.mlqid2)]
# extra <- c(extra, vars.mlqid3[!(vars.mlqid3 %in% vars.mlqid2) & !(vars.mlqid3 %in% extra)]) ## 12 extra variables - too many
extra <- c("sun_spr_5yrmn:tasmax_spr_5yrmn", "tasmax_spr_5yrmn:treehedge", "anc_wood:OS_Terrain_100_WE", "anc_wood:scrub",	"anc_wood:treehedge") ## Manual check. Add interactions that were in the majority of the local best model subset

## Get the terms that **were** included in the best model subset
vars.mlqid2 <- unique(vars.mlqid2[!(vars.mlqid2 %in% c("df", "logLik", "AICc", "delta", "weight", "(Intercept)", "family"))]) ## 20 vars

## Get the maximal set of variables
vars.maximal <- paste(c(vars.mlqid2, extra), collapse = " + ") ## 25 variables

## Dredge the maximal list of variables
mlqi.comb <- glm(as.formula(paste0("occ ~ ", vars.maximal)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data. Fitted probabilities 0 or 1 occurred.
vif(mlqi.comb) ## anc_wood has huge colinearity with variables that interact with it

mlqid.comb <- dredge(mlqi.comb, beta="sd", evaluate=T, trace=2, subset=msubset) 
save(mlqid.comb, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs_best"))
write.csv(mlqid.comb, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs_best.csv"), row.names=F)

mlqid.comb[mlqid.comb$delta<=2,] ## .

##### Add fragmentation variables #####
frag.vars <- c("connect500_treehedge", "clumpy500_treehedge", "pacfrac500_treehedge") ## 

## Get the terms that **were** included in the best model subset
load(paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs_best"))
mlqid.comb <- get.models(mlqid.comb, subset=delta<2)
vars.mlqid.comb <- unique(unlist(lapply(mlqid.comb, function(x) {rownames(summary(x)$coefficients)}))) ## 26 vars
vars.mlqid.comb <- vars.mlqid.comb[vars.mlqid.comb!="(Intercept)"] 

## Get the maximal set of variables
vars.maximal <- paste(c(vars.mlqid.comb, frag.vars[2]), collapse = " + ") ## 28 variables. frag.vars[2] is clumpy500_treehedge

## Dredge the maximal list of variables
mlqi.frag <- glm(as.formula(paste0("occ ~ ", vars.maximal)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data.  Model with all three frag vars contained fitted probabilities 0 or 1 occurred but AIC is 271.54. Higher than best model without frag (284.25). Model with clumpy500_treehedge not as good as model with all (277.84), but still better than with none and no fitted probabilities 0 or 1
AIC(mlqi.frag) ##  Model with pacfrac500 only poorer
vif(mlqi.frag) ## 

mlqid.frag <- dredge(mlqi.frag, beta="sd", evaluate=T, trace=2, subset=msubset) 
save(mlqid.frag, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs_best_frag"))
write.csv(mlqid.frag, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs_best_frag.csv"), row.names=F)
mlqid.frag[mlqid.frag$delta<=0.2,] ## 31 models, 13 Glm fits 0 or 1

#################################################
###### Produce & save final model ###############
#################################################
load(file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled_nddlrd_pabs_best_frag"))

final <- model.avg(mlqid.frag, subset = delta<=2, fit=T) ## Average the models with delta AIC < value and ensure the models are actually fit

test <- get.models(mlqid.frag, subset=delta<=2)

### Output the coefficients in a table
final.coefs <- final$coefficients["full",] ## Average parameter estimates over all models (value will be 0 in models wehre the parameter wasn't included)
vars <- names(final.coefs)
coefs <- as.vector(final.coefs)
final.coefs <- data.frame(Variables=vars, MeanCoefficients=coefs)
rownames(final.coefs) <- vars

## Write outputs
save(final, file=paste0(wd.out,"/mods_final_t_sampled_nddlrd_pabs_frag"))
write.csv(final.coefs, file=paste0(wd.out,"/mods_final_t_sampled_nddlrd_pabs_frag.csv"), row.names=F)

