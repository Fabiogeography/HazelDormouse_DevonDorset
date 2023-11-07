
wd.out <- "~/re259/UoE_U_Drive/DORMOUSE_NE/SDM/100m"

d1 <- read.csv(paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled1.csv"))
d2 <- read.csv(paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled2.csv"))
# d <- read.csv(paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled123.csv"))

dd <- d2[!is.na(d2$combo_broadl_m),] ## rows that have LCM + NFI broadleaf. First model has a delta AIC of 1.77229.but also contains treehedge
ddd <- dd[is.na(dd$treehedge),] ## rows that don't have treehedge. 
ddd[1,]  

## First model in entire combination has a delta AIC of 317.3468 is row 793 in file 123.
## First model in first subset has an overall AIC of 321.3906 (a delta AIC of 2.047) line 37 in file 1 (subset 1) 
## First model in second subset has an overall AIC of 317.3468 (a delta AIC of 4.255955) line 664 in file 2 (subset 2). This matches the first model in the entire combination so use this in mlqid.2

load(paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled2")) ## load first subset (5.2Gb)

m <- get.models(mlqid.2, 664)[[1]] ## Has different coefficients to the model in the csv. The coefficients match the model I made with teh correct dat.c below though, so I'll use it. 
summary(m)
d2[664,]


### Output the coefficients in a table
final.coefs <- m$coefficients ## Average parameter estimates over all models (value will be 0 in models wehre the parameter wasn't included)
vars <- names(final.coefs)
coefs <- as.vector(final.coefs)
final.coefs <- data.frame(Variables=vars, Coefficients=coefs)
rownames(final.coefs) <- vars

## Write outputs
save(final, file=paste0(wd.out,"/mods_final_surv-abs_NFI+LCM"))
write.csv(final.coefs, file=paste0(wd.out,"/mods_final_surv-abs_NFI+LCM.csv"), row.names=F)


##### Construct model from variable table #####
b <- ddd[1,]

### Get the variables
vars <- names(b[,!is.na(b)])
vars <- vars[!(vars %in% c("X.Intercept.","df","logLik","AICc","delta","weight"))]
vars <- gsub("I.", "I(", vars)
vars <- gsub(".2.", "^2)", vars)
vars <- gsub("\\.", ":", vars)
vars <- paste(vars, collapse=" + ")

m2 <- glm(as.formula(paste0("occ ~ ", vars)), family=binomial, dat.c, na.action=na.fail)



##### Add fragmentation variables #####
### Only out of interest - can't include in best NFI+LCM model as frag uses THaW
frag.vars <- c("connect500_treehedge", "clumpy500_treehedge", "pacfrac500_treehedge") ## 

## Get the terms that **were** included in the best model subset
# load(paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled_best"))
# mlqid.comb <- get.models(mlqid.comb, subset=delta<2)
vars.m.comb <- rownames(summary(m)$coefficients)
vars.m.comb <- vars.m.comb[vars.m.comb!="(Intercept)"] 

## Get the maximal set of variables
vars.maximal <- paste(c(vars.m.comb, frag.vars), collapse = " + ") ## 16 variables. frag.vars[2] is clumpy500_treehedge

## Dredge the maximal list of variables
m.frag <- glm(as.formula(paste0("occ ~ ", vars.maximal)), family=binomial, dat.c, na.action=na.fail) ## Make model with calibration data.   
AIC(m.frag) ##  AIC with all three frag vars (314.3), and with clumpy only (316.6) is better than model without (319.9)
vif(m.frag) ## 

## Don't run as fragmentation doesn't improve models
# mlqid.frag <- dredge(mlqi.frag, beta="sd", evaluate=T, trace=2, subset=msubset) 
# save(mlqid.frag, file=paste0(wd.out, "/linear_quadratic_interactions_all_t_sampled_best_frag"))
# write.csv(mlqid.frag, file=paste0(wd.out,"/linear_quadratic_interactions_all_t_sampled_best_frag.csv"), row.names=F)


dat.c <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_frag_c_sampled.csv"))
