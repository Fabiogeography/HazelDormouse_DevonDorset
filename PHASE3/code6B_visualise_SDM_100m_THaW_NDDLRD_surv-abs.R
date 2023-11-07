##############################################################
##### EVALUATE MODELS ########################
##### Models made with NDD+LRD presences and survey absences #####
##### Written by: Regan Early ################################
##### Written on: 28th March 2021 ##############################
##### Modified on: 30th March 2022 #########################
##############################################################

.libPaths("D:/SOFTWARE/R-4.1.2/library")
library(ggplot2) ## Needed to make interact_plot run and needs to be version for R 3.3.0 or later
library(interactions) ## interact_plot
# library(Hmisc) ## rcorr ## load after ggplot2 so it doesn't inadvertently load the wrong version of the package
# library(dplyr) ## bind_rows
# library(corrplot)
# library(raster)
# library(rgdal)
# library(spatialEco)
# library(sp)
# library(pROC)
# library(tidyverse) ## glimpse and useful functions
# library(tidyr) ## gather
# library(MuMIn) # dredge() stdize() model.sel()
# library(car) ## Variance Inflation Factor analysis
# library(snow) ## for running on multiple cores of a cluster. Used for dredge function (best when used on a server rather than a desktop).

wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/100m"

dat.m <- read.csv(paste0(wd.out, "/dat_std_22Feb2022_t_frag_m_sampled_nddlrd.csv"))

load(file=paste0(wd.out,"/mods_final_surv-abs_nddlrd")) ## This model averages all constituent models. Probably not meaningful for the interaction variables selected only a few times. 

##############################################################
##### Make graphs of response curves ######
##############################################################
final.coefs <- read.csv(file=paste0(wd.out,"/mods_final_surv-abs_nddlrd.csv"))
vars.f <- as.character(final.coefs$Variables)
vars.f <- subset(vars.f, !grepl("2)", vars.f)) ## remove quadratic terms
vars.f <- subset(vars.f, !grepl(":", vars.f)) ## remove interaction terms
vars.f <- subset(vars.f, !grepl("Intercept", vars.f)) ## remove Intercept
vars.names <- c("Ancient woodland", "NFI + LCM coniferous change", "Slope %",
                "Aspect WE", "Spring rainfall",
                "Spring sun", "Spring max temp", 
                "Winter min temp", "THaW trees+hedge",
                "Aspect NS", "NFI + LCM conif coverage", "NFI + LCM broadl coverage")


destdize <- read.csv(paste0(wd.out, "/dat_means_sds_23Mar2022_t_frag_sampled_prepped.csv"))

mns <- apply(dat.m[,vars.f], 2, mean)

##### Graph averaged model #####
## y.limits for graphs
ylims <- rep(1.2, length(vars.f))
names(ylims) <- vars.f

png(paste0(wd.out, "/final_response_curves_avg_t_nddlrd.png"), width=700, height=500)
par(mfrow=c(3,4))
par(mar=c(4,4.5,2,1)) # c(bottom, left, top, right)
# i <- vars.f[1]
for (i in vars.f) {
  newdat <- data.frame(matrix(data=mns, 
                              byrow=T, nrow=10000, ncol=length(vars.f), dimnames=(list(NULL,vars.f)))) ## Make a dataframe with mean (standardised) values of each variable
  newdat[,i] <- seq(min(dat.m[,i]), max(dat.m[,i]), length.out=10000) ## Vary the parameter of interest
  p <- predict(final, newdat, type="response", full=T, se.fit=T) ## re.form=NA, 
  
  ### Calculate the raw (destandardised) values of the explanatory variable
  ## Multiply by standard deviation (row 2 of destdize) and add mean (row 1 of destdize)
  x.std <- seq(min(newdat[,i]), max(newdat[,i]), length.out=6) ## Identify six values for plotting
  if(i %in% c("combo_broadl_m", "combo_conif_m")) { ## If variable is the NFI+LCM combined divide the axis display values by two to represent mean rather than sum. Sum was used in the model, but whether mean or sum were used will have no effect on model.
    x.destd <- round(((x.std * destdize[2,i]) + destdize[1,i])/2, 1) ## Multiply standardised value by sd of original data and add the mean
  } else {
    x.destd <- round((x.std * destdize[2,i]) + destdize[1,i], 1) ## Multiply standardised value by sd of original data and add the mean
  }  
  plot(newdat[,i], p$fit, type="n", xlab=vars.names[vars.f==i],
       ylab="Prob. occurrence", xaxt="n", yaxt="n",
       ylim=c(0,ylims[i]),#
       cex.lab=1.5) # p$fit
  axis(side=1, at=x.std, labels=x.destd, cex.axis=1.2) ## Add x axis on original scale
  axis(side=2, at=seq(0,1.5,0.3), cex.axis=1.2) ##, labels=c(0,format(ylims[i], scientific = FALSE)))
  l.ci <- p$fit - (1.96*p$se.fit)
  u.ci <- p$fit + (1.96*p$se.fit)
  polygon(c(newdat[,i], rev(newdat[,i])), c(c(l.ci, rev(u.ci))), col = 'grey80', border = NA)
  lines(newdat[,i], p$fit)
  
}

dev.off()

##### Interaction plots ##### 
## Can't enter a model average object into interact_plot so make a dummy model
## Make a dummy model to alter coefficients
test <- as.character(final.coefs$Variables, as.factor=F)
test <-subset(test, !grepl("Intercept", test))
test <- paste0(test, collapse=" + ")
dummy <- glm(as.formula(paste0("occ ~ ", test)), family=binomial, dat.m, na.action=na.fail) ## Make model with calibration data

## Update dummy model with coefficients from model average
for(n in names(dummy$coefficients)) {
  dummy$coefficients[n] <- final.coefs[final.coefs$Variables==n,"MeanCoefficients"]
}

## Identify the interactors
test <- as.character(final.coefs$Variables, as.factor=F)
interactors <- test[grepl(":", test)]
interactors <- strsplit(interactors, ":") ## eight interactions

## Remember can't trust SD because taken from original dummy model
## run manually, as interactors[[1]][1] (etc) doesn't seem to work

## Start with interactions in all of the best models:
# jpeg(paste0(wd.out,"/int_t_rainf_tasrng.jpg"), width=600, height=400)
# interact_plot(dummy, pred=rainfall_spr_5yrmn, modx=tasrng_win_5yrmn, interval=T, rug=T)
# dev.off()
# 
# jpeg(paste0(wd.out,"/int_t_scrub_tasmin.jpg"), width=600, height=400)
# interact_plot(dummy, pred=scrub, modx=tasmin_win_5yrmn, interval=T, rug=T)
# dev.off()
# 
# jpeg(paste0(wd.out,"/int_t_treehedge_scrub.jpg"), width=600, height=400)
# interact_plot(dummy, pred=treehedge, modx=scrub, interval=T, rug=T)
# dev.off()
# 
# jpeg(paste0(wd.out,"/int_t_treehedge_tasrng.jpg"), width=600, height=400)
# interact_plot(dummy, pred=treehedge, modx=tasrng_win_5yrmn, interval=T, rug=T)
# dev.off() 
# 
# ## in 31 of 50 models:
# jpeg(paste0(wd.out,"/int_t_scrub_tasmax.jpg"), width=600, height=400)
# interact_plot(dummy, pred=scrub, modx=tasmax_spr_5yrmn, interval=T, rug=T)
# dev.off()
# 
# ## Not plotting the other 4 as they appear so rarely. 