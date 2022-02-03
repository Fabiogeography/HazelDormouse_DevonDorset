##############################################################
##### EVALUATE MODELS ########################
##### Written by: Regan Early ################################
##### Written on: 28th March 2021 ##############################
##### Modified on: 11th November 2021  #########################
##############################################################

library(ggplot2, lib.loc="D:/SOFTWARE/R-4.1.1/library") ## Needed to make interact_plot run and needs to be version for R 3.3.0 or later
library(interactions, lib.loc="D:/SOFTWARE/R-4.1.1/library") ## interact_plot

wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/100m"

dat.m <- read.csv(paste0(wd.out, "/dat_std_11Nov2021_m.csv"))

load(file=paste0(wd.out,"/mods_final_noInt_combo"))

##############################################################
##### Make graphs of response curves ######
##############################################################
final.coefs <- read.csv(file=paste0(wd.out,"/mods_final_noInt_combo.csv"))
vars.f <- as.character(final.coefs$Variables)
vars.f <- subset(vars.f, !grepl("2)", vars.f)) ## remove quadratic terms
vars.f <- subset(vars.f, !grepl(":", vars.f)) ## remove interaction terms
vars.f <- subset(vars.f, !grepl("Intercept", vars.f)) ## remove Intercept
vars.names <- c("Change broadleaf coverage", "Broadleaf coverage", "Change coniferous coverage", "Coniferous coverage",
                "NS_aspect", "Slope", "Rainfall autumn", "Sun spring", "Min winter temperature",
                "WE_aspect")
                
destdize <- read.csv(paste0(wd.out, "/dat_means_sds_11Nov2021.csv"))

mns <- apply(dat.m[,vars.f], 2, mean)

##### Graph averaged model #####
## y.limits for graphs
ylims <- rep(1.2, length(vars.f))
names(ylims) <- vars.f

png(paste0(wd.out, "/final_response_curves_avg_noInt_combo.png"), width=600, height=500)
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
  x.destd <- round((x.std * destdize[2,i]) + destdize[1,i], 1) ## Multiply standardised value by sd of original data and add the mean
  
  plot(newdat[,i], p$fit, type="n", xlab=vars.names[vars.f==i],
       ylab="Prob. Settling", xaxt="n", yaxt="n",
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
interactors <- strsplit(interactors, ":")

## Remember can't trust SD because taken from original dummy model
## run manually, as interactors[[1]][1] (etc) doesn't seem to work
jpeg(paste0(wd.out,"/int_conif_all_m_tasmax_smr.jpg"), width=600, height=400)
interact_plot(dummy, pred=conif_all_m, modx=tasmax_smr_5yrmn, interval=T, rug=T)
dev.off()

jpeg(paste0(wd.out,"/int_conif_all_m_tasmin_win.jpg"), width=600, height=400)
interact_plot(dummy, pred=conif_all_m, modx=tasmin_win_5yrmn, interval=T, rug=T)
dev.off()

jpeg(paste0(wd.out,"/int_conif_all_m_sun_spr.jpg"), width=600, height=400)
interact_plot(dummy, pred=conif_all_m, modx=sun_spr_5yrmn, interval=T, rug=T)
dev.off()

jpeg(paste0(wd.out,"/int_tasmax_smr_sun_spr.jpg"), width=600, height=400)
interact_plot(dummy, pred=tasmax_smr_5yrmn, modx=sun_spr_5yrmn, interval=T, rug=T)
dev.off()

jpeg(paste0(wd.out,"/int_broadl_all_m_sun_spr.jpg.jpg"), width=600, height=400)
interact_plot(dummy, pred=broadl_all_m, modx=sun_spr_5yrmn, interval=T, rug=T)
dev.off()
