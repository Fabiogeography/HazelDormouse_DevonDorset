##############################################################
##### Explore data before making SDM ########################
##### In this variant of the code the presences used (occ100mt_sampled) have been sampled to match the years of the survey absences #####
##### Written by: Regan Early ################################
##### Written on: 11th Jan 2021 ##############################
##### Modified on: 23 March 2022  #########################
##############################################################

.libPaths("C:\\SOFTWARE\\R-4.1.2\\library")
library(raster)
library(rgdal)
library(pROC)
library(tidyverse) ## glimpse and useful functions
library(corrplot)
library(Hmisc) ## rcorr

wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
wd.distn <- "E:/NON_PROJECT/DORMOUSE_NE/DISTRIBUTION"
wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/EXPLORE"
wd.clim <- "E:/NON_PROJECT/DORMOUSE_WALES/GIS/100m" # "E:/GIS_DATA/CLIMATE/UK/HadUK-Grid_1km"
OSGB.proj <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'

## Devon Dorset area
# dd <- readOGR(dsn="E:/NON_PROJECT/DORMOUSE_NE/GIS", layer="devondorset_buffer5km", stringsAsFactors = F)
## THaW area
t100 <- raster("E:/GIS_DATA/LandCover/UK/Devon_THaW/thaw_region.tif")

##### OBTAIN DATA #####
##### Dormouse records #####
occ_t <- readOGR(dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="occ100mt_sampled") ## load, so always working with the same dataset

##### CLIMATE DATA #####
### This code is for the new survey absences. Climate values were assigned to the cells wtih presences and pseudo-abseces in PHASE 2.
### Create mask of occupied grid-cells for each year, multiply each year raster by the mask to nullify cells with no records that year, and add all rasters together
clim.vars <- c("rainfall_spr", "rainfall_smr", "rainfall_aut", "sun_spr", "sun_smr", "sun_aut", "tasmax_spr", "tasmax_smr", "tasmax_aut", "tasmin_win", "tasrng_win")

## Restrict to extent of occupancy data, to save time
e <- extent(t100) ## geographic extent of the study region

r <- raster(paste0(wd.clim, "/rainfall_spr_5yrmn_1990.tif")) ## These were created for the dormouse Wales project
r <- crop(r, e)
r <- (r/r) - 1 ## make all values 0

# a <- occ[occ$occ== 0, ]
for(i in min(occ_t$year):max(occ_t$year)) { ## Note that although climate data are available from 2020 now, there are few dormouse records for 2020. 
  occ.yr <- occ_t[occ_t$year==i,]
  
  r.yr <- r ## make the mask grid
  r.yr[occ.yr] <- 1 ## set grid-cells in mask grid that have a presence or absence in the given year to be 1

  s <- stack(paste0(wd.clim, "/",clim.vars,"_5yrmn_",i,".tif"))
  s <- crop(s, e)
  s.yr <- s * r.yr ## Only grid-cells that contain presence or absence data in the given year retain a value. Others are given 0. 
  names(s.yr) <- substr(names(s), 1, nchar(names(s))-5)
  assign(paste0("clim.",i), s.yr)

}

### Adds the climate at each grid-cell where species is present or absent. problem could arise if value of any variable is precisely 0, as this is used as the numm value. Likelihood of this is low in this study region, so impact low. 
clim <- get(paste0("clim.", min(occ_t$year)))

for (y in min(occ_t$year)+1:max(occ_t$year)) { 
  temp <- get(paste0("clim.",y))
  clim <- clim + temp
  print(y)
}

names(clim) <- paste0(substr(names(s), 1, nchar(names(s))-5), "t")

writeRaster(clim, "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m/clim100m_t_surveyabsences_sampled")
clim <- stack("E:/NON_PROJECT/DORMOUSE_NE/GIS/100m/clim100m_t_surveyabsences_sampled")

## Load the env data made in PHASE 2 and unite with climate data for presences and survey absences
env.r <- stack(paste0(wd.env, "/env_100m7_dd_16thDec2021")) 
nms <- names(env.r)
env.r <- crop(env.r, extent(t100))
names(env.r) <- nms

env.r <- stack(env.r, clim)
nms <- names(env.r)
env.r <- env.r*t100
names(env.r) <- nms

##### Add THaW data #####
scrub <- raster("E:\\GIS_DATA\\LandCover\\UK\\Devon_THaW\\scrub_100m.tif")
names(scrub) <- "scrub"
# hedge <- raster("E:\\GIS_DATA\\LandCover\\UK\\Devon_THaW\\hedge_100m.tif") ## Failed to read directory at offset 314272346
# names(hedge) <- "hedge"
treehedge <- raster("E:\\GIS_DATA\\LandCover\\UK\\Devon_THaW\\trees_hedge_100m.tif")
names(treehedge) <- "treehedge"

env.r <- stack(env.r, scrub, treehedge) #, hedge

##### Add THaW fragmentation data #####
frag.layers <- c("clumpy100_treehedge","connect100_treehedge","ta100_treehedge","pacfrac100_treehedge","connect500_treehedge","clumpy500_treehedge","ta500_treehedge","pacfrac500_treehedge")
frag <- stack(paste0(wd.env, "/",frag.layers, ".tif")) ## seems to be a problem with ta500_treehedge
env.r <- stack(env.r, frag) 

## Save raster stack
writeRaster(env.r, paste0(wd.env, "/env_100m7_t_sampled_frag")) 
env.r <- stack(paste0(wd.env, "/env_100m7_t_sampled_frag")) 

##### Extract the environmental data at the presences, original pseudo-absences and survey absences #####
### Created in PHASE 2
occ_dd <- readOGR(dsn=wd.distn, layer="occ100m1990_dd") ## Note that does not have m100 values for pseudo-absences
occ_dd <- occ_dd[!(is.na(occ_dd$m100) & occ_dd$occ==1),] ## remove the once presence without a m100 value
occ_dd$source[occ_dd$occ==0] <- "pseudo-absence"
occ_dd$m100 <- 1 ## add dummy m100 value so data frame can cope with na.omit later on

## Add extra columns to occ_dd
x <- names(occ_t)[!(names(occ_t) %in% names(occ_dd))] ## Add columns in p to a
occ_dd@data[,x] <- NA 

## Merge dd and t data
abs_dd <- occ_dd[occ_dd$source == "pseudo-absence",]
abs_dd$easting <- coordinates(abs_dd)[,1]
abs_dd$northing <- coordinates(abs_dd)[,2]

occ <- rbind(occ_t, abs_dd)
occ$sei[is.na(occ$sei)] <- -9999 ## add dummy value rather than NA so data frame can cope with na.omit later on
# writeOGR(occ, dsn=paste0(wd.distn,"\\SURVEY_ABSENCES"), layer="occ100mt_sampled_dd", driver="ESRI Shapefile") ## load, so always working with the same dataset

env <- cbind(occ, raster::extract(env.r, occ)) ## spatial points data frame with env data. Note 'extract' is also a function in tidyr
dat <- as.data.frame(env) ## Make dataframe for easier modelling
glimpse(dat) 
summary(dat) ## 330 records fall outside THaW range (NA)
## There are 381 NA values for clumpy100 and 572 for pacfrac100. Investigating in ArcMap suggests that strange values were assigned to some of the 100m cells, which were outside teh values for the indices so I converted to NA. Just use 500m variables for the analysis. 

## Testing the 690 records with NA in env values are indeed outside the ThaW region
test <- dat[is.na(dat$lcm_broadl_m),]
plot(test$coords.x1, test$coords.x2)
plot(env.r[[1]], add=T)

test2 <- dat[is.na(dat$clumpy100_treehedge),]
plot(env.r[["clumpy100_treehedge"]])
points(test2$easting, test2$northing)

## Remove the 330 records. Note some NAs will still remain. 
dat <- dat[!is.na(dat$lcm_broadl_m),] 

## Remove one value with a huge outlying value for scrub, of 13%, standardised value of 17.309, which leads to fitted probabilities 0 or 1 in models with quadratic form of scrub.
which(dat$scrub>0.13) ## row 241
dat[240:242,] ## check outlier row
dat <- dat[-241,] 

# write.csv(dat, paste0(wd.env,"/dat_23Mar2022_t_frag.csv"), row.names=F)
dat <- read.csv(paste0(wd.env,"/dat_23Mar2022_t_frag.csv"))

####### visualise ################
## Add the climate values for the ThaW survey absences to the plot data. I didn't do this in the dat dataframe, as some pseudo-absences might overlap with the survey absences, and I didn't want to lose the data for the pseudo-absenes
clim.vars <- c('tasmin_win_5yrmn', 'tasrng_win_5yrmn', 'tasmax_spr_5yrmn', 'sun_spr_5yrmn',
               'rainfall_spr_5yrmn', 'tasmax_smr_5yrmn', 'sun_smr_5yrmn', 'rainfall_smr_5yrmn', 'tasmax_aut_5yrmn', 'sun_aut_5yrmn', 'rainfall_aut_5yrmn')
dat.plot <- dat

for (n in clim.vars) {
  # n <- clim.vars[1]
  nt <- paste0(n,"t")
  dat.plot[dat.plot$source %in% c("NDMP", "LRD+NDD", "PTES", "Non-PTES"), n] <- dat.plot[dat.plot$source %in% c("NDMP", "LRD+NDD", "PTES", "Non-PTES"), nt]
}

# ## Temporarily remove all survey absences from 2020 and 2021 to remove 0s from climate graphs, which prevent us seeing full range of data. Should include climate data for them
# dat.plot <- dat.plot[dat.plot$year<2020,]

## Plot division of all numeric variables between the different types of presence and absence points
vars <- c('combo_broadl_m', 'combo_broadl_c', 'combo_conif_m', 'combo_conif_c', 'anc_wood', 'OS_Terrain_100_NS', 'OS_Terrain_100_WE',
         'OS_Terrain_100_slope_pct', 'treehedge', 'scrub', 'tasmin_win_5yrmn', 'tasrng_win_5yrmn', 'tasmax_spr_5yrmn', 'sun_spr_5yrmn',
         'rainfall_spr_5yrmn', 'tasmax_smr_5yrmn', 'sun_smr_5yrmn', 'rainfall_smr_5yrmn', 'tasmax_aut_5yrmn', 'sun_aut_5yrmn', 'rainfall_aut_5yrmn',
         'clumpy100_treehedge','connect100_treehedge','ta100_treehedge','pacfrac100_treehedge','connect500_treehedge','clumpy500_treehedge','ta500_treehedge','pacfrac500_treehedge')

dat.plot <- cbind(dat.plot$source, dat.plot$year, dat.plot[,vars])
colnames(dat.plot)[1] <- "source"
colnames(dat.plot)[2] <- "year"
dat.plot$source[dat.plot$source=="PTES"] <- "PTES-abs"
dat.plot$source[dat.plot$source=="Non-PTES"] <- "nonPTES-abs"

## Adjust values of coverage for NFI+LCM forest, due to error in creating these layers
dat.plot$combo_broadl_m <- dat.plot$combo_broadl_m/2
dat.plot$combo_conif_m <- dat.plot$combo_conif_m/2

## Set orders of panels and boxplots
dat.plot$source <- factor(dat.plot$source, levels = c("NDMP", "LRD+NDD", "pseudo-absence", "PTES-abs", "nonPTES-abs"))

## Save for modelling
dat.mod <- as.data.frame(cbind(dat[,c("IDorig", "easting", "northing", "m100", "occ", "sei")], dat.plot))
dat.mod <- dat.mod[!is.na(dat.mod$treehedge),] ## Annoying two points!
# write.csv(dat.mod, paste0(wd.env,"/dat_23Mar2022_t_frag_sampled_prepped.csv"), row.names=F)

jpeg(paste0(wd.out,"/var_boxplots_t_sampled.jpg"), width=600, height=800)
dat.plot %>%
  gather(-source, key = "var", value = "value") %>%
  ggplot(aes(x=as.factor(source), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  # facet_grid(vars(vars), scales = "free") + #factor(size, levels=c('50%','100%','150%','200%')
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) 
dev.off()

par(mfrow=c(3,2))
plot(dat.plot$year, dat.plot$tasmin_win)
plot(dat.plot$year, dat.plot$tasmax_spr)
plot(dat.plot$year, dat.plot$sun_spr)
plot(dat.plot$year, dat.plot$sun_smr)
plot(dat.plot$year, dat.plot$rainfall_smr)

### Plot fragmentation values
dat.frag <- dat[,c("source", "treehedge","clumpy100_treehedge","connect100_treehedge","ta100_treehedge","pacfrac100_treehedge","connect500_treehedge","clumpy500_treehedge","ta500_treehedge","pacfrac500_treehedge")]

jpeg(paste0(wd.out,"/var_boxplots_t_frag.jpg"), width=600, height=800)
dat.frag %>%
  gather(-source, key = "var", value = "value") %>%
  ggplot(aes(x=as.factor(source), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  # facet_grid(vars(vars), scales = "free") + #factor(size, levels=c('50%','100%','150%','200%')
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) 
dev.off() ## ta100 and treehedge graphs should be identical, which they are, hooray!


####### CORRELATIONS #####
### NDMP
df.ndmp <- dat.plot[!(dat.plot$source %in% c("pseudo-absence", "LRD+NDD")),] ## NDMP presence and survey absence data only
table(df.ndmp$source) ## 49 presences vs 164 absences
df.ndmp <- df.ndmp[,c("occ", vars)]
# df <- df[, colSums(df)!=0]
df.ndmp <- na.omit(df.ndmp)

corrs <- cor(df.ndmp[-1], method = c("spearman"))
jpeg(paste0(wd.out,"/corrs_t_sampled_ndmp.jpg"), width=600, height=600)
corrplot(corrs)
dev.off()

### NDD+LRD
df.nddlrd <- dat.plot[!(dat.plot$source %in% c("pseudo-absence", "NDMP")),] ## NDD and LRD presence and survey absence data only
table(df.nddlrd$source) # 209 presnces vs 164 absences
df.nddlrd <- df.nddlrd[,c("occ", vars)]
# df <- df[, colSums(df)!=0]
df.nddlrd <- na.omit(df.nddlrd)

corrs <- cor(df.nddlrd[-1], method = c("spearman"))
jpeg(paste0(wd.out,"/corrs_t_sampled_nddlrd.jpg"), width=600, height=600)
corrplot(corrs)
dev.off()

### NDMP and NDD+LRD
df <- dat.plot[!(dat.plot$source %in% c("pseudo-absence")),] ## survey absence data only
table(df$source) # 258 presnces vs 164 absences
df <- df[,c("source", vars)]
# df <- df[, colSums(df)!=0]
df <- na.omit(df)

corrs <- cor(df.nddlrd[-1], method = c("spearman"))
jpeg(paste0(wd.out,"/corrs_t_sampled.jpg"), width=600, height=600)
corrplot(corrs)
dev.off()

#### Extract correlations > 0.7 and tabulate
corrs2 <- rcorr(as.matrix(df))
corrs2$r[lower.tri(corrs2$r)] <- NA
a <- as.data.frame(which(abs(corrs2$r) >= 0.7 & abs(corrs2$r) != 1, arr.ind = TRUE))
##  correlation >= 0.7
r <- rownames(corrs2$r)[a$row]
c <- colnames(corrs2$r)[a$col]
 
corr.out <- data.frame(matrix(nrow=nrow(a), ncol=2, dimnames=list(c(), c("vars", "corrs"))), stringsAsFactors=F)
corr.out$vars <- paste0(r, " & ", c)
corr.out$corrs <- apply(a, MARGIN=1, FUN=function(x) { corrs2$r[x[1], x[2]] }) ## 
# write.csv(corr.out, file=paste0(wd.out, "/corrs_t_sampled.csv"), row.names=F)


##### DECIDE WHICH CORRELATED VARIABLES TO RETAIN #####
### Merge NDMP and NDD+LRD as only 49 points for NDMP. Use source as explanatory variable and interactor. 
### Calibration and validation data
df$id <- 1:nrow(df)
df$occ <- 0
df$occ[df$source=="NDMP"] <- 1

dat.c <- sample_frac(df, size=0.7)
dat.v <- df[!(df$id %in% dat.c$id),]

### Variables to choose between - for variables that are measured in teh same seasons. 
for (c in c("tasmax")) { 
  vars <- paste0(c, c("_spr_5yrmn", "_smr_5yrmn", "_aut_5yrmn"))
  # vars <- c("rainfall_spr_5yrmn", "rainfall_smr_5yrmn", "rainfall_aut_5yrmn")
  # vars <- c("sun_spr_5yrmn", "sun_smr_5yrmn", "sun_aut_5yrmn")
  # vars <- c("rainfall_spr_5yrmn", "sun_spr_5yrmn", "tasmax_aut_5yrmn", "tasmax_smr_5yrmn", "tasmax_spr_5yrmn", "tasmin_win_5yrmn" )
  # vars <- c("rainfall_spr_5yrmn", "sun_spr_5yrmn", "tasmax_spr_5yrmn", "tasmin_win_5yrmn" )
  
  for (v in vars) {
    m <- glm(as.formula(paste0("occ ~ ", v)), family=binomial, dat.c) ## Make model with calibration data
    p <- predict(m, dat.v, re.form=NA, type="response") ## Predict model with validation data
    
    roc_obj <- roc(dat.v$occ, p)
    print(paste(v, round(auc(roc_obj),3)))
  }
}

vars <- c("rainfall_spr_5yrmn", "sun_spr_5yrmn", "tasmax_spr_5yrmn", "tasmin_win_5yrmn" )
corrs <- cor(df[,vars], method = c("spearman"))



