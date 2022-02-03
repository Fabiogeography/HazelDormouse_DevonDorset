##############################################################
##### Explore data before making SDM ########################
##### Written by: Regan Early ################################
##### Written on: 11th Jan 2021 ##############################
##### Modified on: 30th November 2021  #########################
##############################################################

library(raster, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(rgdal, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(pROC, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(tidyverse, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library") ## glimpse and useful functions
library(corrplot, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library")
library(Hmisc, lib.loc="D:\\SOFTWARE\\R-4.1.1\\library") ## rcorr

wd.env <- "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m"
wd.distn <- "E:/NON_PROJECT/DORMOUSE_NE/DISTRIBUTION"
wd.out <- "E:/NON_PROJECT/DORMOUSE_NE/SDM/EXPLORE"
wd.clim <- "E:/NON_PROJECT/DORMOUSE_WALES/GIS/100m" # "E:/GIS_DATA/CLIMATE/UK/HadUK-Grid_1km"
OSGB.proj <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'

## Devon Dorset area
dd <- readOGR(dsn="E:/NON_PROJECT/DORMOUSE_NE/GIS", layer="devondorset_buffer5km", stringsAsFactors = F)

##### OBTAIN DATA #####
##### Dormouse records #####
p <- readOGR(dsn=wd.distn, layer="p100m1990_dd", stringsAsFactors = F)
a <- readOGR(dsn=wd.distn, layer="a100m1990_dd", stringsAsFactors = F)
proj4string(p) <- proj4string(a) <- OSGB.proj ## Already in OSGB so doesn't transform, just ensures projection is written identically to other objects
 
## Change name of a column so it matches columns in p
a$IDorig <- a$OSGB_G_
a$OSGB_G_ <- NULL
 
## Add extra columns in p
x <- names(p)[!(names(p) %in% names(a))] ## Add columns in p to a
a@data[,x] <- NA
 
occ <- rbind(p, a)

## Assign year to absences.See dormouse Wales project for decision process here 
## Make the frequency of years assigned to absences match the frequency of years in presence data
summary(occ$year) ## Check what the earliest and latest years where dormouse recorded in DD and use below
test <- sample(c(2011:2019), table(occ$occ)[1], replace=TRUE, prob=table(occ$year))
par(mfrow=c(1,2))
hist(occ$year, main = "Year when presence most\nrecently recorded", xlab="year",
     breaks=seq(2011,2019,1))
hist(test, main="Year assigned\nto absence", xlab="year",
     breaks=seq(2011,2019,1))

occ$year[is.na(occ$year)] <- test; rm(test)

occ <- occ[occ$year<2020,] ## remove 2020 data - it's annoying!
 
writeOGR(occ, dsn=wd.distn, layer="occ100m1990_dd", driver="ESRI Shapefile")
occ <- readOGR(dsn=wd.distn, layer="occ100m1990_dd") ## load, so always working with the same dataset

##### HABITAT DATA #####
### Proportion of each 1km grid-cell covered by the respective habitat. 
## National Forest Product. England (approximately) and Wales. Note sea has value 0
nfi.layers <- c("broadl", "conif", "coppice", "copp_stds", "yngtree", "lodens", "mix_con", "mix_bdl", "shrub", "felled")
nfi.m <- stack(paste0(wd.env,"/",nfi.layers, "_m.tif")) ## Mean habitat cover between 2011 and 2017 inclusive. Note sea has 0 value
nfi.c <- stack(paste0(wd.env,"/",nfi.layers, "_c.tif")) ## Change in habitat cover between 2011 and 2017 inclusive. Note sea has 0 value

## Ancient woodland. Note sea has value 0
anc.layers <- c("anc_rep_wood", "anc_sn_wood")
anc <- stack(paste0(wd.env,"/",anc.layers,".tif")) ## Think these data were collected between 2006 and 2021 inclusive

## LCM
lcm <- expand.grid(c(1990,2015,2017,2019), c("_broadl", "_conif"))
lcm <- paste0("lcm_", lcm$Var1, lcm$Var2)
lcm <- stack(paste0(wd.env,"/", lcm, ".tif")) ## Mean habitat cover between 2011 and 2017 inclusive

##### CLIMATE DATA #####
### This code is only for checking that all occupancy records fall in a cell with climate data
### Create mask of occupied grid-cells for each year, multiply each year raster by the mask to nullify cells with no records that year, and add all rasters together
clim.vars <- c("rainfall_spr", "rainfall_smr", "rainfall_aut", "sun_spr", "sun_smr", "sun_aut", "tasmax_spr", "tasmax_smr", "tasmax_aut", "tasmin_win", "tasrng_win")

## Restrict to extent of occupancy data, to save time
e <- extent(dd) ## geographic extent of the study region

r <- raster(paste0(wd.clim, "/rainfall_spr_5yrmn_1990.tif")) ## These were created for the dormouse Wales project
r <- crop(r, e)
r <- (r/r) - 1 ## make all values 0

for(i in min(occ$year):max(occ$year)) { ## Note that although climate data are available from 2020 now, there are almost no dormouse records for 2020. 
  occ.yr <- occ[occ$year==i,]
  
  r.yr <- r ## make the mask grid
  r.yr[occ.yr] <- 1 ## set grid-cells in mask grid that have a presence or absence in the given year to be 1

  s <- stack(paste0(wd.clim, "/",clim.vars,"_5yrmn_",i,".tif"))
  s <- crop(s, e)
  s.yr <- s * r.yr ## Only grid-cells that contain presence or absence data in the given year retain a value. Others are given 0. 
  names(s.yr) <- substr(names(s), 1, nchar(names(s))-5)
  assign(paste0("clim.",i), s.yr)

}

### Adds the climate at each grid-cell where species is present or absent. problem could arise if value of any variable is precisely 0, as this is used as the numm value. Likelihood of this is low in this study region, so impact low. 
clim <- get(paste0("clim.", min(occ$year)))
for (y in min(occ$year)+1:max(occ$year)) { 
  a <- get(paste0("clim.",y))
  # a <- crop(a, extent(w))
  clim <- clim + a
  print(y)
}

names(clim) <- substr(names(s), 1, nchar(names(s))-5)

# writeRaster(clim, "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m/clim100m_dd")
clim <- stack("E:/NON_PROJECT/DORMOUSE_NE/GIS/100m/clim100m_dd")

##### TOPOGRAPHY DATA #####
## Aspect cannot not used because it is calculated as circular degrees clockwise from 0 to 360 degrees, and it is therefore difficult to compare quantitatively (e.g. there is only a 2u difference between aspects of 1u and 359u whereas the numerical difference is 358u). The method suggested by Copland (1998) was used to perform the conversion from aspect to NS and WE, as follows: NS=cos(aspect) and WE=sin(aspect) such that values of NS and WE range from -1 to 1 and represent the extent to which a slope faces north (NS=1), south (NS=-1), east (WE=1), or west
## m(WE=-1). (https://www.researchgate.net/figure/DEM-resolutions5m-25m-50m-and-100m-and-statistical-variation-of-terrain_fig3_220650484)

aspNS <- raster(paste0(wd.env,"/OS_Terrain_100_NS.tif")) ## Aspect. Extent to which a slope faces north (NS=1) or south (NS=-1), 
aspWE <- raster(paste0(wd.env,"/OS_Terrain_100_WE.tif")) ## Aspect. Extent to which a slope faces east (WE=1) or west (WE=-1), 
slope <- raster(paste0(wd.env,"/OS_Terrain_100_slope_pct.tif")) ## Slope in %
crs(aspNS) <- crs(aspWE) <- crs(slope) <- crs(anc) ## Slight discrepancy in k value means need to adjust ythe CRS to match hte other rasters to be stacked

## Crop other environmental data to extent (LCM has slightly lower ymin than others, so crop all to anc extent)
anc <- crop(anc, e)
lcm <- crop(lcm, extent(anc))
nfi.m <- crop(nfi.m, extent(anc))
nfi.c <- crop(nfi.c, extent(anc))
aspNS <- crop(aspNS, extent(anc))
aspWE <- crop(aspWE, extent(anc))
slope <- crop(slope, extent(anc))

##### Unite all environmental data and extract to data.frame #####
env.r <- stack(nfi.m, nfi.c, anc, lcm, aspNS, aspWE, slope, clim)

env.r[["broadl_all_m"]] <- env.r[["mix_bdl_m"]] + env.r[["broadl_m"]]
env.r[["conif_all_m"]] <- env.r[["mix_con_m"]] + env.r[["conif_m"]]

env.r[["lcm_broadl_m"]] <- mean(env.r[["lcm_1990_broadl"]], env.r[["lcm_2015_broadl"]], env.r[["lcm_2017_broadl"]], env.r[["lcm_2019_broadl"]])
env.r[["lcm_broadl_c"]] <- env.r[["lcm_2019_broadl"]] - env.r[["lcm_1990_broadl"]]
env.r[["lcm_conif_m"]] <- mean(env.r[["lcm_1990_conif"]], env.r[["lcm_2015_conif"]], env.r[["lcm_2017_conif"]], env.r[["lcm_2019_conif"]])
env.r[["lcm_conif_c"]] <- env.r[["lcm_2019_conif"]] - env.r[["lcm_1990_conif"]]

env.r[["combo_broadl_m"]] <- env.r[["lcm_broadl_m"]] + env.r[["broadl_all_m"]] ## sum NFI and LCM forest
writeRaster(env.r[["combo_broadl_m"]], paste0(wd.env, "/combo_broadl_m.tif"))
env.r[["combo_conif_m"]] <- env.r[["lcm_conif_m"]] + env.r[["conif_all_m"]] ## sum NFI and LCM forest
writeRaster(env.r[["combo_conif_m"]], paste0(wd.env, "/combo_conif_m.tif"))

env.r[["combo_broadl_c"]] <- min(env.r[["lcm_broadl_c"]], env.r[["broadl_c"]]) ## the most negative forest change
writeRaster(env.r[["combo_broadl_c"]], paste0(wd.env, "/combo_broadl_c.tif"))
env.r[["combo_conif_c"]] <- min(env.r[["lcm_conif_c"]], env.r[["conif_c"]]) ## the most negative forest change
writeRaster(env.r[["combo_conif_c"]], paste0(wd.env, "/combo_conif_c.tif"))

env.r[["anc_wood"]] <- env.r[["anc_rep_wood"]] + env.r[["anc_sn_wood"]] ## the most negative forest change

env.r[["broadl_all_c"]] <- env.r[["mix_bdl_c"]] + env.r[["broadl_c"]]
env.r[["conif_all_c"]] <- env.r[["mix_con_c"]] + env.r[["conif_c"]]

# writeRaster(env.r, paste0(wd.env, "/env_100m7_dd_16thDec2021")) ## Values are different to the ones I have on file, below
env.r <- stack(paste0(wd.env, "/env_100m7_dd")) 

env <- cbind(occ, raster::extract(env.r, occ)) ## spatial points data frame with env data. Note 'extract' is also a function in tidyr
# the above line doesn't work when clim is attached
dat <- as.data.frame(env) ## Make dataframe for easier modelling
glimpse(dat) 
summary(dat) ## 2 records fall outside lcm range
dat <- dat[!is.na(dat$lcm_broadl_m),]
dat <- dat[!is.na(dat$tasrng_win_5yrmn),]

# write.csv(dat, paste0(wd.env,"/dat_30Nov2021_dd.csv"), row.names=F)
dat <- read.csv(paste0(wd.env,"/dat_30Nov2021_dd.csv"))

####### visualise ################
### Plot division of all numeric variables between presence and absence points
jpeg(paste0(wd.out,"/var_boxplots1_dd.jpg"), width=600, height=400)
dat.1 <- cbind(dat$source, dat[,grepl("_m", colnames(dat))])
colnames(dat.1)[1] <- "occ"
dat.1 %>%
  gather(-occ, key = "var", value = "value") %>%
  ggplot(aes(x=as.factor(occ), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
dev.off()

jpeg(paste0(wd.out,"/var_boxplots2_dd.jpg"), width=600, height=400)
dat.2 <- cbind(dat$source, dat[,grepl("_c", colnames(dat))])
colnames(dat.2)[1] <- "occ"
dat.2 %>%
  gather(-occ, key = "var", value = "value") %>%
  ggplot(aes(x=as.factor(occ), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
dev.off()

jpeg(paste0(wd.out,"/var_boxplots3_dd.jpg"), width=600, height=400)
dat.3 <- cbind(dat$source, dat[,c("lcm_1990_broadl", "lcm_2015_broadl", "lcm_2017_broadl", "lcm_2019_broadl", "lcm_1990_conif", "lcm_2015_conif", "lcm_2017_conif", "lcm_2019_conif")])
colnames(dat.3)[1] <- "occ"
dat.3 %>%
  gather(-occ, key = "var", value = "value") %>%
  ggplot(aes(x=as.factor(occ), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
dev.off()

dat.4 <- cbind(dat$source, dat[,c("OS_Terrain_100_NS", "OS_Terrain_100_WE", "OS_Terrain_100_slope_pct",
                                  "anc_rep_wood", "anc_sn_wood")])
colnames(dat.4) <- c("occ", "Aspect NS", "Aspect WE", "Slope %",
                     "Ancient replanted wood", "Ancient semi-nat wood")
dat.4$occ[is.na(dat.4$occ)] <- "PA"
jpeg(paste0(wd.out,"/var_boxplots4_dd.jpg"), width=600, height=400)
dat.4 %>%
  gather(-occ, key = "var", value = "value", factor_key = T) %>%
  ggplot(aes(x=as.factor(occ), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
dev.off()

jpeg(paste0(wd.out,"/var_boxplots5_dd.jpg"), width=600, height=400)
dat.5 <- cbind(dat$source, dat[,grepl("_5yrmn", colnames(dat))])
colnames(dat.5)[1] <- "occ"
dat.5 %>%
  gather(-occ, key = "var", value = "value") %>%
  ggplot(aes(x=as.factor(occ), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
dev.off()

jpeg(paste0(wd.out,"/var_boxplots6_dd.jpg"), width=600, height=400)
dat.6 <- dat[,c("source", "lcm_1990_broadl", "lcm_2015_broadl", "lcm_2017_broadl", "lcm_2019_broadl", "lcm_broadl_m", "lcm_broadl_c",
                "lcm_1990_conif", "lcm_2015_conif", "lcm_2017_conif", "lcm_2019_conif", "lcm_conif_m", "lcm_conif_c")]
dat.6 %>%
  gather(-source, key = "var", value = "value") %>%
  ggplot(aes(x=as.factor(source), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
dev.off()

## How many grid-cells with mixed vs 'full' forest?
table(dat$broadl_m>0) ## 881 without, 359 with
table(dat$broadl_all_m>0) ## 865 without, 375 with. Little change. Supports grouping mixed into all broadl forest.
table(dat$conif_m>0) ## 1142 without, 98 with
table(dat$conif_all_m>0) ## 1128 without, 112 with. Little change. Supports grouping mixed into all broadl forest.

## Investigate number of cells with env variable values I'm not sure about
table(dat$yngtree_m>0) ## 1205 without, 35 with
table(dat$lodens_m>0) ## 1235 without, 5 with
table(dat$shrub_m>0) ## 1232 without, 8 with
table(dat$felled_m>0) ## 1218 without, 22 with

## Investigate change in broadleaved and coniferous
jpeg(paste0(wd.out,"/var_boxplots7_dd.jpg"), width=600, height=400)
dat.7 <- dat[,c("source", "broadl_all_c", "conif_all_c", "combo_broadl_c", "combo_conif_c")] 
dat.7 %>%
  gather(-source, key = "var", value = "value") %>%
  ggplot(aes(x=as.factor(source), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
dev.off() 

table(dat$broadl_c>0) ## 1215 without, 25 with
table(dat$conif_c>0) ## 1236 without, 4 with

## Check out non-NFI vegetation types
table(dat$anc_rep_wood>0) ## 1204 without, 36 with
table(dat$anc_sn_wood>0) ## 1179 without, 61 with

####### CORRELATIONS #####
### Forest. LCM and NFI combos not included yet
df <- dat[,c("occ", "broadl_all_m", "broadl_c", "conif_all_m", "conif_c", "lcm_broadl_m", "lcm_broadl_c", "lcm_conif_m", "lcm_conif_c", "anc_rep_wood", "anc_sn_wood")]
df <- df[, colSums(df)!=0]

corrs.forest <- cor(df, method = c("spearman"))
jpeg(paste0(wd.out,"/corrs_forest_dd.jpg"), width=600, height=600)
corrplot(corrs.forest)
dev.off()

#### Extract correlations > 0.7 and tabulate
corrs.forest2 <- rcorr(as.matrix(df))
corrs.forest2$r[lower.tri(corrs.forest2$r)] <- NA
a <- as.data.frame(which(abs(corrs.forest2$r) >= 0.7 & abs(corrs.forest2$r) != 1, arr.ind = TRUE))
## Only one correlation >= 0.7
r <- rownames(corrs.forest2$r)[a$row]
c <- colnames(corrs.forest2$r)[a$col]
 
corr.out <- data.frame(matrix(nrow=nrow(a), ncol=2, dimnames=list(c(), c("vars", "corrs"))), stringsAsFactors=F)
corr.out$vars <- paste0(r, " & ", c)
corr.out$corrs <- apply(a, MARGIN=1, FUN=function(x) { corrs.forest2$r[x[1], x[2]] }) ## broadl_all_m & lcm_broadl_m 0.84 and conif_all_m & lcm_conif_m 0.79
# write.csv(corr.out, file=paste0(wd.out, "/corrs_forest_dd.csv"), row.names=F)

### CLIMATE
df <- dat[,grepl("_5yrmn", colnames(dat))]
df <- df[, colSums(na.omit(df))!=0]
df <- na.omit(df)

corrs.climate <- cor(df, method = c("spearman"))
jpeg(paste0(wd.out,"/corrs_climate_dd.jpg"), width=600, height=600)
corrplot(corrs.climate)
dev.off()

#### Extract correlations > 0.7 and tabulate
corrs.climate2 <- rcorr(as.matrix(df))
corrs.climate2$r[lower.tri(corrs.climate2$r)] <- NA
a <- as.data.frame(which(abs(corrs.climate2$r) >= 0.7 & abs(corrs.climate2$r) != 1, arr.ind = TRUE))
r <- rownames(corrs.climate2$r)[a$row]
c <- colnames(corrs.climate2$r)[a$col]

corr.out <- data.frame(matrix(nrow=nrow(a), ncol=2, dimnames=list(c(), c("vars", "corrs"))), stringsAsFactors=F)
corr.out$vars <- paste0(r, " & ", c)
corr.out$corrs <- apply(a, MARGIN=1, FUN=function(x) { corrs.climate2$r[x[1], x[2]] })
# write.csv(corr.out, file=paste0(wd.out, "/corrs_climate_dd.csv"), row.names=F)

##### DECIDE WHICH CORRELATED VARIABLES TO RETAIN #####
### Calibration and validation data
dat.c <- sample_frac(dat, size=0.7)
dat.v <- dat[!(dat$IDorig %in% dat.c$IDorig),]

### Variables to choose between - for vasriables that are measured in teh same seasons. 
for (c in c("tasmax", "sun")) { 
  vars <- paste0(c, c("_spr_5yrmn", "_smr_5yrmn", "_aut_5yrmn"))
  # vars <- c("rainfall_spr_5yrmn", "rainfall_smr_5yrmn", "rainfall_aut_5yrmn")
  # vars <- c("sun_spr_5yrmn", "sun_smr_5yrmn", "sun_aut_5yrmn")
  
  for (v in vars) {
    m <- glm(as.formula(paste0("occ ~ ", v)), family=binomial, dat.c) ## Make model with calibration data
    p <- predict(m, dat.v, re.form=NA, type="response") ## Predict model with validation data
    
    roc_obj <- roc(dat.v$occ, p)
    print(paste(v, round(auc(roc_obj),3)))
  }
}

### Variables to choose between - for vasriables that are not measured in teh same seasons. 
c <- c("tasmax", "tasrng")
vars <- c(paste0(c[1], c("_spr_5yrmn", "_smr_5yrmn", "_aut_5yrmn")), 
          paste0(c[2], c("_win_5yrmn")))

for (v in vars) {
  m <- glm(as.formula(paste0("occ ~ ", v)), family=binomial, dat.c) ## Make model with calibration data
  p <- predict(m, dat.v, re.form=NA, type="response") ## Predict model with validation data
  
  roc_obj <- roc(dat.v$occ, p)
  print(paste(v, round(auc(roc_obj),3)))
}

##### CHECK FOR CORRELATIONS BETWEEN CLIMATE AND FOREST VARIABLES
## Version with NFI and LCM seperate
df <- dat[,c("occ", "broadl_all_m", "broadl_c", "conif_all_m", "conif_c", 
             "lcm_broadl_m", "lcm_broadl_c", "lcm_conif_m", "lcm_conif_c",
             "rainfall_smr_5yrmn",	"sun_aut_5yrmn", "tasmin_win_5yrmn", "tasmax_aut_5yrmn", "tasrng_win_5yrmn",
             "OS_Terrain_100_NS", "OS_Terrain_100_WE", "OS_Terrain_100_slope_pct")]

## Version with the combined data from NFI and LCM
# df <- dat[,c("occ", "combo_broadl_m", "combo_broadl_c", "combo_conif_m", "combo_conif_c",
#              "rainfall_aut_5yrmn",	"sun_spr_5yrmn", "tasmin_win_5yrmn",
#              "OS_Terrain_100_NS", "OS_Terrain_100_WE", "OS_Terrain_100_slope_pct")]
# df <- df[, colSums(na.omit(df))!=0]
# df <- na.omit(df)

corrs <- cor(df, method = c("spearman"))
jpeg(paste0(wd.out,"/corrs_climforestterrain_dd.jpg"), width=600, height=600)
corrplot(corrs)
dev.off()

#### Extract correlations > 0.7 and tabulate
corrs <- rcorr(as.matrix(df))
corrs$r[lower.tri(corrs$r)] <- NA
a <- as.data.frame(which(abs(corrs$r) >= 0.7 & abs(corrs$r) != 1, arr.ind = TRUE))

## No correlations >= 0.7 
r <- rownames(corrs$r)[a$row]
c <- colnames(corrs$r)[a$col]

corr.out <- data.frame(matrix(nrow=nrow(a), ncol=2, dimnames=list(c(), c("vars", "corrs"))), stringsAsFactors=F)
corr.out$vars <- paste0(r, " & ", c)
corr.out$corrs <- apply(a, MARGIN=1, FUN=function(x) { corrs$r[x[1], x[2]] }) ## only the nfi and lcn broadleaved and coniferous forest as before

jpeg(paste0(wd.out,"/var_boxplots_selected_dd.jpg"), width=600, height=400)
dat.8 <- dat[,c("source", "broadl_all_m", "broadl_c", "conif_all_m", "conif_c", 
                "lcm_broadl_m", "lcm_broadl_c", "lcm_conif_m", "lcm_conif_c",
                "rainfall_smr_5yrmn",	"sun_aut_5yrmn", "tasmin_win_5yrmn", "tasmax_aut_5yrmn", "tasrng_win_5yrmn",
                "OS_Terrain_100_NS", "OS_Terrain_100_WE", "OS_Terrain_100_slope_pct")] 
dat.8 %>%
  gather(-source, key = "var", value = "value") %>%
  ggplot(aes(x=as.factor(source), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
dev.off()

##### For report #####

vars <- c("coppice_m","copp_stds_m","yngtree_m","lodens_m","shrub_m","felled_m",
          "coppice_c","copp_stds_c","yngtree_c","lodens_c","shrub_c","felled_c",
          "broadl_m","mix_bdl_m","broadl_all_m","conif_m","mix_con_m","conif_all_m",
          "broadl_c","mix_bdl_c","broadl_all_c","conif_c","mix_con_c","conif_all_c",
          "lcm_broadl_m","lcm_conif_m","lcm_broadl_c","lcm_conif_c",
          "combo_broadl_m","combo_conif_m","combo_broadl_c","combo_conif_c")

var.names <- c("NFI coppice", "NFI coppice stds", "NFI yngtree", "NFi lodens", "NFI shrub", "NFI felled",
               "NFI coppice change", "NFI coppice stds change", "NFI yngtree change", "NFi lodens change", "NFI shrub change", "NFI felled change",
               "NFI broadl", "NFI broadl mix", "NFI broadl+mix", "NFI conif", "NFI conif mix", "NFI conif+mix",
               "NFI broadl change", "NFI broadl mix change", "NFI broadl+mix change", "NFI conif change", "NFI conif mix change", "NFI conif+mix change",
               "LCM broadl", "LCM conif", "LCM broadl change", "LCM conif change", 
               "LCM NFI broadl+mix", "LCM NFI conif+mix", "LCM NFI broadl+mix change", "LCM NFI conif+mix change"
                )

dat.1 <- cbind(dat$source, dat[,vars])

## Convert combined NFI and LCM variables to averages
dat.1$combo_broadl_m <- dat.1$combo_broadl_m/2
dat.1$combo_conif_m <- dat.1$combo_conif_m/2
dat.1$combo_broadl_c <- dat.1$combo_broadl_c/2
dat.1$combo_conif_c <- dat.1$combo_conif_c/2

colnames(dat.1) <- c("occ", var.names)
dat.1$occ[is.na(dat.1$occ)] <- "PA"
jpeg(paste0(wd.out,"/var_boxplots_forest_report_dd_new.jpg"), width=1000, height=600)
dat.1 %>%
  gather(-occ, key = "var", value = "value", factor_key=T) %>% ## factor_key means order of columns is preserved
  ggplot(aes(x=as.factor(occ), y=value)) +
  geom_boxplot() +
  # ylim(-1,1) +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
dev.off()

