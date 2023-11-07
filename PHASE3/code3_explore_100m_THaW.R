##############################################################
##### Explore data before making SDM ########################
##### Written by: Regan Early ################################
##### Written on: 11th Jan 2021 ##############################
##### Modified on: 21 feb 2022  #########################
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

writeRaster(clim, "E:/NON_PROJECT/DORMOUSE_NE/GIS/100m/clim100m_t_surveyabsences")
clim <- stack("E:/NON_PROJECT/DORMOUSE_NE/GIS/100m/clim100m_t_surveyabsences")

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

## Save raster stack
writeRaster(env.r, paste0(wd.env, "/env_100m7_t")) 
env.r <- stack(paste0(wd.env, "/env_100m7_t")) 

##### Extract the environmental data at the presences, original pseudo-absences and survey absences #####
### Created in PHASE 2
occ_dd <- readOGR(dsn=wd.distn, layer="occ100m1990_dd") ## Note that does not have m100 values for pseudo-absences
occ_dd <- occ_dd[!(is.na(occ_dd$m100) & occ_dd$occ==1),] ## remove the once presence without a m100 value
occ_dd$source[occ_dd$occ==0] <- "pseudo-absence"

## Add extra columns to occ_dd
x <- names(occ_t)[!(names(occ_t) %in% names(occ_dd))] ## Add columns in p to a
occ_dd@data[,x] <- NA

## Merge dd and t data
occ <- rbind(occ_t, occ_dd)

env <- cbind(occ, raster::extract(env.r, occ)) ## spatial points data frame with env data. Note 'extract' is also a function in tidyr
dat <- as.data.frame(env) ## Make dataframe for easier modelling
glimpse(dat) 
summary(dat) ## 690 records fall outside THaW range. 

## Testing the 690 records with NA in env values are indeed outside the ThaW region
test <- dat[is.na(dat$lcm_broadl_m),]
plot(test$easting, test1$northing)
plot(env.r[["lcm_broadl_m"]], add=T)

## Remove the 690 records
dat <- dat[!is.na(dat$pacfrac500_treehedge),] 

# write.csv(dat, paste0(wd.env,"/dat_15Feb2022_t.csv"), row.names=F)
dat <- read.csv(paste0(wd.env,"/dat_15Feb2022_t.csv"))

####### visualise ################
## Add the climate values for the ThaW survey absences to the plot data. I didn't do this in the dat dataframe, as some pseudo-absences might overlap with the survey absences, and I didn't want to lose the data for the pseudo-absenes
clim.vars <- c('tasmin_win_5yrmn', 'tasrng_win_5yrmn', 'tasmax_spr_5yrmn', 'sun_spr_5yrmn',
               'rainfall_spr_5yrmn', 'tasmax_smr_5yrmn', 'sun_smr_5yrmn', 'rainfall_smr_5yrmn', 'tasmax_aut_5yrmn', 'sun_aut_5yrmn', 'rainfall_aut_5yrmn')
dat.plot <- dat.old

for (n in clim.vars) {
  # n <- clim.vars[1]
  nt <- paste0(n,"t")
  dat.plot[dat.plot$source %in% c("PTES", "Non-PTES"), n] <- dat.plot[dat.plot$source %in% c("PTES", "Non-PTES"), nt]
}

## Temporarily remove all survey absences from 2020 and 2021 to remove 0s from climate graphs, which prevent us seeing full range of data. Should include climate data for them
dat.plot <- dat.plot[dat.plot$year<2020,]

## Plot division of all numeric variables between the different types of presence and absence points
vars <- c('combo_broadl_m', 'combo_broadl_c', 'combo_conif_m', 'combo_conif_c', 'anc_wood', 'OS_Terrain_100_NS', 'OS_Terrain_100_WE',
         'OS_Terrain_100_slope_pct', 'treehedge', 'scrub', 'tasmin_win_5yrmn', 'tasrng_win_5yrmn', 'tasmax_spr_5yrmn', 'sun_spr_5yrmn',
         'rainfall_spr_5yrmn', 'tasmax_smr_5yrmn', 'sun_smr_5yrmn', 'rainfall_smr_5yrmn', 'tasmax_aut_5yrmn', 'sun_aut_5yrmn', 'rainfall_aut_5yrmn')

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

# jpeg(paste0(wd.out,"/var_boxplots_t.jpg"), width=600, height=800)
dat.plot %>%
  gather(-source, key = "var", value = "value") %>%
  ggplot(aes(x=as.factor(source), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  # facet_grid(vars(vars), scales = "free") + #factor(size, levels=c('50%','100%','150%','200%')
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) 
dev.off()

plot(dat.plot$year, dat.plot$sun_aut_5yrmn)
plot(dat.plot$year, dat.plot$tasmax_aut_5yrmn)

dat.plot.2010 <- dat.plot[dat.plot$year>=2010,]
jpeg(paste0(wd.out,"/var_boxplots_t_2010.jpg"), width=600, height=800)
dat.plot.2010 %>%
  gather(-source, key = "var", value = "value") %>%
  ggplot(aes(x=as.factor(source), y=value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  # facet_grid(vars(vars), scales = "free") + #factor(size, levels=c('50%','100%','150%','200%')
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) 
dev.off()

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

