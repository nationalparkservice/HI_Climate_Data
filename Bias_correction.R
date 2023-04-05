library(stars)
library(raster)
library(dplyr)
rm(list=ls())

monthly.index<-read.csv(paste0(here::here('data/Monthly_indexing.csv')))
Obs.RF.dir <- "D:/HI_Data/OBS/RFMonthYr_Rasters_Ma_in_1920_2012/Month_Rasters_Ma_in/"
Obs.Tmean.dir <- "D:/HI_Data/OBS/Tair_month_raster/"

scens<-c("present","rcp45","rcp85")
RFpresent <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'RF.monthly-',scens[1]))
RFRCP45 <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'RF.monthly-',scens[2]))
RFRCP85 <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'RF.monthly-',scens[3]))

boundary <- st_read('./data/HALE/HALE_boundary.shp')
boundary <- st_transform(boundary, st_crs(RFpresent))

#Annual Precip - ANNPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"ANN_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent, c("x", "y"), FUN=function(x) mean(x)*12)
rcp45 <- st_apply(RFRCP45, c("x", "y"), FUN=function(x) mean(x)*12)
rcp85 <- st_apply(RFRCP85, c("x", "y"), FUN=function(x) mean(x)*12)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'ANNPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'ANNPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'ANNPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'ANNPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'ANNPrecipBC_rcp85'))

#Jan Precip - JanPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"JAN_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,1])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,1])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,1])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'JanPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JanPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JanPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JanPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JanPrecipBC_rcp85'))

#Feb Precip - FebPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"FEB_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,2])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,2])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,2])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'FebPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'FebPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'FebPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'FebPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'FebPrecipBC_rcp85'))

#Mar Precip - MarPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"MAR_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,3])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,3])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,4])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'MarPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'MarPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'MarPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'MarPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'MarPrecipBC_rcp85'))


#Apr Precip - AprPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"APR_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,4])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,4])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,4])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'AprPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'AprPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'AprPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'AprPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'AprPrecipBC_rcp85'))


#May Precip - MayPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"MAY_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,5])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,5])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,5])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'MayPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'MayPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'MayPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'MayPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'MayPrecipBC_rcp85'))


#Jun Precip - JunPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"JUN_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,6])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,6])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,6])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'JunPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JunPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JunPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JunPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JunPrecipBC_rcp85'))


#Jul Precip - JulPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"JUL_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,7])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,7])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,7])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'JulPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JulPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JulPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JulPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JulPrecipBC_rcp85'))


#Aug Precip - AugPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"AUG_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,8])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,8])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,8])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'AugPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'AugPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'AugPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'AugPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'AugPrecipBC_rcp85'))


#Sep Precip - SepPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"SEP_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,9])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,9])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,9])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'SepPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'SepPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'SepPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'SepPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'SepPrecipBC_rcp85'))


#Oct Precip - OctPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"OCT_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,10])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,10])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,10])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'OctPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'OctPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'OctPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'OctPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'OctPrecipBC_rcp85'))


#Nov Precip - NovPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"NOV_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,11])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,11])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,11])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'NovPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'NovPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'NovPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'NovPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'NovPrecipBC_rcp85'))


#Dec Precip - DecPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"DEC_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent[,,,c(monthly.index[,12])],c("x","y"),mean)
rcp45 <- st_apply(RFRCP45[,,,c(monthly.index[,12])],c("x","y"),mean)
rcp85 <- st_apply(RFRCP85[,,,c(monthly.index[,12])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'DecPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'DecPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'DecPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'DecPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'DecPrecipBC_rcp85'))

rm(RFpresent,RFRCP45,RFRCP85,present,rcp45,rcp85,BC_rcp45,BC_rcp85,Delta_rcp45,Delta_rcp85,Obs)




"tair_ann"
#Extract January data

P<-readRDS(paste0(here::here('data/Output/Data-files//'),'RF.monthly-present'))
R45<-readRDS(paste0(here::here('data/Output/Data-files//'),'RF.monthly-rcp45'))

#Annual Precip - ANNPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"ANN_mean_1990-2009.tif")); Obs<- setNames(Obs, "inch") #Change input
st_crop(Obs,boundary) %>%  
  mutate(mean = inch*25.4) %>% select(mean) -> Obs

present <- st_apply(RFpresent, c("x", "y"), FUN=function(x) mean(x)*12)
rcp45 <- st_apply(RFRCP45, c("x", "y"), FUN=function(x) mean(x)*12)
rcp85 <- st_apply(RFRCP85, c("x", "y"), FUN=function(x) mean(x)*12)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'ANNPrecip-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'ANNPrecipDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'ANNPrecipDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'ANNPrecipBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'ANNPrecipBC_rcp85'))




