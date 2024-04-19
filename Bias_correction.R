library(stars)
library(raster)
library(dplyr)
rm(list=ls())

monthly.index<-read.csv(paste0(here::here('data/Monthly_indexing.csv')))
Obs.RF.dir <- "D:/HI_Data/OBS/StateRFGrids_mm/"
Obs.Tmean.dir <- "D:/HI_Data/OBS/Tair_month_raster/"

scens<-c("present","rcp45","rcp85")
RFpresent <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'RF.monthly-',scens[1]))
RFRCP45 <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'RF.monthly-',scens[2]))
RFRCP85 <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'RF.monthly-',scens[3]))

# boundary <- st_read('./data/HALE/HALE_boundary.shp')
boundary <- st_read("C:/Users/arunyon/OneDrive - DOI/Documents/GIS/HAVO_Kilauea_Summit_Wet_Dry_Zones/HAVO_Kilauea_Summit_Wet_Dry_Zones.shp")
boundary <- st_transform(boundary, st_crs(RFpresent))

#Annual Precip - ANNPrecip
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mmann")) ; Obs<- setNames(Obs, "mean") #Units in mm -- same as projections; #Change input to match projections
Obs <- st_crop(Obs,boundary)

# # Code from HALE -- to change names and units -- not needed for state raster, b/c in mm
# Obs<- setNames(Obs, "inch") #Change input
# st_crop(Obs,boundary) %>%  
#   mutate(mean = inch*25.4) %>% select(mean) -> Obs

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm01")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm02")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm03")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm04")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm05")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm06")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm07")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm08")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm09")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm10")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm11")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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
Obs <- read_stars(paste0(Obs.RF.dir,"staterf_mm12")) ; Obs<- setNames(Obs, "mean") 
Obs <- st_crop(Obs,boundary)

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


##################
### Tmean

Tmeanpresent <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'Tmean.monthly-',scens[1]))
TmeanRCP45 <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'Tmean.monthly-',scens[2]))
TmeanRCP85 <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'Tmean.monthly-',scens[3]))

#Annual Tmean - ANNTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_ann")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs # projections in TmeanF

present <- st_apply(Tmeanpresent, c("x", "y"), FUN=function(x) mean(x))
rcp45 <- st_apply(TmeanRCP45, c("x", "y"), FUN=function(x) mean(x))
rcp85 <- st_apply(TmeanRCP85, c("x", "y"), FUN=function(x) mean(x))

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'ANNTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'ANNTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'ANNTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'ANNTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'ANNTmeanBC_rcp85'))

#Jan Tmean - JanTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_jan")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,1])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,1])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,1])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'JanTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JanTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JanTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JanTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JanTmeanBC_rcp85'))


#Feb Tmean - FebTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_feb")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,2])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,2])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,2])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'FebTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'FebTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'FebTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'FebTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'FebTmeanBC_rcp85'))


#Mar Tmean - MarTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_mar")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,3])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,3])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,3])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'MarTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'MarTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'MarTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'MarTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'MarTmeanBC_rcp85'))


#Apr Tmean - AprTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_apr")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,4])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,4])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,4])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'AprTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'AprTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'AprTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'AprTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'AprTmeanBC_rcp85'))


#May Tmean - MayTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_may")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,5])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,5])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,5])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'MayTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'MayTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'MayTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'MayTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'MayTmeanBC_rcp85'))


#Jun Tmean - JunTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_jun")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,6])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,6])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,6])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'JunTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JunTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JunTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JunTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JunTmeanBC_rcp85'))


#Jul Tmean - JulTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_jul")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,7])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,7])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,7])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'JulTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JulTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JulTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'JulTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'JulTmeanBC_rcp85'))


#Aug Tmean - AugTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_aug")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,8])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,8])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,8])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'AugTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'AugTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'AugTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'AugTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'AugTmeanBC_rcp85'))


#Sep Tmean - SepTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_sep")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,9])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,9])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,9])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'SepTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'SepTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'SepTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'SepTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'SepTmeanBC_rcp85'))


#Oct Tmean - OctTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_oct")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,10])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,10])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,10])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'OctTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'OctTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'OctTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'OctTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'OctTmeanBC_rcp85'))


#Nov Tmean - NovTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_nov")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,11])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,11])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,11])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'NovTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'NovTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'NovTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'NovTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'NovTmeanBC_rcp85'))


#Dec Tmean - DecTmean
Obs <- read_stars(paste0(Obs.Tmean.dir,"tair_dec")); Obs<- setNames(Obs, "TmeanC") #Change input
st_crop(Obs,boundary) %>%  
  mutate(TmeanF = (TmeanC*9/5)+32) %>% select(TmeanF) -> Obs

present <- st_apply(Tmeanpresent[,,,c(monthly.index[,12])],c("x","y"),mean)
rcp45 <- st_apply(TmeanRCP45[,,,c(monthly.index[,12])],c("x","y"),mean)
rcp85 <- st_apply(TmeanRCP85[,,,c(monthly.index[,12])],c("x","y"),mean)

Obs<- st_warp(Obs, present)
saveRDS(Obs, file = paste0(here::here('data/Output/Data-files//'),'DecTmean-Obs')) #Change this line

Delta_rcp45 <- rcp45 - present
Delta_rcp85 <- rcp85 - present

BC_rcp45 <- Delta_rcp45 + Obs
BC_rcp85 <- Delta_rcp85 + Obs

saveRDS(Delta_rcp45, file = paste0(here::here('data/Output/Data-files//'),'DecTmeanDelta_rcp45'))
saveRDS(Delta_rcp85, file = paste0(here::here('data/Output/Data-files//'),'DecTmeanDelta_rcp85'))
saveRDS(BC_rcp45, file = paste0(here::here('data/Output/Data-files//'),'DecTmeanBC_rcp45'))
saveRDS(BC_rcp85, file = paste0(here::here('data/Output/Data-files//'),'DecTmeanBC_rcp85'))


## Wet/Dry season delta maps
#Wet: Oct-Apr
#Dry: May-Sep
data.dir<-"C:/Users/arunyon/3D Objects/Local-files/Git-repos/HI_Climate_Data/data/Output/Data-files/"

## Wet Precip rcp45
Nov<- readRDS(paste0(data.dir,"NovPrecipDelta_rcp45"))
Dec<-readRDS(paste0(data.dir,"DecPrecipDelta_rcp45"))
Jan<-readRDS(paste0(data.dir,"JanPrecipDelta_rcp45"))
Feb<-readRDS(paste0(data.dir,"FebPrecipDelta_rcp45"))
Mar<-readRDS(paste0(data.dir,"MarPrecipDelta_rcp45"))
Apr<-readRDS(paste0(data.dir,"AprPrecipDelta_rcp45"))

wet.pr<-(Nov+Dec+Jan+Feb+Mar+Apr)/6 #Avg Wet season delta
saveRDS(wet.pr, file = paste0(here::here('data/Output/Data-files//'),'WetPrecipDelta_rcp45'))

## Dry Precip rcp45
May<-readRDS(paste0(data.dir,"MayPrecipDelta_rcp45")) 
Jun<- readRDS(paste0(data.dir,"JunPrecipDelta_rcp45"))
Jul<-readRDS(paste0(data.dir,"JulPrecipDelta_rcp45"))
Aug<-readRDS(paste0(data.dir,"AugPrecipDelta_rcp45"))
Sep<-readRDS(paste0(data.dir,"SepPrecipDelta_rcp45"))
Oct<-readRDS(paste0(data.dir,"OctPrecipDelta_rcp45")) 

dry.pr<-(May+Jun+Jul+Aug+Sep+Oct)/6 #Avg Wet season delta
saveRDS(dry.pr, file = paste0(here::here('data/Output/Data-files//'),'DryPrecipDelta_rcp45'))

## Wet Precip rcp85
Nov<- readRDS(paste0(data.dir,"NovPrecipDelta_rcp85"))
Dec<-readRDS(paste0(data.dir,"DecPrecipDelta_rcp85"))
Jan<-readRDS(paste0(data.dir,"JanPrecipDelta_rcp85"))
Feb<-readRDS(paste0(data.dir,"FebPrecipDelta_rcp85"))
Mar<-readRDS(paste0(data.dir,"MarPrecipDelta_rcp85"))
Apr<-readRDS(paste0(data.dir,"AprPrecipDelta_rcp85"))

wet.pr<-(Nov+Dec+Jan+Feb+Mar+Apr)/6 #Avg Wet season delta
saveRDS(wet.pr, file = paste0(here::here('data/Output/Data-files//'),'WetPrecipDelta_rcp85'))

## Dry Precip rcp85
May<-readRDS(paste0(data.dir,"MayPrecipDelta_rcp85")) 
Jun<- readRDS(paste0(data.dir,"JunPrecipDelta_rcp85"))
Jul<-readRDS(paste0(data.dir,"JulPrecipDelta_rcp85"))
Aug<-readRDS(paste0(data.dir,"AugPrecipDelta_rcp85"))
Sep<-readRDS(paste0(data.dir,"SepPrecipDelta_rcp85"))
Oct<-readRDS(paste0(data.dir,"OctPrecipDelta_rcp85")) 

dry.pr<-(May+Jun+Jul+Aug+Sep+Oct)/6 #Avg Wet season delta
saveRDS(dry.pr, file = paste0(here::here('data/Output/Data-files//'),'DryPrecipDelta_rcp85'))

## Wet Tmean rcp45
Nov<- readRDS(paste0(data.dir,"NovTmeanDelta_rcp45"))
Dec<-readRDS(paste0(data.dir,"DecTmeanDelta_rcp45"))
Jan<-readRDS(paste0(data.dir,"JanTmeanDelta_rcp45"))
Feb<-readRDS(paste0(data.dir,"FebTmeanDelta_rcp45"))
Mar<-readRDS(paste0(data.dir,"MarTmeanDelta_rcp45"))
Apr<-readRDS(paste0(data.dir,"AprTmeanDelta_rcp45"))

wet.tmean<-(Nov+Dec+Jan+Feb+Mar+Apr)/6 #Avg Wet season delta
saveRDS(wet.tmean, file = paste0(here::here('data/Output/Data-files//'),'WetTmeanDelta_rcp45'))

## Dry Tmean rcp45
May<-readRDS(paste0(data.dir,"MayTmeanDelta_rcp45")) 
Jun<- readRDS(paste0(data.dir,"JunTmeanDelta_rcp45"))
Jul<-readRDS(paste0(data.dir,"JulTmeanDelta_rcp45"))
Aug<-readRDS(paste0(data.dir,"AugTmeanDelta_rcp45"))
Sep<-readRDS(paste0(data.dir,"SepTmeanDelta_rcp45"))
Oct<-readRDS(paste0(data.dir,"OctTmeanDelta_rcp45")) 

dry.tmean<-(May+Jun+Jul+Aug+Sep+Oct)/6 #Avg Wet season delta
saveRDS(dry.tmean, file = paste0(here::here('data/Output/Data-files//'),'DryTmeanDelta_rcp45'))

## Wet Tmean rcp85
Nov<- readRDS(paste0(data.dir,"NovTmeanDelta_rcp85"))
Dec<-readRDS(paste0(data.dir,"DecTmeanDelta_rcp85"))
Jan<-readRDS(paste0(data.dir,"JanTmeanDelta_rcp85"))
Feb<-readRDS(paste0(data.dir,"FebTmeanDelta_rcp85"))
Mar<-readRDS(paste0(data.dir,"MarTmeanDelta_rcp85"))
Apr<-readRDS(paste0(data.dir,"AprTmeanDelta_rcp85"))

wet.tmean<-(Nov+Dec+Jan+Feb+Mar+Apr)/6 #Avg Wet season delta
saveRDS(wet.tmean, file = paste0(here::here('data/Output/Data-files//'),'WetTmeanDelta_rcp85'))

## Dry Tmean rcp85
May<-readRDS(paste0(data.dir,"MayTmeanDelta_rcp85")) 
Jun<- readRDS(paste0(data.dir,"JunTmeanDelta_rcp85"))
Jul<-readRDS(paste0(data.dir,"JulTmeanDelta_rcp85"))
Aug<-readRDS(paste0(data.dir,"AugTmeanDelta_rcp85"))
Sep<-readRDS(paste0(data.dir,"SepTmeanDelta_rcp85"))
Oct<-readRDS(paste0(data.dir,"OctTmeanDelta_rcp85")) 

dry.tmean<-(May+Jun+Jul+Aug+Sep+Oct)/6 #Avg Wet season delta
saveRDS(dry.tmean, file = paste0(here::here('data/Output/Data-files//'),'DryTmeanDelta_rcp85'))


