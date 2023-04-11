library(tidyr)
library(dplyr)
library(stars)
rm(list=ls())

scens<-c("present","rcp45","rcp85")
RFpresent <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'RF.monthly-',scens[1]))
RFRCP45 <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'RF.monthly-',scens[2]))
RFRCP85 <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'RF.monthly-',scens[3]))
Tmeanpresent <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'Tmean.monthly-',scens[1]))
TmeanRCP45 <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'Tmean.monthly-',scens[2]))
TmeanRCP85 <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'Tmean.monthly-',scens[3]))

# For each climate zone extract monthly ts
zones <- st_read('./data/HALE/HALE_Ecoregions_Split.shp')
zones <- st_transform(zones, st_crs(RFpresent))
climate_zones <- zones$Short_Name

pcpdf <- data.frame()
tmndf <- data.frame()
for (i in 1:length(scens)){
    RF <- readRDS(file = paste0(here::here('data/Output/Data-files//'),'RF.monthly-',scens[i])) %>% 
      mutate(PrecipIn = Rainfall_mm/25.4) #%>% select(PrecipIn) 
    Tmean<- readRDS(file = paste0(here::here('data/Output/Data-files//'),'Tmean.monthly-',scens[i]))
   
     for(j in 1:length(climate_zones)){
      zones1 <-  zones %>% filter(Short_Name == climate_zones[j]) #subset zone
      RF1 <- st_crop(RF,zones1) #crop RF and Tmean to zone
      Tmean1 <- st_crop(Tmean,zones1)
      
      #summarize RF and Tmean by time
      pcp.time <- st_apply((RF1 %>% dplyr::select(PrecipIn)), c("date"),mean,na.rm=TRUE, rename=FALSE)
      pcpdf1 <- data.frame(pcp.time) %>% mutate(scen=scens[i],zone=climate_zones[j])
      pcpdf<-rbind(pcpdf,pcpdf1)
      
      tmean.time <- st_apply((Tmean1 %>% dplyr::select(TmeanF)), c("date"),mean,na.rm=TRUE, rename=FALSE)
      tmndf1 <-data.frame(tmean.time) %>% mutate(scen=scens[i],zone=climate_zones[j]) 
      tmndf<-rbind(tmndf,tmndf1)
  rm(pcpdf1,tmndf1)
  }
}

df <- merge(pcpdf,tmndf,by=c("date","scen","zone"))
df <- df[with(df, order(zone, scen,date)),]

write.csv(df,paste0(here::here('data/Output/Data-files//',"Zone-Monthly-data.csv")),row.names=F)



