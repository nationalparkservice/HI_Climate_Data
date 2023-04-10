library(tidyr)
library(dplyr)
library(ncdf4)
library(stars)
library(raster)
library(here)
library(ggplot2)
library(ggthemes)
library(stringr)
library(zoo)
library(SPEI)
library(viridis)
rm(list=ls())

Monthly_dir_RF<-"D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/DynDS_6MonthlyTIF_State_250m/"
scens<-c("present","rcp45","rcp85")

### SPI calcululations *Longman et al. cannot calculate PET in HI with these data
for (i in 1:length(scens)){
  file.list = list.files(path = paste0(Monthly_dir_RF,scens[1],"/"), pattern = '.tif', full.names = TRUE)
  RF_stack<-stack(file.list[1:12]) # Read in as stack - cldn't fig out how to do calc on stars obj.
  cropstack<- crop(RF_stack, boundary) #crop to boundary
  cropstack<-stack(cropstack)
  
  dates<-as.Date(paste0(gsub("_","-",str_extract(file.list,"\\d{4}[_]\\d{2}")),"-01"),format="%Y-%m-%d")
  cropstack <- setZ(cropstack, dates)
  names(cropstack) <- as.yearmon(getZ(cropstack)) #zoo package
  
  r.mat <- as.matrix(cropstack)
  
  # Run spei()
  funSPImat <- function(x, sc, start, end, na.rm=TRUE,...) {
    dat <- ts(x, start = c(1990, 1), end = c(1990, 12), frequency = 12) #need to set dates of object. Could automate but not worth it now
    as.numeric((spi(dat, sc, ref.start = start, ref.end = end, na.rm=na.rm, ...))$fitted) 
  }
  
  fitted.mat <- t(apply(r.mat, 1, funSPImat, sc = 1, start = c(1990, 1), #dates of object again
                        end = c(1990, 12)))
  
  # Convert back to raster brick
  spi <- setValues(cropstack, fitted.mat)
  dates <- as.yearmon(getZ(cropstack))
  names(spi) <- as.yearmon(dates)
  
  spi.stars<-st_as_stars(spi) #convert back to stars
  spi.stars1<-st_crop(spi.stars,boundary) #crop to boundary
  saveRDS(spi.stars1, file = paste0(here::here('data/Output/Data-files//'),'Tmean.monthly-',scens[i]))
  rm(spi.stars,spi.stars1,spi,cropstack,RF_stack)
}
dates[1]
