library(tidyr)
library(dplyr)
library(ncdf4)
library(stars)
library(raster)
library(here)

DataDir <- "D:/HI_Data/Xue/"

Xue <- read_ncdf(paste0(DataDir,"hist_mean_rain.nc"))

nc<-nc_open(paste0(DataDir,"hist_mean_rain.nc"))
nc$dim

