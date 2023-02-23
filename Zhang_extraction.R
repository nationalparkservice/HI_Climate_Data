library(tidyr)
library(dplyr)
library(ncdf4)
library(stars)
library(raster)
library(here)

DataDir <- "D:/HI_Data/Zhang/"

test.tif <- raster("D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/DynDS_7MonthlyTIF_byIsland_250m/Maui/present/ma_rf_pres_1990_01_250m.tif")
plot(test.tif)
crs(test.tif)

boundary <- st_read('./data/HALE/HALE_boundary.shp')
boundary <- st_transform(boundary, st_crs(test.tif))

plot(test.tif)
plot(boundary[1])

t<-mask(test.tif,boundary)
plot(t)

Zhang <- read_stars("D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/DynDS_7MonthlyTIF_byIsland_250m/Maui/present/ma_rf_pres_1990_01_250m.tif")
Zhang
plot(Zhang)
