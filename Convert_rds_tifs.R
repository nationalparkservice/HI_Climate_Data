#Write .rds files to .gdb for park use
library(terra)
library(dplyr)
library(stars)
rm(list=ls())

DataDir <- "C:/Users/arunyon/DOI/HALE Climate Change Scenario Planning - 05.0 Climate futures development/Output/Data-files/"
OutDir <- "C:/Users/arunyon/3D Objects/Local-files/HALE-CFs"

file.list <- setdiff(list.files(path=DataDir),list.files(path=DataDir,pattern = ".csv")) #excludes files with .csv extension
file.list <- setdiff(file.list,list.files(path=DataDir,pattern = "monthly")) #get rid of monthly, multiband objs
#Tmean is in F, precip is in mm, SPI is an index
temp.list <- setdiff(file.list,list.files(path=DataDir,pattern = "Precip"))
precip.list <- setdiff(file.list,list.files(path=DataDir,pattern = "Tmean"))
spi.list <- list.files(path=DataDir,pattern = "SPI")

#write out temp files as is
for(i in 1:length(temp.list)){
rds <- readRDS(paste0(DataDir,temp.list[i]))
r<-rast(rds) 
writeRaster(r,file.path(OutDir,paste0(temp.list[i],".tif")))
}

#Convert to inches and write precip files
for(i in 1:length(precip.list)){
  rds <- readRDS(paste0(DataDir,precip.list[i])) #|>   #change var name
  rds<- rds |> mutate(pcp.in = !!rlang::sym(names(rds))/25.4) |>  select(pcp.in) #Some var are named mean and otehrs are Rainfall
  r<-rast(rds) 
  writeRaster(r,file.path(OutDir,paste0(precip.list[i],".tif")),overwrite=TRUE)
}
 
#SPI - need to read in monthly and create averages
# Average rasters
scens<-c("present","rcp45","rcp85")

present.spi <- readRDS(paste0(DataDir,'SPI.monthly-',scens[1]))
present.spi <- setNames(present.spi,"SPI.6")
present.avg <-st_apply(present.spi, c("x", "y"), mean,na.rm=TRUE,rename=FALSE)
r <- rast(present.avg)
writeRaster(r,file.path(OutDir,"SPI6-Obs.tif"),overwrite=TRUE)

spi45 <- readRDS(paste0(DataDir,'SPI.monthly-',scens[2]))
spi45 <- setNames(spi45,"SPI.6")
spi45.avg <-st_apply(spi45, c("x", "y"), mean,na.rm=TRUE,rename=FALSE)
r <- rast(spi45.avg)
writeRaster(r,file.path(OutDir,"SPI6-rcp45.tif"),overwrite=TRUE)

spi85 <- readRDS(paste0(DataDir,'SPI.monthly-',scens[3]))
spi85 <- setNames(spi85,"SPI.6")
spi85.avg <-st_apply(spi85, c("x", "y"), mean,na.rm=TRUE,rename=FALSE)
r <- rast(spi85.avg)
writeRaster(r,file.path(OutDir,"SPI6-rcp85.tif"),overwrite=TRUE)

delta45 <- spi45.avg - present.avg
r <- rast(delta45);writeRaster(r,file.path(OutDir,"SPI6-rcp45-Delta.tif"),overwrite=TRUE)
delta85 <- spi85.avg - present.avg
r <- rast(delta85);writeRaster(r,file.path(OutDir,"SPI6-rcp85-Delta.tif"),overwrite=TRUE)

