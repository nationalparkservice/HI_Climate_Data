library(tidyr)
library(dplyr)
library(ncdf4)
library(stars)
library(raster)
library(here)
library(ggplot2)
library(ggthemes)
library(stringr)
rm(list=ls())

# ggplot() + # Resolution is course
#   # geom_raster(data = topo_df ,aes(x = x, y = y,alpha=HYP_HR_SR_W_1), show.legend=FALSE) +
#   geom_stars(data = Zhang[1,], alpha = 0.8) +
#   # facet_wrap("time") +
#   # scale_fill_viridis() + 
#   #coord_equal() + 
#   geom_sf(data = boundary, aes(), fill = NA) + 
#   theme_map() +
#   theme(legend.position = "bottom",
#         legend.key.width = unit(6, "cm"),
#         legend.key.height = unit(.3, "cm"),
#         legend.justification = "center",
#         plot.title=element_text(size=12,face="bold",hjust=0.5)) +
#   # plot.background = element_rect(colour = col, fill=NA, size=5)) +
#   # labs(fill = "Water Balance") +
#   scale_colour_manual(values = c(rgb(207, 31, 46, maxColorValue = 255)), "#ffda85")

Monthly_dir_RF<-"D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/DynDS_6MonthlyTIF_State_250m/"
scens<-c("present","rcp45","rcp85")

#Create list of file in directory
file.list = list.files(path = paste0(Monthly_dir_RF,scens[1],"/"), pattern = '.tif', full.names = TRUE)
file.list[1]
monthly_RF<-read_stars(c(file.list[1],file.list[2]),along="band")

monthly_RF<-setNames(monthly_RF,"Rainfall_mm")

date<-as.Date(paste0(gsub("_","-",str_extract(file.list[1:2],"\\d{4}[_]\\d{2}")),"-01"),format="%Y-%m-%d")
mystar_RF = st_set_dimensions(monthly_RF, 3, values =date, names = "date")

boundary <- st_read('./data/HALE/HALE_boundary.shp')
boundary <- st_transform(boundary, st_crs(mystar_RF))
RF_crop <-st_crop(mystar_RF,boundary)

RF_avg<-st_apply(RF_crop, c("x", "y"), mean) #mean of all time periods
RF_max <- st_apply(RF_crop, c("x", "y"), max)
RF_range <- st_apply(RF_crop, c("x", "y"), range) #min and max bands for each pixel

pcp.time <- st_apply((RF_crop %>% dplyr::select(Rainfall_mm)), c("date"),mean,na.rm=TRUE, rename=FALSE)

df <- data.frame(pcp.time)


Monthly_dir_Tmean<-"D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_TemperatureTIF_AllFolders/DynDS_6MonthlyTIF_State_250m/"

#Create list of file in directory
file.list = list.files(path = paste0(Monthly_dir_Tmean,scens[1],"/"), pattern = '.tif', full.names = TRUE)
file.list<-file.list[1:2]
monthly_Tmean<-read_stars(file.list,along="band")

monthly_Tmean<-setNames(monthly_Tmean,"TmeanK")

date<-as.Date(paste0(gsub("_","-",str_extract(file.list,"\\d{4}[_]\\d{2}")),"-01"),format="%Y-%m-%d")
mystar_Tmean = st_set_dimensions(monthly_Tmean, 3, values =date, names = "date")

Tmean_crop <-st_crop(mystar_Tmean,boundary)

Tmean_crop %>% mutate(TmeanF = (9/5*(TmeanK - 273.15)+32)) ->T


Tmean_avg<-st_apply(Tmean_crop, c("x", "y"), mean) #mean of all time periods

tmean.time <- st_apply((Tmean_crop %>% dplyr::select(TmeanF)), c("date"),mean,na.rm=TRUE, rename=FALSE)

df <- data.frame(tmean.time)

### SPI calcululations *Longman et al. cannot calculate PET in HI with these data
library(SPEI)
library(zoo)
file.list = list.files(path = paste0(Monthly_dir_RF,scens[1],"/"), pattern = '.tif', full.names = TRUE)
RF_stack<-stack(file.list[1:12]) # Read in as stack - cldn't fig out how to do calc on stars obj.
cropstack<- crop(RF_stack, boundary) #crop to boundary
cropstack<-stack(cropstack)

dates<-as.Date(paste0(gsub("_","-",str_extract(file.list[1:12],"\\d{4}[_]\\d{2}")),"-01"),format="%Y-%m-%d")
cropstack <- setZ(cropstack, dates)
names(cropstack) <- as.yearmon(getZ(cropstack)) #zoo package

r.mat <- as.matrix(cropstack)

# Run spei()
funSPImat <- function(x, sc, start, end, na.rm=TRUE,...) {
  dat <- ts(x, start = c(1990, 1), end = c(1990, 12), frequency = 1) #need to set dates of object. Could automate but not worth it now
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
