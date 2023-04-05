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

tif <- read_stars("D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/DynDS_6MonthlyTIF_State_250m/rcp85/State_RF_rcp85_2080_01_250m.tif",along=band)
boundary <- st_read('./data/HALE/HALE_boundary.shp')
boundary <- st_transform(boundary, st_crs(tif))
rm(tif)

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

## Monthly Rainfall
Monthly_dir_RF<-"D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/DynDS_6MonthlyTIF_State_250m/"
scens<-c("present","rcp45","rcp85")
df <- data.frame()
#Create list of file in directory
for(i in 1:length(scens)){
file.list = list.files(path = paste0(Monthly_dir_RF,scens[i],"/"), pattern = '.tif', full.names = TRUE)
monthly_RF<-read_stars(file.list,along="band")

monthly_RF<-setNames(monthly_RF,"Rainfall_mm")

date<-as.Date(paste0(gsub("_","-",str_extract(file.list,"\\d{4}[_]\\d{2}")),"-01"),format="%Y-%m-%d")
mystar_RF = st_set_dimensions(monthly_RF, 3, values =date, names = "date")

RF_crop <-st_crop(mystar_RF,boundary)
saveRDS(RF_crop, file = paste0(here::here('data/Output/Data-files//'),'RF.monthly-',scens[i]))

# RF_avg<-st_apply(RF_crop, c("x", "y"), mean) #mean of all time periods
# RF_max <- st_apply(RF_crop, c("x", "y"), max)
# RF_range <- st_apply(RF_crop, c("x", "y"), range) #min and max bands for each pixel

pcp.time <- st_apply((RF_crop %>% dplyr::select(Rainfall_mm)), c("date"),mean,na.rm=TRUE, rename=FALSE)

df1 <- data.frame(pcp.time)
df1$scen=scens[i]
df<-rbind(df,df1)
rm(df1)
}
write.csv(df, file=paste0(here::here('data/Output/Data-files//'),'RF.monthly.HALE.csv'))


## Monthly Temperature
Monthly_dir_Tmean<-"D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_TemperatureTIF_AllFolders/DynDS_6MonthlyTIF_State_250m/"
df <- data.frame()

#Create list of file in directory
for (i in 1:length(scens)){
file.list = list.files(path = paste0(Monthly_dir_Tmean,scens[i],"/"), pattern = '.tif', full.names = TRUE)
monthly_Tmean<-read_stars(file.list,along="band")

monthly_Tmean<-setNames(monthly_Tmean,"TmeanK")

date<-as.Date(paste0(gsub("_","-",str_extract(file.list,"\\d{4}[_]\\d{2}")),"-01"),format="%Y-%m-%d")
mystar_Tmean = st_set_dimensions(monthly_Tmean, 3, values =date, names = "date")

Tmean_crop <-st_crop(mystar_Tmean,boundary)
Tmean_crop %>% mutate(TmeanF = (9/5*(TmeanK - 273.15)+32))  ->Tmean_crop
Tmean_crop["TmeanF"]->Tmean_crop

saveRDS(Tmean_crop, file = paste0(here::here('data/Output/Data-files//'),'Tmean.monthly-',scens[i]))
# Tmean_avg<-st_apply(Tmean_crop, c("x", "y"), mean) #mean of all time periods

tmean.time <- st_apply((Tmean_crop %>% dplyr::select(TmeanF)), c("date"),mean,na.rm=TRUE, rename=FALSE)
df1 <- data.frame(tmean.time)
df1$scen=scens[i]
df<-rbind(df,df1)
rm(df1)
}
write.csv(df, file=paste0(here::here('data/Output/Data-files//'),'Tmean.monthly.HALE.csv'))

rm(monthly_RF,mystar_RF,monthly_Tmean,mystar_Tmean,pcp.time,RF_avg,RF_crop,RF_max,RF_range,Tmean_crop,tmean.time)
gc()


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
