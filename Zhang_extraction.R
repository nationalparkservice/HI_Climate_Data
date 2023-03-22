library(tidyr)
library(dplyr)
library(ncdf4)
library(stars)
library(raster)
library(here)
library(ggplot2)
library(ggthemes)
library(stringr)

DataDir <- "D:/HI_Data/Zhang/"

prcp.dir <- "D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/"
tmean.dir <- "D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_TemperatureTIF_AllFolders/"

tif <- read_stars("D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/DynDS_7MonthlyTIF_byIsland_250m/Maui/present/ma_rf_pres_1990_01_250m.tif",
                  dims=c("x","y","time"))
boundary <- st_read('./data/HALE/HALE_boundary.shp')
boundary <- st_transform(boundary, st_crs(tif))
Zhang <-st_crop(tif,boundary)

# Plots

# topo <- stack("C:/Users/achildress/Downloads/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W.tif") # read in as stack so can see RBG layers
# topo <- crop(topo, extent(boundary))
# # plotRGB(ak)
# topo_df  <- as.data.frame(topo, xy = TRUE) # this step is important to get it to plot in ggplot



ggplot() + # Resolution is course
  # geom_raster(data = topo_df ,aes(x = x, y = y,alpha=HYP_HR_SR_W_1), show.legend=FALSE) +
  geom_stars(data = Zhang[1,], alpha = 0.8) +
  # facet_wrap("time") +
  # scale_fill_viridis() + 
  #coord_equal() + 
  geom_sf(data = boundary, aes(), fill = NA) + 
  theme_map() +
  theme(legend.position = "bottom",
        legend.key.width = unit(6, "cm"),
        legend.key.height = unit(.3, "cm"),
        legend.justification = "center",
        plot.title=element_text(size=12,face="bold",hjust=0.5)) +
  # plot.background = element_rect(colour = col, fill=NA, size=5)) +
  # labs(fill = "Water Balance") +
  scale_colour_manual(values = c(rgb(207, 31, 46, maxColorValue = 255)), "#ffda85")

Monthly_dir<-"D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/DynDS_6MonthlyTIF_State_250m/rcp45/"

#Create list of file in directory
file.list = list.files(path = Monthly_dir, pattern = '.tif', full.names = TRUE)
file.list[1]
z<-read_stars(c(file.list[1],file.list[2]),along="band")

#crop list
l <- list() # Create a list to put the stars objects into

for(i in 1:length(hist_filelist)){
  suppressMessages(
    l[[i]] <- read_ncdf(hist_filelist[i], curvilinear = c("longitude", "latitude")) # need to read in as ncdf or coordinate system does not translate (not sure why)
  )
}

# Crop

cropped_hist <- list() # create list for cropped stars objects

for(i in 1:length(l)){ # add cropped stars objects to a new list
  nc = l[[i]]
  nc = st_transform(nc, st_crs(shp))
  nc_crop = nc[shp]
  cropped_hist[[i]] = nc_crop
}

cropped_st_hist <- list()

for(i in 1:length(cropped_hist)){
  cropped_st_hist[[i]] <- st_as_stars(cropped_hist[[i]])
}
#stack with year-month converted as time dimension
#Bias correction









path = "C:/Users/achildress/Documents/HI_test/Zhang_monthly/present/"
file.list = list.files(path = path, pattern = '.tif', full.names = TRUE)
z<-read_stars(file.list,along="band")
z<-setNames(z,"Rainfall")

date<-as.Date(paste0(gsub("_","-",str_extract(file.list,"\\d{4}[_]\\d{2}")),"-01"),format="%Y-%m-%d")
mystar = st_set_dimensions(z, 3, values =date, names = "date")

boundary <- st_read('./data/HALE/HALE_boundary.shp')
boundary <- st_transform(boundary, st_crs(mystar))
Zhang <-st_crop(mystar,boundary)


st_apply(Zhang, c("x", "y"), mean) #mean of all time periods

pcp.time <- st_apply((Zhang %>% dplyr::select(Rainfall)), c("date"),mean,na.rm=TRUE, rename=FALSE)

df <- data.frame(pcp.time)

