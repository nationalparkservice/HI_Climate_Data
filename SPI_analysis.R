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

data.dir <- "C:/Users/achildress/Documents/Git-repos/HI_Climate_Data/data/Output/Data-files/"
plot.dir <- "C:/Users/achildress/Documents/Git-repos/HI_Climate_Data/data/Output/Plots/"

Monthly_dir_RF<-"D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/DynDS_6MonthlyTIF_State_250m/"
scens<-c("present","rcp45","rcp85")

tif <- read_stars("D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/DynDS_6MonthlyTIF_State_250m/rcp85/State_RF_rcp85_2080_01_250m.tif",along=band)
boundary <- st_read('./data/HALE/HALE_boundary.shp')
boundary <- st_transform(boundary, st_crs(tif))


zones <- st_read('./data/HALE/HALE_Ecoregions_Split.shp')
zones <- st_transform(zones, st_crs(tif))
climate_zones <- zones$Short_Name
rm(tif)

### SPI calcululations *Longman et al. cannot calculate PET in HI with these data
  file.list = list.files(path = paste0(Monthly_dir_RF,scens[1],"/"), pattern = '.tif', full.names = TRUE)
  RF_stack<-stack(file.list) # Read in as stack - cldn't fig out how to do calc on stars obj.
  cropstack<- crop(RF_stack, boundary) #crop to boundary
  cropstack<-stack(cropstack)
  
  dates<-as.Date(paste0(gsub("_","-",str_extract(file.list,"\\d{4}[_]\\d{2}")),"-01"),format="%Y-%m-%d")
  cropstack <- setZ(cropstack, dates)
  names(cropstack) <- as.yearmon(getZ(cropstack)) #zoo package
  
  r.mat <- as.matrix(cropstack)
  
  # Run spei()
  funSPImat <- function(x, sc, start, end, na.rm=TRUE,...) {
    dat <- ts(x, start = c(1990, 1), end = c(2009, 12), frequency = 12) #need to set dates of object. Could automate but not worth it now
    as.numeric((spi(dat, sc, ref.start = start, ref.end = end, na.rm=na.rm, ...))$fitted) 
  }
  
  fitted.mat <- t(apply(r.mat, 1, funSPImat, sc = 6, start = c(1990, 1), #dates of object again
                        end = c(2009, 12)))
  
  # Convert back to raster brick
  spi <- setValues(cropstack, fitted.mat)
  dates <- as.yearmon(getZ(cropstack))
  names(spi) <- as.yearmon(dates)
  
  spi.stars<-st_as_stars(spi) #convert back to stars
  spi.stars1<-st_crop(spi.stars,boundary) #crop to boundary
  spi.stars1 <- setNames(spi.stars1,"SPI-6")
  saveRDS(spi.stars1, file = paste0(data.dir,'SPI.monthly-',scens[1]))
  
  rm(spi.stars,spi.stars1,spi,cropstack,RF_stack)


for (i in 2:length(scens)){
  file.list = list.files(path = paste0(Monthly_dir_RF,scens[i],"/"), pattern = '.tif', full.names = TRUE)
  RF_stack<-stack(file.list) # Read in as stack - cldn't fig out how to do calc on stars obj.
  cropstack<- crop(RF_stack, boundary) #crop to boundary
  cropstack<-stack(cropstack)
  
  dates<-as.Date(paste0(gsub("_","-",str_extract(file.list,"\\d{4}[_]\\d{2}")),"-01"),format="%Y-%m-%d")
  cropstack <- setZ(cropstack, dates)
  names(cropstack) <- as.yearmon(getZ(cropstack)) #zoo package
  
  r.mat <- as.matrix(cropstack)
  
  # Run spei()
  funSPImat <- function(x, sc, start, end, na.rm=TRUE,...) {
    dat <- ts(x, start = c(2080, 1), end = c(2099, 12), frequency = 12) #need to set dates of object. Could automate but not worth it now
    as.numeric((spi(dat, sc, ref.start = start, ref.end = end, na.rm=na.rm, ...))$fitted) 
  }
  
  fitted.mat <- t(apply(r.mat, 1, funSPImat, sc = 6, start = c(2080, 1), #dates of object again
                        end = c(2099, 12)))
  
  # Convert back to raster brick
  spi <- setValues(cropstack, fitted.mat)
  dates <- as.yearmon(getZ(cropstack))
  names(spi) <- as.yearmon(dates)
  
  spi.stars<-st_as_stars(spi) #convert back to stars
  spi.stars1<-st_crop(spi.stars,boundary) #crop to boundary
  spi.stars1 <- setNames(spi.stars1,"SPI-6")
  saveRDS(spi.stars1, file = paste0(data.dir,'SPI.monthly-',scens[i]))
  rm(spi.stars,spi.stars1,spi,cropstack,RF_stack)
}
  
  
present.spi <- readRDS(paste0(data.dir,'SPI.monthly-',scens[1])) 
present.spi <- setNames(present.spi,"SPI.6")
present.avg <-st_apply(present.spi, c("x", "y"), mean,na.rm=TRUE)

spi45 <- readRDS(paste0(data.dir,'SPI.monthly-',scens[2])) 
spi45 <- setNames(spi45,"SPI.6")
spi45.avg <-st_apply(spi45, c("x", "y"), mean,na.rm=TRUE)

spi85 <- readRDS(paste0(data.dir,'SPI.monthly-',scens[3])) 
spi85 <- setNames(spi85,"SPI.6")
spi85.avg <-st_apply(spi85, c("x", "y"), mean,na.rm=TRUE)

df=data.frame()
for (i in 1:length(climate_zones)){
zone1 <- zones %>% filter(Short_Name == climate_zones[i])

present.spi.1 <- st_crop(present.spi, zone1)
present.spi.time <- st_apply((present.spi.1 %>% dplyr::select(SPI.6)), c("band"),mean,na.rm=TRUE, rename=FALSE)
df1 <- data.frame(present.spi.time)
df1$scen=scens[3]
df1$zone=climate_zones[i]
df <- rbind(df,df1)

spi45.1 <- st_crop(spi45, zone1)
spi45.time <- st_apply((spi45.1 %>% dplyr::select(SPI.6)), c("band"),mean,na.rm=TRUE, rename=FALSE)
df1 <- data.frame(spi45.time)
df1$scen=scens[3]
df1$zone=climate_zones[i]
df <- rbind(df,df1)

spi85.1 <- st_crop(spi85, zone1)
spi85.time <- st_apply((spi85.1 %>% dplyr::select(SPI.6)), c("band"),mean,na.rm=TRUE, rename=FALSE)
df1 <- data.frame(spi85.time)
df1$scen=scens[3]
df1$zone=climate_zones[i]
df <- rbind(df,df1)

}

PlotTheme = theme(axis.text=element_text(size=14),    #Text size for axis tick mark labels
                  axis.title.x=element_blank(),               #Text size and alignment for x-axis label
                  axis.title.y=element_text(size=18, vjust=0.5,  margin=margin(t=20, r=20, b=20, l=20)),              #Text size and alignment for y-axis label
                  plot.title=element_text(size=20,face="bold",hjust=0.5, margin=margin(t=20, r=20, b=20, l=20)),      #Text size and alignment for plot title
                  legend.title=element_text(size=18),                                                                    #Text size of legend category labels
                  legend.text=element_text(size=16),                                                                   #Text size of legend title
                  legend.position = "bottom")  

SPEI_annual_bar <- function(data,title,CFmethod=""){
  ggplot(data = data, aes(x=as.numeric(as.character(Year)), y=SPI.6,fill = col)) + 
    geom_bar(stat="identity",aes(fill=col),col="black") + 
    geom_hline(yintercept=-.5,linetype=2,colour="black",size=1) +
    scale_fill_manual(name="",values =c("turquoise2","orange1"),drop=FALSE) +
    labs(title = title, 
         x = "Date", y = "SPI") +
    guides(color=guide_legend(override.aes = list(size=7))) + PlotTheme +
    facet_wrap(~zone, nrow = 2, ncol=2,scales = "free_x")
}
df$Year<-as.numeric(str_sub(df$band,-4,-1))
df1 <- aggregate(SPI.6~Year+scen+zone,df,mean)
df1$col[df1$SPI.6>=0]<-"above average"
df1$col[df1$SPI.6<0]<-"below average"
df1$col<-factor(df1$col, levels=c("above average","below average"))

# Present
CFP<-subset(df1, scen == scens[1] )
SPEI_annual_bar(CFP,title="SPEI values for Present") 
ggsave("SPI-present-Annual-bar.png", path = plot.dir, width = 15, height = 9)

# CF1
CF1<-subset(df1, scen == scens[2] )
SPEI_annual_bar(CF1,title="SPEI values for Climate Future 1") 
ggsave("SPI-CF1-Annual-bar.png", path = plot.dir, width = 15, height = 9)

# CF2
CF2<-subset(df1, scen == scens[3] )
SPEI_annual_bar(CF2,title="SPEI values for Climate Future 2") 
ggsave("SPI-CF2-Annual-bar.png", path = plot.dir, width = PlotWidth, height = PlotHeight)


############ Drought characteristics ############
