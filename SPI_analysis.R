library(tidyr)
library(plyr)
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
library(lemon);library(ggpubr);library(gridExtra);library(grid);library(gtable)
rm(list=ls())

data.dir <- "C:/Users/achildress/Documents/Git-repos/HI_Climate_Data/data/Output/Data-files/"
plot.dir <- "C:/Users/achildress/Documents/Git-repos/HI_Climate_Data/data/Output/Plots/"

Monthly_dir_RF<-"D:/HI_Data/Zhang/Processed/HI_Downscaling_Share/DynDS_RainfallTIF_AllFolders/DynDS_6MonthlyTIF_State_250m/"
scens<-c("present","rcp45","rcp85")

tif <- read_stars("./data/HALE/tif.tif",along=band)
boundary <- st_read('./data/HALE/HALE_boundary.shp')
boundary <- st_transform(boundary, st_crs(tif))


zones <- st_read('./data/HALE/HALE_Ecoregions_Split.shp')
zones <- st_transform(zones, st_crs(tif))
climate_zones <- zones$Short_Name
rm(tif)

cols <- c("#9A9EE5","#E10720")

# ### SPI calcululations *Longman et al. cannot calculate PET in HI with these data
# Commented out after .rds files saved
#   file.list = list.files(path = paste0(Monthly_dir_RF,scens[1],"/"), pattern = '.tif', full.names = TRUE)
#   RF_stack<-stack(file.list) # Read in as stack - cldn't fig out how to do calc on stars obj.
#   cropstack<- crop(RF_stack, boundary) #crop to boundary
#   cropstack<-stack(cropstack)
#   
#   dates<-as.Date(paste0(gsub("_","-",str_extract(file.list,"\\d{4}[_]\\d{2}")),"-01"),format="%Y-%m-%d")
#   cropstack <- setZ(cropstack, dates)
#   names(cropstack) <- as.yearmon(getZ(cropstack)) #zoo package
#   
#   r.mat <- as.matrix(cropstack)
#   
#   # Run spei()
#   funSPImat <- function(x, sc, start, end, na.rm=TRUE,...) {
#     dat <- ts(x, start = c(1990, 1), end = c(2009, 12), frequency = 12) #need to set dates of object. Could automate but not worth it now
#     as.numeric((spi(dat, sc, ref.start = start, ref.end = end, na.rm=na.rm, ...))$fitted) 
#   }
#   
#   fitted.mat <- t(apply(r.mat, 1, funSPImat, sc = 6, start = c(1990, 1), #dates of object again
#                         end = c(2009, 12)))
#   
#   # Convert back to raster brick
#   spi <- setValues(cropstack, fitted.mat)
#   dates <- as.yearmon(getZ(cropstack))
#   names(spi) <- as.yearmon(dates)
#   
#   spi.stars<-st_as_stars(spi) #convert back to stars
#   spi.stars1<-st_crop(spi.stars,boundary) #crop to boundary
#   spi.stars1 <- setNames(spi.stars1,"SPI-6")
#   saveRDS(spi.stars1, file = paste0(data.dir,'SPI.monthly-',scens[1]))
#   
#   h.cropstack <- cropstack
#   
#   rm(spi.stars,spi.stars1,spi,cropstack,RF_stack)
# 
# 
# for (i in 2:length(scens)){
#   file.list = list.files(path = paste0(Monthly_dir_RF,scens[i],"/"), pattern = '.tif', full.names = TRUE)
#   RF_stack<-stack(file.list) # Read in as stack - cldn't fig out how to do calc on stars obj.
#   cropstack<- crop(RF_stack, boundary) #crop to boundary
#   cropstack<-stack(h.cropstack,cropstack)
#   
#   h.dates <- as.Date(seq(as.Date("2060-01-01"), as.Date("2079-12-01"), by="months"))
#   dates<-as.Date(paste0(gsub("_","-",str_extract(file.list,"\\d{4}[_]\\d{2}")),"-01"),format="%Y-%m-%d")
#   dates<-c(h.dates,dates)
#   cropstack <- setZ(cropstack, dates)
#   names(cropstack) <- as.yearmon(getZ(cropstack)) #zoo package
#   
#   r.mat <- as.matrix(cropstack)
#   
#   # Run spei()
#   funSPImat <- function(x, sc, start, end, na.rm=TRUE,...) {
#     dat <- ts(x, start = c(2060, 1), end = c(2099, 12), frequency = 12) #need to set dates of object. Could automate but not worth it now
#     as.numeric((spi(dat, sc, ref.start = start, ref.end = c(2079,12), na.rm=na.rm, ...))$fitted) 
#   }
#   
#   fitted.mat <- t(apply(r.mat, 1, funSPImat, sc = 6, start = c(2060, 1), #dates of object again
#                         end = c(2099, 12)))
#   
#   # Convert back to raster brick
#   spi <- setValues(cropstack, fitted.mat)
#   dates <- as.yearmon(getZ(cropstack))
#   names(spi) <- as.yearmon(dates)
#   
#   spi.stars<-st_as_stars(spi) #convert back to stars
#   spi.stars1<-st_crop(spi.stars[,,,-c(1:240)],boundary) #crop to boundary
#   spi.stars1 <- setNames(spi.stars1,"SPI-6")
#   saveRDS(spi.stars1, file = paste0(data.dir,'SPI.monthly-',scens[i]))
#   rm(spi.stars,spi.stars1,spi,cropstack,RF_stack)
# }

######################
## Average maps
# 
# present.spi <- readRDS(paste0(data.dir,'SPI.monthly-',scens[1])) 
# present.spi <- setNames(present.spi,"SPI.6")
# present.avg <-st_apply(present.spi, c("x", "y"), mean,na.rm=TRUE)
# 
# spi45 <- readRDS(paste0(data.dir,'SPI.monthly-',scens[2])) 
# spi45 <- setNames(spi45,"SPI.6")
# spi45.avg <-st_apply(spi45, c("x", "y"), mean,na.rm=TRUE)
# 
# spi85 <- readRDS(paste0(data.dir,'SPI.monthly-',scens[3])) 
# spi85 <- setNames(spi85,"SPI.6")
# spi85.avg <-st_apply(spi85, c("x", "y"), mean,na.rm=TRUE)
# 
# topo <- stack('./data/HALE/HALENatEa1.tif')
# topo_df  <- as.data.frame(topo, xy = TRUE) 
# 
# scale.min = min(c(present.avg$mean, spi45.avg$mean, spi85.avg$mean),na.rm=TRUE)
# scale.max = max(c(present.avg$mean, spi45.avg$mean, spi85.avg$mean),na.rm=TRUE)
# 
# map.plot <- function(data, title,xaxis,metric,col){
#   ggplot() + 
#     geom_raster(data = topo_df ,aes(x = x, y = y,alpha=HALENatEa1_1), show.legend=FALSE) +
#     geom_stars(data = data, alpha = 0.8) + 
#     geom_sf(data = boundary, aes(), fill = NA) +
#     geom_sf(data = zones, aes(), fill = NA) +
#     scale_fill_viridis(direction=-1, option = "turbo", limits = c(scale.min, scale.max),begin=.45,  
#                        guide = guide_colorbar(title.position = "top", title.hjust = 0.5),oob = scales::squish) + #mako for WB delta
#     labs(title = title) +
#     theme_map() +
#     theme(legend.position = "bottom",
#           legend.key.width = unit(6, "cm"),
#           legend.key.height = unit(.3, "cm"),
#           legend.justification = "center",
#           plot.title=element_text(size=12,face="bold",hjust=0.5),
#           plot.background = element_rect(colour = col, fill=NA, linewidth = 5)) + 
#     labs(fill = metric)
# }
# cfp <- map.plot(data=present.avg,title="Recent",metric="Average SPI-6",col="grey")
# cf1 <- map.plot(data=spi45.avg,title="Climate Future 1",metric="Average SPI-6",col=cols[1])
# cf2 <- map.plot(data=spi85.avg,title="Climate Future 2",metric="Average SPI-6",col=cols[2])
# 
# maps <- grid_arrange_shared_legend(cfp,cf1, cf2,ncol = 3, nrow = 1, position = "bottom", 
#                                    top = textGrob(paste0("Average Standardized Precipitation Index"),
#                                                   gp=gpar(fontface="bold", col="black", fontsize=16)))
# ggsave("SPI-maps.png", plot=maps,width = 15, height = 9, path = plot.dir,bg="white")
# 
# df=data.frame()
# for (i in 1:length(climate_zones)){
# zone1 <- zones %>% filter(Short_Name == climate_zones[i])
# 
# present.spi.1 <- st_crop(present.spi, zone1)
# present.spi.time <- st_apply((present.spi.1 %>% dplyr::select(SPI.6)), c("band"),mean,na.rm=TRUE, rename=FALSE)
# df1 <- data.frame(present.spi.time)
# df1$scen=scens[1]
# df1$zone=climate_zones[i]
# df <- rbind(df,df1)
# 
# spi45.1 <- st_crop(spi45, zone1)
# spi45.time <- st_apply((spi45.1 %>% dplyr::select(SPI.6)), c("band"),mean,na.rm=TRUE, rename=FALSE)
# df1 <- data.frame(spi45.time)
# df1$scen=scens[2]
# df1$zone=climate_zones[i]
# df <- rbind(df,df1)
# 
# spi85.1 <- st_crop(spi85, zone1)
# spi85.time <- st_apply((spi85.1 %>% dplyr::select(SPI.6)), c("band"),mean,na.rm=TRUE, rename=FALSE)
# df1 <- data.frame(spi85.time)
# df1$scen=scens[3]
# df1$zone=climate_zones[i]
# df <- rbind(df,df1)
# 
# }
# 
# write.csv(df,paste0(plot.dir,"SPI-by-zone.csv"),row.names = FALSE)

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
    guides(color=guide_legend(override.aes = list(linewidth=7))) + PlotTheme +
    # theme(plot.background = element_rect(colour = col, fill=NA, linewidth = 5)) +
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
ggsave("SPI-CF2-Annual-bar.png", path = plot.dir, width = 15, height = 9)


############ Drought characteristics ############

df <- read.csv(paste0(plot.dir,"SPI-by-zone.csv"))
df$Year<-as.numeric(str_sub(df$band,-4,-1))

Drought_characteristics <- data.frame()
for (z in 1:length(climate_zones)){
df1 <- subset(df, zone=climate_zones[z])
df2<-aggregate(SPI.6~Year+scen,df1,mean)
# Split into periods
drt3<-subset(df2, scen==scens[1])
min(drt3$SPI.6)

Future.drt<-subset(df2, scen %in% scens[2:3])
min(Future.drt$SPI.6)

# Calculate drought characteristics
drt3$Drought=0
drt3$Drought[which(drt3$SPI.6 < -0.5)] <- 1

# Drought Duration calculation
# 1 Create var for beginnign drought and var for end drought, then count months between
head(drt3)

# Create count of years within CF
length(drt3$Year)
drt3$count<-seq(1, length(drt3$Year),1) 

drt3$length<-0
drt3$length <- drt3$Drought * unlist(lapply(rle(drt3$Drought)$lengths, seq_len))
mean(drt3$length[drt3$length>0])

# To get duration, now just remove those that are not droughts and do calculations on length

# Give each drought period an ID
D<-which(drt3$length==1)
HistoricalDrought<-data.frame()
HistoricalDrought<-setNames(data.frame(matrix(ncol=9,nrow=length(D))),c("DID","Start","End","Year","per","CF","duration","severity","peak"))
HistoricalDrought$Start = Sys.time(); HistoricalDrought$End = Sys.time()
HistoricalDrought$per<-as.factor("H")


# Calculate variables for each drought period
for (i in 1:length(D)){
  HistoricalDrought$DID[i]<-i
  HistoricalDrought$Start[i]<-as.Date(paste0(drt3$Year[D[i]],"-01-01"),format="%Y-%m-%d")
  HistoricalDrought$Year[i]<-drt3$Year[D[i]]
}

ND<- which((drt3$length == 0) * unlist(lapply(rle(drt3$length)$lengths, seq_len)) == 1)
if(ND[1]==1) ND<-ND[2:length(ND)]
if(drt3$Drought[length(drt3$Drought)]==1) ND[length(ND)+1]<-length(drt3$length)

###### !!!!!!!!!!! 
# If last row in drought df is a drought period - use next line of code. Otherwies proceed.
# ND[length(ND)+1]<-length(drt3$length) #had to add this step because last drought went until end of df so no end in ND

#Duration # months SPEI < truncation; Severity # Sum(SPEI) when SPEI < truncation; Peak # min(SPEI) when SPEI < truncation

for (i in 1:length(ND)){
  HistoricalDrought$End[i]<-as.Date(paste0(drt3$Year[ND[i]],"-01-01"),format="%Y-%m-%d")
  HistoricalDrought$duration[i]<-drt3$length[ND[i]-1]
  HistoricalDrought$severity[i]<-sum(drt3$SPI[D[i]:(ND[i]-1)])
  HistoricalDrought$peak[i]<-min(drt3$SPI[D[i]:(ND[i]-1)])
}

## Freq
d<-which(drt3$length==1)
nd<-which((drt3$length == 0) * unlist(lapply(rle(drt3$length)$lengths, seq_len)) == 1)
if(length(nd)>length(d)) {nd=nd[2:length(nd)]}
for (j in 1:length(d)){
  HistoricalDrought$freq[which(HistoricalDrought$Year==drt3$Year[d[j]])] <-
    drt3$count[d[j+1]]-drt3$count[nd[j]]
}

####### Future
# Calculate drought characteristics
Future.drt$Drought=0
Future.drt$Drought[which(Future.drt$SPI.6 < -.5)] <- 1
colnames(Future.drt)[colnames(Future.drt) == "scen"] ="CF"

# Drought Duration calculation
# 1 Create var for beginnign drought and var for end drought, then count months between
head(Future.drt)

# Create count of months within CF
length(Future.drt$CF)/length(unique(Future.drt$CF))


# Give each drought period an ID
FutureDrought<-data.frame()
FutureDrought<-setNames(data.frame(matrix(ncol=10,nrow=0)),c("DID","Start","End","Year","per","CF","duration","severity","peak","freq"))
FutureDrought.i <- FutureDrought

# Future.drt$CF <- droplevels(Future.drt$CF)
CF.split<-split(Future.drt,Future.drt$CF)

# Calculate drought characteristics for each CF -- have to split by CF to avoid mixing up counts
for (c in 1:length(CF.split)){
  name=as.character(unique(CF.split[[c]]$CF))
  
  
  CF.split[[c]]$count<-rep(seq(1, length(CF.split[[c]]$CF)/length(unique(CF.split[[c]]$CF)), 
                               1),length(unique(CF.split[[c]]$CF))) # repeat # of CFs 
  
  CF.split[[c]]$length<-0
  CF.split[[c]]$length <- CF.split[[c]]$Drought * unlist(lapply(rle(CF.split[[c]]$Drought)$lengths, seq_len))
  mean(CF.split[[c]]$length[CF.split[[c]]$length>0])
  
  D<-which(CF.split[[c]]$length==1)
  fd <- FutureDrought.i
  if(length(D)>0) {
    fd[length(D),] <- NA
  } else {
    fd[1,] <- 0
  }
  
  fd$per<-as.factor("F")
  fd$CF = name
  
  for (i in 1:length(D)){
    fd$DID[i]<-i
    fd$Start[i]<-as.character(as.Date(paste0(CF.split[[c]]$Year[D[i]],"-01-01"),format="%Y-%m-%d"))
    fd$Year[i]<-CF.split[[c]]$Year[D[i]]
  }
  
  ND<- which((CF.split[[c]]$length == 0) * unlist(lapply(rle(CF.split[[c]]$length)$lengths, seq_len)) == 1)
  if(ND[1]==1 & length(D)>0) ND<-ND[2:length(ND)]
  if(CF.split[[c]]$Drought[length(CF.split[[c]]$Drought)]==1) ND[length(ND)+1]<-length(CF.split[[c]]$length)
  
  for (i in 1:length(ND)){
    fd$End[i]<-as.character(as.Date(paste0(CF.split[[c]]$Year[ND[i]],"-01-01"),format="%Y-%m-%d"))
    if(length(D) > 0){
      fd$severity[i]<-sum(CF.split[[c]]$SPI[D[i]:(ND[i]-1)])
      fd$peak[i]<-min(CF.split[[c]]$SPI[D[i]:(ND[i]-1)])
      fd$duration[i]<-CF.split[[c]]$length[ND[i]-1]
    } else {
      fd$severity[i]<-0
      fd$peak[i]<-0
      fd$duration[i]<-0
    }
  }
  
  ## Freq
  d<-which(CF.split[[c]]$length==1)
  nd<-which((CF.split[[c]]$length == 0) * unlist(lapply(rle(CF.split[[c]]$length)$lengths, seq_len)) == 1)
  if(length(nd)>length(d)) {nd=nd[2:length(nd)]}
  for (j in 1:length(d)){
    fd$freq[which(fd$Year==CF.split[[c]]$Year[d[j]])] <-
      CF.split[[c]]$count[d[j+1]]-CF.split[[c]]$count[nd[j]]
  }
  FutureDrought <- rbind(FutureDrought,fd)
  rm(fd)
}
rm(FutureDrought.i)
FutureDrought$CF = factor(FutureDrought$CF, levels = scens[2:3])
Future.drt <- ldply(CF.split, data.frame) #convert back to df

########### Merge
HistoricalDrought$CF="Recent"
head(HistoricalDrought)
head(FutureDrought)
Drought<-rbind(HistoricalDrought,FutureDrought)
write.csv(Drought,paste0(plot.dir,"Drt.all.csv"),row.names=FALSE)  # csv with all drought events

Hist_char<-setNames(data.frame(matrix(ncol=6,nrow=1)),c("CF","per","Duration","Severity","Intensity","Frequency"))
Hist_char$CF<-"Recent"
Hist_char$per<-"R"
Hist_char$Frequency<-mean(HistoricalDrought$freq,na.rm=TRUE)
Hist_char$Duration<-mean(HistoricalDrought$duration)
Hist_char$Severity<-mean(HistoricalDrought$severity)
Hist_char$Intensity<-mean(HistoricalDrought$peak)
Hist_char$Drt.Free <- mean(rle(drt3$length)$lengths[which(rle(drt3$length)$values==0)])


Drought_char<-setNames(data.frame(matrix(ncol=6,nrow=length(levels(FutureDrought$CF)))),c("CF","per","Duration","Severity","Intensity","Frequency"))
Drought_char$CF<-levels(FutureDrought$CF)
Drought_char$per<-"F"
for (i in 1:length(Drought_char$CF)){
  name<-Drought_char$CF[i]
  Drought_char$Frequency[i]<-mean(FutureDrought$freq[which(FutureDrought$CF == name)],na.rm=TRUE)
  Drought_char$Duration[i]<-mean(FutureDrought$duration[which(FutureDrought$CF == name)])
  Drought_char$Severity[i]<-mean(FutureDrought$severity[which(FutureDrought$CF == name)])
  Drought_char$Intensity[i]<-mean(FutureDrought$peak[which(FutureDrought$CF == name)])
  Drought_char$Drt.Free[i]<-mean(rle(subset(Future.drt,CF==name)$length)$lengths[which(rle(subset(Future.drt,CF==name)$length)$values==0)])
}

Drought_char<-rbind(Hist_char,Drought_char) 
Drought_char$zone=climate_zones[z]
Drought_characteristics <-rbind(Drought_char,Drought_characteristics)
# rm(Drought_char)
}

# csv for averages for each CF for hist and future periods
write.csv(Drought_characteristics,paste0(plot.dir,"Drought_characteristics.csv"),row.names=FALSE)


########################################### BAR PLOTS ###############################################
BarPlotTheme = theme(axis.text.x=element_text(size=16),    #Text size for axis tick mark labels
                     axis.text.y=element_text(size=14),
                     axis.title.x=element_blank(),               #Text size and alignment for x-axis label
                     axis.title.y=element_text(size=18, vjust=0.5,  margin=margin(t=20, r=20, b=20, l=20)),              #Text size and alignment for y-axis label
                     plot.title=element_text(size=20,face="bold",hjust=0.5, margin=margin(t=20, r=20, b=20, l=20)),      #Text size and alignment for plot title
                     legend.position = "none") 

var_bar_plot <- function(data,var, cols, title, ylab){ 
  At<-aggregate(eval(parse(text=var))~CF,data=data,mean); 
  names(At)<-c("CF",var) 
  p<-ggplot(At, aes(x=CF,y=(eval(parse(text=var))),fill=CF)) + 
    geom_bar(stat="identity",position="dodge",colour="black") + 
    BarPlotTheme + 
    # coord_cartesian(ylim=c(0, 40)) + 
    labs(title = title,  
         y = ylab)  + 
    scale_fill_manual(name="",values = cols)  
  # annotate(geom="text", x=Inf, y=-Inf, label=CFmethod,color="black",vjust=-1,hjust=1)
  if(min(eval(parse(text=paste("At$",var,sep=""))))<20) {p = p + coord_cartesian(ylim = c(0, max(eval(parse(text=paste("At$",var,sep=""))))))}  
  else{p= p + coord_cartesian(ylim = c(min(eval(parse(text=paste("At$",var,sep=""))))*.9, max(eval(parse(text=paste("At$",var,sep=""))))))} 
  p 
} 


SPEI_annual_bar <- function(data,title,CFmethod=""){
  ggplot(data = data, aes(x=as.numeric(as.character(Year)), y=SPI.6,fill = col)) +
    geom_bar(stat="identity",aes(fill=col),col="black") +
    geom_hline(yintercept=-.5,linetype=2,colour="black",size=1) +
    scale_fill_manual(name="",values =c("turquoise2","orange1"),drop=FALSE) +
    labs(title = title,
         x = "Date", y = "SPI") +
    guides(color=guide_legend(override.aes = list(linewidth=7))) + PlotTheme 
  # theme(plot.background = element_rect(colour = col, fill=NA, linewidth = 5)) +
}

#Drought duration barplot
Drt.all = Drought_characteristics
Drt.all$CF = factor(Drt.all$CF, levels = c("Recent",scens[2:3]))

for(z in 1:length(climate_zones)){
Drought_all <- subset(Drt.all,zone == climate_zones[z])

#Drought duration barplot
var_bar_plot(Drought_all,"Duration", c("grey",cols), paste0(climate_zones[z], "-Average Drought Duration"), "Years")
ggsave(paste0(climate_zones[z],"-DroughtDuration-Bar.png"), path = plot.dir, height=6, width=9)

#Drought severity barplot
var_bar_plot(Drought_all,"Severity", c("grey",cols), paste0(climate_zones[z], "-Average Drought Severity"), 
             "Severity (Intensity * Duration)") + coord_cartesian(ylim = c(0, min(Drought_all$Severity)))
ggsave(paste0(climate_zones[z],"-DroughtSeverity-Bar.png"), path = plot.dir, height=6, width=9)

#Drought intensity barplot
var_bar_plot(Drought_all,"Intensity", c("grey",cols), paste0(climate_zones[z], "-Average Drought Intensity"), 
             "Intensity (Minimum SPEI values)") + coord_cartesian(ylim = c(0, min(Drought_all$Intensity)))
ggsave(paste0(climate_zones[z],"-DroughtIntensity-Bar.png"), path = plot.dir, height=6, width=9)

#Drought-free interval barplot
var_bar_plot(Drought_all,"Drt.Free", c("grey",cols), paste0(climate_zones[z], "-Average Drought-Free Interval"), 
             "Years")
ggsave(paste0(climate_zones[z],"-DroughtFrequency-Bar.png"), path = plot.dir, height=6, width=9)


####################################### REPORT FIGURES ##############################################
# Option 1
df <- read.csv(paste0(plot.dir,"SPI-by-zone.csv"))
df$Year<-as.numeric(str_sub(df$band,-4,-1))
df1 <- subset(df, zone==climate_zones[z])
df2<-aggregate(SPI.6~Year+scen+zone,df1,mean)
df2$col[df2$SPI.6>=0]<-"above average"
df2$col[df2$SPI.6<0]<-"below average"
df2$col<-factor(df2$col, levels=c("above average","below average"))
CFP<-subset(df2, scen == scens[1] )
CF1<-subset(df2, scen == scens[2] )
CF2<-subset(df2, scen == scens[3] )

a <- SPEI_annual_bar(CF1, title="Climate Future 1") + coord_cartesian(ylim = c(min(df2$SPI), max(df2$SPI)))
b <- SPEI_annual_bar(CF2, title="Climate Future 2") +  coord_cartesian(ylim = c(min(df2$SPI), max(df2$SPI)))

c <- var_bar_plot(Drought_all,"Duration", c("grey",cols), "Duration", "Years")
d <- var_bar_plot(Drought_all,"Drt.Free", c("grey",cols), "Drought-free\ninterval", 
                  "Years")
e<- var_bar_plot(Drought_all,"Severity", c("grey",cols), "Severity", 
                 "Severity \n(Intensity * Duration)")+ coord_cartesian(ylim = c(0, min(Drought_all$Severity)))

spei.time <- grid_arrange_shared_legend(a + rremove("ylab") + rremove("x.text"),b +  rremove("ylab"),
                                        nrow=2,ncol=1,position="bottom")

spei.time <- annotate_figure(spei.time, left = textGrob("SPI", rot = 90, vjust = 1, gp = gpar(cex = 2)))

drt.char <- grid.arrange(c+rremove("x.text"),d+rremove("x.text"),e,nrow=3,
                         top = textGrob(paste0(climate_zones[z], "-Average \ndrought characteristics"),gp=gpar(fontface="bold", col="black", fontsize=20,hjust=0.5)))

g <- grid.arrange(spei.time, drt.char,ncol = 2, clip = FALSE)
ggsave(plot=g,paste0(climate_zones[z],"DroughtCharacteristics-1-Panel.png"),path = plot.dir, height=9, width=12,bg = 'white')


# Option 2
c <- c+ theme(legend.title=element_text(size=24),legend.text=element_text(size=22),legend.position = "bottom")
d <- d+ theme(legend.title=element_text(size=24),legend.text=element_text(size=22),legend.position = "bottom")
e <- e+ theme(legend.title=element_text(size=24),legend.text=element_text(size=22),legend.position = "bottom")

drt.char <-grid_arrange_shared_legend(c+ rremove("x.text"),d+ rremove("x.text"),e+ rremove("x.text"),
                                      ncol=3,nrow=1,position="bottom",
                                      top = textGrob(paste0(climate_zones[z], "-Average drought characteristics"),gp=gpar(fontface="bold", col="black", fontsize=26,hjust=0.5)))
g <- grid.arrange(spei.time, drt.char,nrow=2,ncol = 1, clip = FALSE)
ggsave(plot=g,paste0(climate_zones[z],"DroughtCharacteristics-2-Panel.png"),path = plot.dir, height=9, width=12,bg = 'white')
}
