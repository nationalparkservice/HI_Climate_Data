library(tidyr)
library(dplyr)
library(stars)
library(raster)
library(here)
library(ggplot2)
library(ggthemes)
library(zoo)
library(viridis)
library(ggbreak)
library(lemon)
library(ggpubr);library(gridExtra);library(grid);library(gtable)

rm(list=ls())

data.dir <- here::here('data/Output/Data-files//')
plot.dir <- here::here('data/Output//Plots//')

#ANNPrecip -- load .rds files
CF1.rds <- readRDS(paste0(data.dir,"WetPrecipDelta_rcp45")) %>%  #change var name
  mutate(pcp.in = mean /25.4) #%>% extract(pcp.in) 
CF1.rds <- CF1.rds["pcp.in"]
CF2.rds <- readRDS(paste0(data.dir,"WetPrecipDelta_rcp85")) %>%  #change var name
  mutate(pcp.in = mean /25.4)# %>% #select(pcp.in) 
CF2.rds <- CF2.rds["pcp.in"]

boundary <-st_read('C:/Users/arunyon/OneDrive - DOI/Documents/GIS/HAVO_Kilauea_Summit_Wet_Dry_Zones/HAVO_Kilauea_Summit_Wet_Dry_Zones.shp')
boundary <- st_transform(boundary, st_crs(CF1.rds))
CF_GCM <- data.frame(CF=c("Climate Future 1", "Climate Future 2"), scen=c("rcp45","rcp85"))
cols <- c("#9A9EE5","#E10720")

var = "PrecipIn" #change to name of var in df
long.title = "average wet-season precipitation (in/month)" #change to be legend
delta.var<- "PrecipIn"
scale="viridis"
seas="wet"

# insert topo
topo <- stack('./data/data/natehii0100a.tif')
topo <- projectRaster(topo,crs = crs(boundary)); topo <- crop(topo, boundary)
topo_df  <- as.data.frame(topo, xy = TRUE) 

# Generate sample data for ts plot
df = read.csv(paste0(plot.dir,"Seasonal-delta.csv")) 
df = subset(df, season == seas & CF!="Recent")

scale.min = min(c(CF1.rds$pcp.in, CF2.rds$pcp.in),na.rm=TRUE) #change names
scale.max = max(c(CF1.rds$pcp.in, CF2.rds$pcp.in),na.rm=TRUE) #change names

# ggplot
map.plot <- function(data, title,xaxis,metric,col){
  ggplot() +
    geom_raster(data = topo_df ,aes(x = x, y = y,alpha=natehii0100a_1), show.legend=FALSE) +
    geom_stars(data = data, alpha = 0.8) + 
    geom_sf(data = boundary, aes(), fill = NA,colour="black") + 
    scale_fill_viridis(direction=-1, option = scale, 
                       limits = c(scale.min, scale.max), oob = scales::squish) + 
    labs(title = title) +
    theme_map() +
    theme(legend.position = "bottom",
          legend.key.width = unit(6, "cm"),
          legend.key.height = unit(.3, "cm"),
          legend.justification = "center",
          # plot.title=element_blank(),
          plot.title=element_text(size=12,face="bold",hjust=0.5),
          plot.background = element_rect(colour = col, fill=NA, size=5),
          plot.margin = unit(c(.5,0,0,0), "cm")) + 
    labs(fill = metric)
}

# CF1
cf1.plot <- map.plot(data=CF1.rds, title="Climate Future 1",metric=paste0("Average ",long.title),col=cols[1])
cf2.plot <- map.plot(data=CF2.rds, title="Climate Future 2",metric=paste0("Average ",long.title),col=cols[2])

###########################################################
##### NONE OF THESE WORK ##################################
# Merge into one plot

maps <- grid_arrange_shared_legend(cf1.plot, cf2.plot,ncol = 1, nrow = 2, position = "bottom", 
                                   top = textGrob(paste0("Change in ",long.title),
                                                  gp=gpar(fontface="bold", col="black", fontsize=16)))


################################### MONTHLY DOT PLOT ##################


dotplot <- ggplot(df, aes(x=Precip.delta,y=zone,fill=CF)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black") + 
  geom_point(stat="identity",size=8,colour="black",aes(fill = factor(CF), shape = factor(CF))) +
  theme(axis.text=element_text(size=16),    #Text size for axis tick mark labels
        axis.title.x=element_blank(),               #Text size and alignment for x-axis label
        plot.title=element_blank(),
        # axis.title.y=element_text(size=16, vjust=0.5,  margin=margin(t=20, r=20, b=20, l=20)),              #Text size and alignment for y-axis label
        # plot.title=element_text(size=20, vjust=0.5, face="bold", margin=margin(t=20, r=20, b=20, l=20)),      #Text size and alignment for plot title
        legend.title=element_text(size=16),                                                                    #Text size of legend category labels
        legend.text=element_text(size=14),                                                                   #Text size of legend title
        legend.position = "bottom")  +
  labs(title = paste0("Change in seasonal ",long.title), 
       x = "Change (°F)", y = "") +
  scale_fill_manual(name="",values =cols) +
  scale_shape_manual(name="",values = c(21,22)) +
  scale_y_discrete(limits=rev)
dotplot

g <- grid.arrange(maps, dotplot,ncol = 2, widths = c(4, 4), clip = FALSE)

annotate_figure(g, top = text_grob(paste0("Change in seasonal ",long.title, "; 1990-2009 vs 2080-2099"), 
                                   face = "bold", size = 20))

ggsave(paste0(seas,"_season_",var,".png"), width = 15, height = 9, path = plot.dir,bg="white")
