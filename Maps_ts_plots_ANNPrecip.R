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

data.dir <- "C:/Users/achildress/Documents/Git-repos/HI_Climate_Data/data/Output/Data-files/"
plot.dir <- "C:/Users/achildress/Documents/Git-repos/HI_Climate_Data/data/Output/Plots/"

#ANNPrecip -- load .rds files
CF1.rds <- readRDS(paste0(data.dir,"ANNPrecipDelta_rcp45")) %>% #change var name
  mutate(pcp.in = Rainfall_mm /25.4) %>% select(pcp.in) 
CF2.rds <- readRDS(paste0(data.dir,"ANNPrecipDelta_rcp85")) %>% #change var name
  mutate(pcp.in = Rainfall_mm /25.4) %>% select(pcp.in) 
  
boundary <-boundary <- st_read('./data/HALE/HALE_boundary.shp')
boundary <- st_transform(boundary, st_crs(CF1.rds))
CF_GCM <- data.frame(CF=c("Climate Future 1", "Climate Future 2"), scen=c("rcp45","rcp85"))
cols <- c("#9A9EE5","#E10720")

var = "PrecipIn" #change to name of var in df
long.title = "total precipitation (in/year)" #change to be legend
delta.var = "PrecipIn"
scale="viridis"

# insert topo
topo <- stack('./data/HALE/HALENatEa1.tif')
topo_df  <- as.data.frame(topo, xy = TRUE) 

# Generate sample data for ts plot
df = read.csv(paste0(data.dir,"RF.monthly.HALE.csv")) %>% mutate(PrecipIn = Rainfall_mm/25.4)
df = merge(df, CF_GCM,by="scen",all=TRUE)
df$CF[which(is.na((df$CF)))] = "Recent"
# df$CF_col[which(is.na((df$CF_col)))] = "grey"
df$CF = factor(df$CF, levels=c("Recent",CF_GCM$CF))
df$Year = as.Date(format(as.Date(df$date,format="%Y-%m-%d"),"%Y"),format="%Y")
DF=aggregate(PrecipIn~Year+CF,df,sum)
DF$period = factor(ifelse(DF$Year<"2010-01-01","Past","Future"),levels=c("Past","Future"))
  
  means <- DF %>% group_by(CF) %>%
    summarize(var = mean(eval(parse(text=var)))) 
  
  scale.min = min(c(CF1.rds$pcp.in, CF2.rds$pcp.in),na.rm=TRUE)
  scale.max = max(c(CF1.rds$pcp.in, CF2.rds$pcp.in),na.rm=TRUE)
  
  # ggplot
  map.plot <- function(data, title,xaxis,metric,col){
    ggplot() + 
      geom_raster(data = topo_df ,aes(x = x, y = y,alpha=HALENatEa1_1), show.legend=FALSE) +
      geom_stars(data = data, alpha = 0.8) + 
      geom_sf(data = boundary, aes(), fill = NA) +
      scale_fill_viridis(direction=-1, option = scale, limits = c(scale.min, scale.max),  
                         guide = guide_colorbar(title.position = "top", title.hjust = 0.5),oob = scales::squish) + #mako for WB delta
      labs(title = title) +
      theme_map() +
      theme(legend.position = "bottom",
            legend.key.width = unit(6, "cm"),
            legend.key.height = unit(.3, "cm"),
            legend.justification = "center",
            plot.title=element_text(size=12,face="bold",hjust=0.5),
            plot.background = element_rect(colour = col, fill=NA, linewidth = 5)) + 
      labs(fill = metric)
  }
  cf1 <- map.plot(data=CF1.rds,title="Climate Future 1",metric=long.title,col=cols[1])
  cf2 <- map.plot(data=CF2.rds,title="Climate Future 2",metric=long.title,col=cols[2])
  
 ts<-ggplot(DF, aes(x=Year, y=(eval(parse(text=var))), group=CF, colour = CF)) +
    geom_line(colour = "black",size=2.5, stat = "identity") +
    geom_line(size = 2, stat = "identity") +
    geom_point(colour= "black", size=4, aes(fill = factor(CF), shape = factor(CF))) +
    theme(axis.text=element_text(size=14),
          # axis.text.x=element_blank(),
          axis.title.x=element_text(size=16,vjust=1.0),
          axis.title.y=element_text(size=16,vjust=1.0),
          plot.title=element_blank(),
          legend.text=element_text(size=14), legend.title=element_text(size=14),
          legend.position = "bottom") +
    labs(title = paste0("Change in annual ",long.title), 
         x = "Year", y = long.title) +
    scale_color_manual(name="",values = c("grey",cols)) +
    scale_fill_manual(name="",values = c("grey",cols)) +
    scale_shape_manual(name="",values = c(21,22,23)) +
    facet_wrap(~period, nrow = 1, ncol=3,scales = "free_x") 
    # coord_fixed(ratio = .3) 
  ts
  
  
  #### Just maps and ts plot
  maps <- grid_arrange_shared_legend(cf1, cf2,ncol = 2, nrow = 1, position = "bottom", 
                                     top = textGrob(paste0("Change in ",long.title),
                                                    gp=gpar(fontface="bold", col="black", fontsize=16)))
  # g <- ggarrange(maps,ts, nrow=2)
  # g
  
  #### Maps, ts, table
  delta.var <- means
  delta.var$var[2:3] <- delta.var$var[2:3] - delta.var$var[1]
  delta.var$var <- signif(delta.var$var, digits = 1)
  
  table <- tableGrob(delta.var, rows = NULL,cols=NULL)
  # table <- gtable_add_grob(table, grobs = rectGrob(gp = gpar(fill=NA, lwd=2)), #library(gtable)
  #                          t=4,b=nrow(table),l=1,r=ncol(table))
  table <- annotate_figure(table,
                           top = text_grob("Recent = absolute value \n CFs = change values", color = "black",
                                           face = "italic", size = 12))
  tsplots <- grid.arrange(ts, table,ncol = 2, widths = c(4, 1), clip = FALSE)
  
  
  g <- ggarrange(maps,tsplots, nrow=2)+bgcolor("white") 
  g
  

  ggsave(paste0(var,"_ANN.png"), width = 20, height = 9, path = plot.dir)
  
  