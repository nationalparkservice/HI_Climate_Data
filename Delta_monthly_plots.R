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

CF_GCM <- data.frame(CF=c("Climate Future 1", "Climate Future 2"), scen=c("rcp45","rcp85"))
cols <- c("#9A9EE5","#E10720")
scens<-c("present","rcp45","rcp85")
MonthLabels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

# Generate sample data for ts plot
df = read.csv(paste0(data.dir,"Zone-Monthly-data.csv")) 
df = merge(df, CF_GCM,by="scen",all=TRUE) 
df$CF[which(is.na((df$CF)))] = "Recent"
df$CF = factor(df$CF, levels=c("Recent",CF_GCM$CF))
df$Month = as.numeric(format(as.Date(df$date,format="%Y-%m-%d"),"%m"))
df$season<-ifelse(df$Month<5| df$Month>10,"wet","dry")

df.present <- subset(df, CF=="Recent" & season == "dry")
df.present <- aggregate(cbind(TmeanF,PrecipIn)~CF+zone ,df.present,mean) #change var
df.future <- subset(df, CF!="Recent" & season == "dry")
df.future <- aggregate(cbind(TmeanF,PrecipIn)~CF+zone ,df.future,mean) #change var

df.future %>% group_by(CF) %>% mutate(Tmean.delta = TmeanF - df.present$TmeanF,
                                          Precip.delta = PrecipIn - df.present$PrecipIn,
                                      season = "dry")-> df.future.dry
df.present.dry <- df.present %>% mutate(Tmean.delta=NA,Precip.delta=NA,season="dry")

df.present <- subset(df, CF=="Recent" & season == "wet")
df.present <- aggregate(cbind(TmeanF,PrecipIn)~CF+zone ,df.present,mean) #change var
df.future <- subset(df, CF!="Recent" & season == "wet")
df.future <- aggregate(cbind(TmeanF,PrecipIn)~CF+zone ,df.future,mean) #change var

df.future %>% group_by(CF) %>% mutate(Tmean.delta = TmeanF - df.present$TmeanF,
                                      Precip.delta = PrecipIn - df.present$PrecipIn,
                                      season="wet") -> df.future.wet
df.present.wet <- df.present %>% mutate(Tmean.delta=NA,Precip.delta=NA,season="wet")

DF = rbind(df.future.dry, df.present.dry, df.future.wet,df.present.wet)
write.csv(DF,file=paste0(plot.dir,"Seasonal-delta.csv"),row.names=F)


df.monthly <- aggregate(cbind(TmeanF, PrecipIn)~CF+zone+Month,df,mean)
df.monthly

monthly.present <- subset(df.monthly, CF=="Recent")
monthly.future1<- subset(df.monthly, CF=="Climate Future 1")
monthly.future2<- subset(df.monthly, CF=="Climate Future 2")

monthly.future1$DeltaTmean <- monthly.future1$TmeanF - monthly.present$TmeanF
monthly.future2$DeltaTmean <- monthly.future2$TmeanF - monthly.present$TmeanF

monthly.future1$DeltaPrecip <- monthly.future1$PrecipIn - monthly.present$PrecipIn
monthly.future2$DeltaPrecip <- monthly.future2$PrecipIn - monthly.present$PrecipIn

monthly.future.delta <- rbind(monthly.future1,monthly.future2)

write.csv(monthly.future.delta,file=paste0(plot.dir,"Monthly-delta.csv"),row.names=F)

monthly.future.delta <- read.csv(paste0(plot.dir,"Monthly-delta.csv"))

PlotTheme = theme(axis.text=element_text(size=14),    #Text size for axis tick mark labels
                  axis.title.x=element_blank(),               #Text size and alignment for x-axis label
                  axis.title.y=element_text(size=18, vjust=0.5,  margin=margin(t=20, r=20, b=20, l=20)),              #Text size and alignment for y-axis label
                  plot.title=element_text(size=20,face="bold",hjust=0.5, margin=margin(t=20, r=20, b=20, l=20)),      #Text size and alignment for plot title
                  legend.title=element_text(size=18),                                                                    #Text size of legend category labels
                  legend.text=element_text(size=16),                                                                   #Text size of legend title
                  legend.position = "bottom")  

Month_bar_plot <- function(data, xvar, yvar, grp, cols, title,xlab, ylab,labels,CFmethod=""){
  ggplot(data, aes(x={{xvar}},y={{yvar}},fill={{grp}})) +
    geom_bar(stat="identity",position="dodge",colour="black") +
    PlotTheme +
    labs(title = title, 
         x = xlab, y = ylab) +
    scale_fill_manual(name="",values = cols) +
    scale_x_discrete(labels = labels) +
    facet_wrap(~zone, nrow = 1, ncol=4,scales = "free_x") 
  # annotate(geom="text", x=Inf, y=-Inf, label=CFmethod,color="black",vjust=-1,hjust=1)
}

Month_line_plot <- function(data, xvar, yvar, grp, cols, title,xlab, ylab){
  ggplot(data, aes(x={{xvar}}, y={{yvar}}, group={{grp}}, colour = {{grp}})) +
    geom_line(linewidth = 2, stat = "identity",colour="black") + 
    geom_line(linewidth = 1.5, stat = "identity") +
    geom_point(colour= "black", size=4, aes(fill = factor({{grp}}), shape = factor({{grp}}))) +
    PlotTheme +
    labs(title = title,
         x = xlab, y = ylab) +
    scale_color_manual(name="",values = cols) +
    scale_fill_manual(name="",values = cols) +
    scale_shape_manual(name="",values = c(seq(21,21+length(cols)-1,1))) +
    # scale_y_continuous(limits=c(0, ceiling(max(eval(parse(text=paste0(data,"$",yvar))))))) +
    scale_x_discrete(labels = MonthLabels) +
    facet_wrap(~zone, nrow = 1, ncol=4,scales = "free_x")
  # annotate(geom="text", x=Inf, y=-Inf, label=CFmethod,color="black",vjust=-1,hjust=1) 
}

Month_bar_plot(monthly.future.delta,xvar=Month,yvar=DeltaPrecip,grp=CF,cols=cols,
               title="Change in avg monthly precipitation \nin 2090 vs 2000",
               xlab = "Monthly", ylab="Change in precipitation (in)",label=MonthLabels)

ggsave("Monthly.precip.bar.png", width = 15, height = 9, path = plot.dir)


Month_line_plot(monthly.future.delta,xvar=Month,yvar=DeltaTmean,grp=CF,cols=cols,
               title="Change in avg monthly temperature \nin 2090 vs 2000",
               xlab = "Monthly", ylab="Change in temperature (\u00B0F)")

ggsave("Monthly.tmean.bar.png", width = 15, height = 9, path = plot.dir)
