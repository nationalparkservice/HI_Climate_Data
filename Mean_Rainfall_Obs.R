# Calculates ANN and Monthly mean rainfall from Frazier 2016 dataset from years 1990-2009
## This is hacked together and could be written 100% better but just need job done

library(raster)
library(dplyr)
library(stars)

rm(list=ls())

Obs.data<-"D:/HI_Data/OBS/RFMonthYr_Rasters_State_mm_1920_2012/Month_Rasters_State_mm/"

#Annual
names<-list.dirs(path = paste0(Obs.data,"Annual/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009

ANN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(ANN_mean)
writeRaster(ANN_mean, paste0(Obs.data,"/ANN_mean_1990-2009"), format = "GTiff")

#Jan
names<-list.dirs(path = paste0(Obs.data,"State_01Jan/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009

JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/JAN_mean_1990-2009"), format = "GTiff")

#Feb
names<-list.dirs(path = paste0(Obs.data,"State_02Feb/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009

JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/FEB_mean_1990-2009"), format = "GTiff")

#Mar
names<-list.dirs(path = paste0(Obs.data,"State_03Mar/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009

JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/MAR_mean_1990-2009"), format = "GTiff")

#Apr
names<-list.dirs(path = paste0(Obs.data,"State_04Apr/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009

JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/APR_mean_1990-2009"), format = "GTiff")

#May
names<-list.dirs(path = paste0(Obs.data,"State_05May/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009

JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/MAY_mean_1990-2009"), format = "GTiff")

#Jun
names<-list.dirs(path = paste0(Obs.data,"State_06Jun/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009

JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/JUN_mean_1990-2009"), format = "GTiff")

#Jul
names<-list.dirs(path = paste0(Obs.data,"State_07Jul/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009
names[1]
JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/JUL_mean_1990-2009"), format = "GTiff")

#Aug
names<-list.dirs(path = paste0(Obs.data,"State_08Aug/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009
names[1]
JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/AUG_mean_1990-2009"), format = "GTiff")

#Sep
names<-list.dirs(path = paste0(Obs.data,"State_09Sep/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009
names[1]
JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/SEP_mean_1990-2009"), format = "GTiff")

#Oct
names<-list.dirs(path = paste0(Obs.data,"State_10Oct/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009
names[1]
JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/OCT_mean_1990-2009"), format = "GTiff")

#Nov
names<-list.dirs(path = paste0(Obs.data,"State_11Nov/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009
names[1]
JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/NOV_mean_1990-2009"), format = "GTiff")

#Dec
names<-list.dirs(path = paste0(Obs.data,"State_12Dec/"), full.names = TRUE, recursive = F) # 72-91 are years 1990-2009
names[1]
JAN_mean<-mean(raster(names[72]),raster(names[73]),raster(names[74]),raster(names[75]),raster(names[76]),raster(names[77]),
               raster(names[78]),raster(names[79]),raster(names[80]),raster(names[81]),raster(names[82]),raster(names[83]),
               raster(names[84]),raster(names[85]),raster(names[86]),raster(names[87]),raster(names[88]),raster(names[89]),
               raster(names[90]),raster(names[91]))
plot(JAN_mean)
writeRaster(JAN_mean, paste0(Obs.data,"/DEC_mean_1990-2009"), format = "GTiff")
