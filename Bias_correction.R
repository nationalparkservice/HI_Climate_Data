#Extract January data

P<-readRDS(paste0(here::here('data/Output/Data-files//'),'RF.monthly-present'))
R45<-readRDS(paste0(here::here('data/Output/Data-files//'),'RF.monthly-rcp45'))

monthly.index<-read.csv(paste0(here::here('data/Monthly_indexing.csv')))

R45.1<-st_apply(R45[,,,c(monthly.index[,1])],c("x","y"),min) 
plot(R45.1)



