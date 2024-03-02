rm(list=ls())

wd<-"D:/PhD Project/Feedback paper_Albedo/Data and codes/"

if(!require(raster)){ install.packages("raster");library(raster)}
if(!require(rgdal)){ install.packages("rgdal");library(rgdal)} #for spTransform
if(!require(ncdf4)){ install.packages("ncdf4");library(ncdf4)}
if(!require(reshape2)){ install.packages("reshape2");library(reshape2)}

source("D:/PhD Project/Feedback paper_Albedo/Data and codes/Feedback paper_Albedo_Functions.R", encoding = 'UTF-8')

#########################################################################################
############################### Load data
#########################################################################################

############################### albedo

albedo_raster<-brick(paste(wd,"Input data/Decomposition of albedo/monthly albedo 2005 0.1Deg.nc",sep=""))
albedo_raster <- rasterToPoints(albedo_raster, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
albedo_raster <- spTransform(albedo_raster, CRS(albedo_raster@proj4string@projargs))# reproject sp object
albedo.all <- data.frame(albedo_raster@data, lon=coordinates(albedo_raster)[,1],lat=coordinates(albedo_raster)[,2])    # Assign coordinates to @data slot                     
colnames(albedo.all)<-c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12","lon","lat")
albedo.all$lon<-ifelse(albedo.all$lon<=180,albedo.all$lon,albedo.all$lon-360)

albedo_monthly<-albedo.all
albedo_monthly$lon<-ceiling(albedo_monthly$lon/5.625) *5.625
albedo_monthly$lat<-ceiling(albedo_monthly$lat/5.625) *5.625
albedo_monthly<-aggregate(. ~ lon + lat, data = albedo_monthly, FUN = mean)

albedo_seasonal<-cbind.data.frame(lon=albedo_monthly[,"lon"],
                                  lat=albedo_monthly[,"lat"],
                                  season1=rowSums(albedo_monthly[,c("m12","m1","m2")])/3,
                                  season2=rowSums(albedo_monthly[,c("m3","m4","m5")])/3,
                                  season3=rowSums(albedo_monthly[,c("m6","m7","m8")])/3,
                                  season4=rowSums(albedo_monthly[,c("m9","m10","m11")])/3)
albedo_seasonal<-melt(albedo_seasonal,id.vars = c("lon","lat"))
colnames(albedo_seasonal)=c("lon","lat","season","albedo")

############################ snow_cover

snow_cover_raster<-brick(paste(wd,"Input data/Decomposition of albedo/monthly snow cover 2005 0.1Deg.nc",sep=""))
snow_cover_raster <- rasterToPoints(snow_cover_raster, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
snow_cover_raster <- spTransform(snow_cover_raster, CRS(snow_cover_raster@proj4string@projargs))# reproject sp object
snow_cover.all <- data.frame(snow_cover_raster@data, lon=coordinates(snow_cover_raster)[,1],lat=coordinates(snow_cover_raster)[,2])    # Assign coordinates to @data slot                     
colnames(snow_cover.all)<-c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12","lon","lat")
snow_cover.all$lon<-ifelse(snow_cover.all$lon<=180,snow_cover.all$lon,snow_cover.all$lon-360)
snow_cover.all[,c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12")]<-snow_cover.all[,c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12")]/100

snow_cover_monthly<-snow_cover.all
snow_cover_monthly$lon<-ceiling(snow_cover_monthly$lon/5.625) *5.625
snow_cover_monthly$lat<-ceiling(snow_cover_monthly$lat/5.625) *5.625
snow_cover_monthly<-aggregate(. ~ lon + lat, data = snow_cover_monthly, FUN = mean)

snow_cover_seasonal<-cbind.data.frame(lon=snow_cover_monthly[,"lon"],
                                      lat=snow_cover_monthly[,"lat"],
                                      season1=rowSums(snow_cover_monthly[,c("m12","m1","m2")])/3,
                                      season2=rowSums(snow_cover_monthly[,c("m3","m4","m5")])/3,
                                      season3=rowSums(snow_cover_monthly[,c("m6","m7","m8")])/3,
                                      season4=rowSums(snow_cover_monthly[,c("m9","m10","m11")])/3)
snow_cover_seasonal<-melt(snow_cover_seasonal,id.vars = c("lon","lat"))
colnames(snow_cover_seasonal)=c("lon","lat","season","snow_cover")

##################################### fapar

load(paste(wd,"Input data/Decomposition of albedo/fapar95percentile 1982~2016 0.5Deg.rda",sep=""))
fapar_mean<-apply(fapar95percentile, c(1,2), mean) #mean annual maximum fapar over the period 1982~2016
fapar_mean<-as.data.frame(fapar_mean)
lon<-seq(-179.75,179.75,by=0.5)
lat<- -seq(-89.75,89.75,by=0.5)
colnames(fapar_mean)<-lon
fapar_mean$lat<-lat

if(!require(reshape2)){ install.packages("reshape2");library(reshape2)}
fapar<-melt(fapar_mean,id.vars = c("lat"))
colnames(fapar)<-c("lat","lon","fapar")
fapar$lon<-as.numeric(as.character(fapar$lon))
fapar$lon<-ceiling(fapar$lon/5.625) *5.625
fapar$lat<-ceiling(fapar$lat/5.625) *5.625
fapar<-aggregate(. ~ lon + lat, data = fapar, FUN = mean)

##################################### vegetation height

hveg_raster<-brick(paste(wd,"Input data/Decomposition of albedo/Tree height 2005 0.5Deg.nc",sep=""))
hveg_raster <- rasterToPoints(hveg_raster, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
hveg_raster <- spTransform(hveg_raster, CRS(hveg_raster@proj4string@projargs))# reproject sp object
hveg.all <- data.frame(hveg_raster@data, lon=coordinates(hveg_raster)[,1],lat=coordinates(hveg_raster)[,2])    # Assign coordinates to @data slot                     
colnames(hveg.all)<-c("hveg","lon","lat")

hveg<-hveg.all
hveg$lon<-ceiling(hveg$lon/5.625) *5.625
hveg$lat<-ceiling(hveg$lat/5.625) *5.625
hveg<-aggregate(. ~ lon + lat, data = hveg, FUN = mean) 

#################################  Build data frame

identical(albedo_seasonal[,c("lon","lat","season")],snow_cover_seasonal[,c("lon","lat","season")])

df<-cbind.data.frame(albedo_seasonal,snow_cover=snow_cover_seasonal[,"snow_cover"]) 
df<-merge(df,fapar,by=c("lon","lat"))
df<-merge(df,hveg,by=c("lon","lat"))

df<-na.omit(df)
write.csv(df,paste(wd,"Output data/Decomposition of albedo/Decomposition data.csv",sep=""))



#########################################################################################
############################### Decomposition of albedo
#########################################################################################

df<-read.csv(paste(wd,"Output data/Decomposition of albedo/Decomposition data.csv",sep=""),row.names=1)
df$logit_albedo<-log(df$albedo/(1-df$albedo))

f <- function(fapar,snow_cover,hveg,b,kf,ks,kh) {
  b+kf*(1-fapar)+ks*snow_cover+kh*hveg
}
fit <- nls(logit_albedo ~ f(fapar,snow_cover,hveg, b,kf,ks,kh), data=df,start = list(b=1,kf=1,ks=1,kh=-1))
summary(fit)
fitted<-f(df$fapar,df$snow_cover,df$hveg,b=coef(fit)[1],kf=coef(fit)[2],ks=coef(fit)[3],kh=coef(fit)[4])
parameter_info_albedo=summary(fit)[["coefficients"]][,c("Estimate","Std. Error")]
write.csv(parameter_info_albedo,paste(wd,"Output data/Decomposition of albedo/parameter_info_albedo.csv",sep=""))

#plot the model

pred1<-cbind.data.frame(df$lon,df$lat,df$logit_albedo,fitted);
colnames(pred1)<-c("lon","lat","observed_logit_albedo","fitted_logit_albedo")
p1<-ggplot(data=pred1,aes(observed_logit_albedo,fitted_logit_albedo))+theme_bw()+geom_point(size=0.8)+
  geom_abline(slope=1,intercept=0,col="red")+labs(x="logit albedo",y="fitted logit albedo")+geom_bin2d(bins=100)+
  scale_fill_gradientn(colours  = c("purple","blue","dodgerblue2","green","khaki","yellow","orange","red"),
                       na.value = "grey90")+theme(legend.position="bottom")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

pred2<-cbind.data.frame(df$lon,df$lat,df$albedo,1/(1+exp(-fitted)));
colnames(pred2)<-c("lon","lat","observed_albedo","fitted_albedo")
p2<-ggplot(data=pred2,aes(observed_albedo,fitted_albedo))+theme_bw()+geom_point(size=0.8)+
  geom_abline(slope=1,intercept=0,col="red")+labs(x="albedo",y="fitted albedo")+geom_bin2d(bins=100)+
  scale_fill_gradientn(colours  = c("purple","blue","dodgerblue2","green","khaki","yellow","orange","red"),
                       na.value = "grey90")+theme(legend.position="bottom")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

p<-ggarrange(p1,p2,ncol=2,labels=c("(a)","(b)"))

ggsave(file=paste(wd,"Output data/Decomposition of albedo/Fig.1 Decomposition of albedo.png",sep=""),
       p,width=8,height=4.8)

