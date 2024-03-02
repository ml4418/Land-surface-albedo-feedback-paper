rm(list=ls())

wd<-"D:/PhD Project/Feedback paper_Albedo/Data and codes/"

if(!require(raster)){ install.packages("raster");library(raster)}
if(!require(rgdal)){ install.packages("rgdal");library(rgdal)}# for spTransform
if(!require(ncdf4)){ install.packages("ncdf4");library(ncdf4)}
if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(egg)){ install.packages("egg");library(egg)}
if(!require(sf)){ install.packages("sf");library(sf)}
if(!require(reshape2)){ install.packages("reshape2");library(reshape2)}

source("D:/PhD Project/Feedback paper_Albedo/Data and codes/Feedback paper_Albedo_Functions.R", encoding = 'UTF-8')

#######################################################################################
#####################  Load age adjustment for LOVECLIM simulations
#######################################################################################

#get the start dates identified in Liu et al. (2022) for DO 5~12
DO_timing <- read.csv(paste(wd,"Input data/Timing/store_GHG paper.csv",sep=""), row.names=1, stringsAsFactors = FALSE)
DO_timing<-DO_timing$tmin_T

#get the official start dates provided by Wolff et al. (2010) for DO 5~12
DO_timing1 <- read.csv(paste(wd,"Input data/Timing/DO_timing1.csv",sep=""), row.names=1, stringsAsFactors = FALSE)
DO_timing1<-DO_timing1[5:12,]

#get the relationship between the two timings
plot(DO_timing1~DO_timing)
sum<-summary(lm(DO_timing1~DO_timing));print(sum)
b0_t<-sum[["coefficients"]]["(Intercept)","Estimate"]
b1_t<-sum[["coefficients"]]["DO_timing","Estimate"]
abline(a=b0_t,b=b1_t)
round(b0_t,digits=0)
round(b1_t,digits=2)

########################################################################################
########################################################################################
#######################                                           ######################
#######################   1. LOVECLIM temperature simulations     ######################
#######################                                           ######################
########################################################################################
########################################################################################

####################### Load T
LOVECLIM_T<-brick(paste(wd,"Input data/LOVECLIM/LOVECLIM temperature.nc",sep=""))
LOVECLIM_T <- rasterToPoints(LOVECLIM_T, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
LOVECLIM_T <- spTransform(LOVECLIM_T, CRS(LOVECLIM_T@proj4string@projargs))# reproject sp object
LOVECLIM_T.all <- data.frame(lon=coordinates(LOVECLIM_T)[,1],lat=coordinates(LOVECLIM_T)[,2],LOVECLIM_T@data)    # Assign coordinates to @data slot                     
age<-c(50000:29801)
colnames(LOVECLIM_T.all)<-c("lon","lat",age)

LOVECLIM_T.all <- melt(LOVECLIM_T.all,id.vars = c("lon","lat"))
colnames(LOVECLIM_T.all)<-c("lon","lat","age","T")
LOVECLIM_T<-LOVECLIM_T.all
LOVECLIM_T$lon<-as.numeric(LOVECLIM_T$lon)
LOVECLIM_T$lat<-as.numeric(LOVECLIM_T$lat)
LOVECLIM_T$age<-as.numeric(as.character(LOVECLIM_T$age))
LOVECLIM_T$age<-b0_t+b1_t*LOVECLIM_T$age #adjust age scale

####################### Get land and ocean part
if(!require(maptools)){ install.packages("maptools");library(maptools)}
data(wrld_simpl)
set.seed(0)#Create a SpatialPoints object
points<-unique(LOVECLIM_T[,c("lon","lat")])
wrld_simpl@proj4string@projargs<-"+proj=longlat +datum=WGS84 +no_defs"
pts <- SpatialPoints(points, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
wrld<-over(pts, wrld_simpl)# Find which points fall over land except Antarctica and Greenland
wrld<-cbind.data.frame(points,wrld)

#get icefree land part
LOVECLIM_T_land<-merge(LOVECLIM_T,wrld,by=c("lon","lat"))
LOVECLIM_T_land<-LOVECLIM_T_land[!is.na(LOVECLIM_T_land$FIPS),]#get land points
LOVECLIM_T_land<-LOVECLIM_T_land[(LOVECLIM_T_land$NAME!="Greenland")&(LOVECLIM_T_land$NAME!="Antarctica"),]#exclude Greenland and Antarctica
LOVECLIM_T_land<-LOVECLIM_T_land[,colnames(LOVECLIM_T)]

######################### Bin in 25 years and 5.625 degree

LOVECLIM_T$age<-ceiling(LOVECLIM_T$age/25)*25 #bin in 25 years
LOVECLIM_T$lon<-ceiling(LOVECLIM_T$lon/5.625)*5.625 #bin in 5.625 degree
LOVECLIM_T$lat<-ceiling(LOVECLIM_T$lat/5.625)*5.625 #bin in 5.625 degree
LOVECLIM_T_binned<-aggregate(data=LOVECLIM_T,.~lon+lat+age,FUN=mean, na.action = na.omit) #get the mean in each bin
LOVECLIM_T_binned$sd_T<-aggregate(data=LOVECLIM_T,.~lon+lat+age,FUN=sd, na.action = na.omit)[,"T"] #get the sd in each bin
write.csv(LOVECLIM_T_binned,paste(wd,"Output data/LOVECLIM/LOVECLIM_T_binned.csv",sep=""))

LOVECLIM_T_land$age<-ceiling(LOVECLIM_T_land$age/25)*25 #bin in 25 years
LOVECLIM_T_land$lon<-ceiling(LOVECLIM_T_land$lon/5.625)*5.625 #bin in 5.625 degree
LOVECLIM_T_land$lat<-ceiling(LOVECLIM_T_land$lat/5.625)*5.625 #bin in 5.625 degree
LOVECLIM_T_land_binned<-aggregate(data=LOVECLIM_T_land,.~lon+lat+age,FUN=mean, na.action = na.omit) #get the mean in each bin
LOVECLIM_T_land_binned$sd_T<-aggregate(data=LOVECLIM_T_land,.~lon+lat+age,FUN=sd, na.action = na.omit)[,"T"] #get the sd in each bin
write.csv(LOVECLIM_T_land_binned,paste(wd,"Output data/LOVECLIM/LOVECLIM_T_binned_land.csv",sep=""))
write.csv(LOVECLIM_T_land_binned,paste(wd,"Output data/3D var/3D var input/LOVECLIM_T_binned_land.csv",sep=""))

##################################################################################################
##################################################################################################
################                                                                  ################
################   2. LOVECLIM surface solar radiation simulations (seasonal)     ################
################                                                                  ################
##################################################################################################
##################################################################################################

gc()
options(future.globals.maxSize = 6 * 1024^8)

############################ Load ssr

ssr1<-brick(paste(wd,"Input data/LOVECLIM/LOVECLIM_seasonal surface solar radiation_50-416.nc",sep=""))
ssr1 <- rasterToPoints(ssr1, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
ssr1 <- spTransform(ssr1, CRS(ssr1@proj4string@projargs))# reproject sp object
ssr1.all <- data.frame(lon=coordinates(ssr1)[,1],lat=coordinates(ssr1)[,2],ssr1@data)    # Assign coordinates to @data slot                     
age1<-rep(c(50000:(41600+1)),each=4)
age_miss<-seq(from=28800,to=50000,by=200)
age1<-age1[!(age1%in%age_miss)]
season<-rep(c(1:4),times=length(age1)/4)
age1<-paste(age1,season)
colnames(ssr1.all)<-c("lon",'lat',age1)
ssr1.all$lon<-ifelse(ssr1.all$lon<=180,ssr1.all$lon,ssr1.all$lon-360)

ssr2<-brick(paste(wd,"Input data/LOVECLIM/LOVECLIM_seasonal surface solar radiation_417-30.nc",sep=""))
ssr2 <- rasterToPoints(ssr2, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
ssr2 <- spTransform(ssr2, CRS(ssr2@proj4string@projargs))# reproject sp object
ssr2.all <- data.frame(lon=coordinates(ssr2)[,1],lat=coordinates(ssr2)[,2],ssr2@data)    # Assign coordinates to @data slot                     
age2<-rep(c(41700:(28900+1)),each=4)
age_miss<-seq(from=28800,to=50000,by=200)
age2<-age2[!(age2%in%age_miss)]
season<-rep(c(1:4),times=length(age2)/4)
age2<-paste(age2,season)
colnames(ssr2.all)<-c("lon","lat",age2)
ssr2.all$lon<-ifelse(ssr2.all$lon<=180,ssr2.all$lon,ssr2.all$lon-360)

replicate<-rep(c(41700:41601),each=4)
replicate<-paste(replicate,rep(c(1:4),times=length(replicate)/4))
ssr<-cbind.data.frame(ssr1.all[,-which(colnames(ssr1.all)%in%replicate)],ssr2.all[,-c(1,2)])
ssr <- melt(ssr,id.vars = c("lon","lat"))
colnames(ssr)<-c("lon","lat","age_and_season","ssr")
ssr$age_and_season<-as.character(ssr$age_and_season)
d_ssr<-strsplit(ssr$age_and_season, split = " ")
dd_ssr<-data.frame(matrix(unlist(d_ssr), ncol=2, byrow=TRUE))
ssr$age<-as.numeric(dd_ssr$X1)
ssr$season<-dd_ssr$X2

ssr$season<-as.numeric(ssr$season)
ssr<-ssr[,c("lon","lat","age","season","ssr")]
ssr_cast=dcast(ssr,lon+lat+age~season, value.var = "ssr") 

# Get land part
if(!require(maptools)){ install.packages("maptools");library(maptools)}
data(wrld_simpl)
set.seed(0)#Create a SpatialPoints object
points<-unique(ssr_cast[,c("lon","lat")])
wrld_simpl@proj4string@projargs<-"+proj=longlat +datum=WGS84 +no_defs"
pts <- SpatialPoints(points, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
wrld<-over(pts, wrld_simpl)# Find which points fall over land except Antarctica and Greenland
wrld<-cbind.data.frame(points,wrld)

#get icefree land part
ssr_cast_land<-merge(ssr_cast,wrld,by=c("lon","lat"))
ssr_cast_land<-ssr_cast_land[!is.na(ssr_cast_land$FIPS),]#get land points
ssr_cast_land<-ssr_cast_land[(ssr_cast_land$NAME!="Greenland")&(ssr_cast_land$NAME!="Antarctica"),]#exclude Greenland and Antarctica
ssr_cast_land<-ssr_cast_land[,colnames(ssr_cast)]

#Adjuste age and bin
LOVECLIM_seasonal_ssr<-ssr_cast_land
LOVECLIM_seasonal_ssr$age<-b0_t+b1_t*LOVECLIM_seasonal_ssr$age #adjust the age scale

LOVECLIM_seasonal_ssr$age<-ceiling(LOVECLIM_seasonal_ssr$age/25)*25 #bin in 25 years
LOVECLIM_seasonal_ssr$lon<-ceiling(LOVECLIM_seasonal_ssr$lon/5.625)*5.625 #bin in 5.625 degree
LOVECLIM_seasonal_ssr$lat<-ceiling(LOVECLIM_seasonal_ssr$lat/5.625)*5.625 #bin in 5.625 degree
mean_snow<-aggregate(data=LOVECLIM_seasonal_ssr,.~lon+lat+age,FUN=mean, na.action = na.omit) #get the mean in each bin
err_snow<-aggregate(data=LOVECLIM_seasonal_ssr,.~lon+lat+age,FUN=sd, na.action = na.omit) #get the mean in each bin

identical(mean_snow[,c("lon","lat","age")],err_snow[,c("lon","lat","age")])
LOVECLIM_seasonal_ssr_binned<-cbind.data.frame(mean_snow,err_snow[,-c(1:3)])
colnames(LOVECLIM_seasonal_ssr_binned)<-c("lon","lat","age","ssr_1","ssr_2","ssr_3","ssr_4",
                                          "err_ssr_1","err_ssr_2","err_ssr_3","err_ssr_4")

write.csv(LOVECLIM_seasonal_ssr_binned,paste(wd,"Output data/LOVECLIM/LOVECLIM_seasonal_ssr_binned_land.csv",sep=""))

##################################################################################################
##################################################################################################
#######################                                                     ######################
#######################       3. LOVECLIM albedo simulations (seasonal)     ######################
#######################                                                     ######################
##################################################################################################
##################################################################################################

gc()
options(future.globals.maxSize = 6 * 1024^8)

############################ Load albedo

albedo1<-brick(paste(wd,"Input data/LOVECLIM/LOVECLIM_seasonal albedo_50-416.nc",sep=""))
albedo1 <- rasterToPoints(albedo1, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
albedo1 <- spTransform(albedo1, CRS(albedo1@proj4string@projargs))# reproject sp object
albedo1.all <- data.frame(lon=coordinates(albedo1)[,1],lat=coordinates(albedo1)[,2],albedo1@data)    # Assign coordinates to @data slot                     
age1<-rep(c(50000:(41600+1)),each=4)
age_miss<-seq(from=28800,to=50000,by=200)
age1<-age1[!(age1%in%age_miss)]
season<-rep(c(1:4),times=length(age1)/4)
age1<-paste(age1,season)
colnames(albedo1.all)<-c("lon",'lat',age1)
albedo1.all$lon<-ifelse(albedo1.all$lon<=180,albedo1.all$lon,albedo1.all$lon-360)

albedo2<-brick(paste(wd,"Input data/LOVECLIM/LOVECLIM_seasonal albedo_417-30.nc",sep=""))
albedo2 <- rasterToPoints(albedo2, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
albedo2 <- spTransform(albedo2, CRS(albedo2@proj4string@projargs))# reproject sp object
albedo2.all <- data.frame(lon=coordinates(albedo2)[,1],lat=coordinates(albedo2)[,2],albedo2@data)    # Assign coordinates to @data slot                     
age2<-rep(c(41700:(28900+1)),each=4)
age_miss<-seq(from=28800,to=50000,by=200)
age2<-age2[!(age2%in%age_miss)]
season<-rep(c(1:4),times=length(age2)/4)
age2<-paste(age2,season)
colnames(albedo2.all)<-c("lon","lat",age2)
albedo2.all$lon<-ifelse(albedo2.all$lon<=180,albedo2.all$lon,albedo2.all$lon-360)

replicate<-rep(c(41700:41601),each=4)
replicate<-paste(replicate,rep(c(1:4),times=length(replicate)/4))
albedo<-cbind.data.frame(albedo1.all[,-which(colnames(albedo1.all)%in%replicate)],albedo2.all[,-c(1,2)])
albedo <- melt(albedo,id.vars = c("lon","lat"))
colnames(albedo)<-c("lon","lat","age_and_season","albedo")
albedo$age_and_season<-as.character(albedo$age_and_season)
d_albedo<-strsplit(albedo$age_and_season, split = " ")
dd_albedo<-data.frame(matrix(unlist(d_albedo), ncol=2, byrow=TRUE))
albedo$age<-as.numeric(dd_albedo$X1)
albedo$season<-dd_albedo$X2

albedo$season<-as.numeric(albedo$season)
albedo<-albedo[,c("lon","lat","age","season","albedo")]
albedo_cast=dcast(albedo,lon+lat+age~season, value.var = "albedo") 

# Get land part
if(!require(maptools)){ install.packages("maptools");library(maptools)}
data(wrld_simpl)
set.seed(0)#Create a SpatialPoints object
points<-unique(albedo_cast[,c("lon","lat")])
wrld_simpl@proj4string@projargs<-"+proj=longlat +datum=WGS84 +no_defs"
pts <- SpatialPoints(points, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
wrld<-over(pts, wrld_simpl)# Find which points fall over land except Antarctica and Greenland
wrld<-cbind.data.frame(points,wrld)

#get icefree land part
albedo_cast_land<-merge(albedo_cast,wrld,by=c("lon","lat"))
albedo_cast_land<-albedo_cast_land[!is.na(albedo_cast_land$FIPS),]#get land points
albedo_cast_land<-albedo_cast_land[(albedo_cast_land$NAME!="Greenland")&(albedo_cast_land$NAME!="Antarctica"),]#exclude Greenland and Antarctica
albedo_cast_land<-albedo_cast_land[,colnames(albedo_cast)]

#Adjuste age and bin
LOVECLIM_seasonal_albedo<-albedo_cast_land
LOVECLIM_seasonal_albedo$age<-b0_t+b1_t*LOVECLIM_seasonal_albedo$age #adjust the age scale

LOVECLIM_seasonal_albedo$age<-ceiling(LOVECLIM_seasonal_albedo$age/25)*25 #bin in 25 years
LOVECLIM_seasonal_albedo$lon<-ceiling(LOVECLIM_seasonal_albedo$lon/5.625)*5.625 #bin in 5.625 degree
LOVECLIM_seasonal_albedo$lat<-ceiling(LOVECLIM_seasonal_albedo$lat/5.625)*5.625 #bin in 5.625 degree
mean_snow<-aggregate(data=LOVECLIM_seasonal_albedo,.~lon+lat+age,FUN=mean, na.action = na.omit) #get the mean in each bin
err_snow<-aggregate(data=LOVECLIM_seasonal_albedo,.~lon+lat+age,FUN=sd, na.action = na.omit) #get the mean in each bin

identical(mean_snow[,c("lon","lat","age")],err_snow[,c("lon","lat","age")])
LOVECLIM_seasonal_albedo_binned<-cbind.data.frame(mean_snow,err_snow[,-c(1:3)])
colnames(LOVECLIM_seasonal_albedo_binned)<-c("lon","lat","age","albedo_1","albedo_2","albedo_3","albedo_4",
                                             "err_albedo_1","err_albedo_2","err_albedo_3","err_albedo_4")

write.csv(LOVECLIM_seasonal_albedo_binned,paste(wd,"Output data/LOVECLIM/LOVECLIM_seasonal_albedo_binned_land.csv",sep=""))



##################################################################################################
##################################################################################################
#######################                                                     ######################
#######################   4. LOVECLIM snow depth simulations (seasonal)     ######################
#######################                                                     ######################
##################################################################################################
##################################################################################################

gc()
options(future.globals.maxSize = 6 * 1024^8)

############################ Load snow_depth

snow_depth1<-brick(paste(wd,"Input data/LOVECLIM/LOVECLIM_seasonal snow depth_50-416.nc",sep=""))
snow_depth1 <- rasterToPoints(snow_depth1, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
snow_depth1 <- spTransform(snow_depth1, CRS(snow_depth1@proj4string@projargs))# reproject sp object
snow_depth1.all <- data.frame(lon=coordinates(snow_depth1)[,1],lat=coordinates(snow_depth1)[,2],snow_depth1@data)    # Assign coordinates to @data slot                     
age1<-rep(c(50000:(41600+1)),each=4)
age_miss<-seq(from=28800,to=50000,by=200)
age1<-age1[!(age1%in%age_miss)]
season<-rep(c(1:4),times=length(age1)/4)
age1<-paste(age1,season)
colnames(snow_depth1.all)<-c("lon",'lat',age1)
snow_depth1.all$lon<-ifelse(snow_depth1.all$lon<=180,snow_depth1.all$lon,snow_depth1.all$lon-360)

snow_depth2<-brick(paste(wd,"Input data/LOVECLIM/LOVECLIM_seasonal snow depth_417-30.nc",sep=""))
snow_depth2 <- rasterToPoints(snow_depth2, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
snow_depth2 <- spTransform(snow_depth2, CRS(snow_depth2@proj4string@projargs))# reproject sp object
snow_depth2.all <- data.frame(lon=coordinates(snow_depth2)[,1],lat=coordinates(snow_depth2)[,2],snow_depth2@data)    # Assign coordinates to @data slot                     
age2<-rep(c(41700:(28900+1)),each=4)
age_miss<-seq(from=28800,to=50000,by=200)
age2<-age2[!(age2%in%age_miss)]
season<-rep(c(1:4),times=length(age2)/4)
age2<-paste(age2,season)
colnames(snow_depth2.all)<-c("lon","lat",age2)
snow_depth2.all$lon<-ifelse(snow_depth2.all$lon<=180,snow_depth2.all$lon,snow_depth2.all$lon-360)

replicate<-rep(c(41700:41601),each=4)
replicate<-paste(replicate,rep(c(1:4),times=length(replicate)/4))
snow_depth<-cbind.data.frame(snow_depth1.all[,-which(colnames(snow_depth1.all)%in%replicate)],snow_depth2.all[,-c(1,2)])
snow_depth <- melt(snow_depth,id.vars = c("lon","lat"))
colnames(snow_depth)<-c("lon","lat","age_and_season","snow_depth")
snow_depth$age_and_season<-as.character(snow_depth$age_and_season)
d_snow_depth<-strsplit(snow_depth$age_and_season, split = " ")
dd_snow_depth<-data.frame(matrix(unlist(d_snow_depth), ncol=2, byrow=TRUE))
snow_depth$age<-as.numeric(dd_snow_depth$X1)
snow_depth$season<-dd_snow_depth$X2

snow_depth$season<-as.numeric(snow_depth$season)
snow_depth<-snow_depth[,c("lon","lat","age","season","snow_depth")]
snow_depth_cast=dcast(snow_depth,lon+lat+age~season, value.var = "snow_depth") 

# Get land part
if(!require(maptools)){ install.packages("maptools");library(maptools)}
data(wrld_simpl)
set.seed(0)#Create a SpatialPoints object
points<-unique(snow_depth_cast[,c("lon","lat")])
wrld_simpl@proj4string@projargs<-"+proj=longlat +datum=WGS84 +no_defs"
pts <- SpatialPoints(points, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
wrld<-over(pts, wrld_simpl)# Find which points fall over land except Antarctica and Greenland
wrld<-cbind.data.frame(points,wrld)

#get icefree land part
snow_depth_cast_land<-merge(snow_depth_cast,wrld,by=c("lon","lat"))
snow_depth_cast_land<-snow_depth_cast_land[!is.na(snow_depth_cast_land$FIPS),]#get land points
snow_depth_cast_land<-snow_depth_cast_land[(snow_depth_cast_land$NAME!="Greenland")&(snow_depth_cast_land$NAME!="Antarctica"),]#exclude Greenland and Antarctica
snow_depth_cast_land<-snow_depth_cast_land[,colnames(snow_depth_cast)]

#Adjuste age and bin
LOVECLIM_seasonal_snow_depth<-snow_depth_cast_land
LOVECLIM_seasonal_snow_depth$age<-b0_t+b1_t*LOVECLIM_seasonal_snow_depth$age #adjust the age scale

LOVECLIM_seasonal_snow_depth$age<-ceiling(LOVECLIM_seasonal_snow_depth$age/25)*25 #bin in 25 years
LOVECLIM_seasonal_snow_depth$lon<-ceiling(LOVECLIM_seasonal_snow_depth$lon/5.625)*5.625 #bin in 5.625 degree
LOVECLIM_seasonal_snow_depth$lat<-ceiling(LOVECLIM_seasonal_snow_depth$lat/5.625)*5.625 #bin in 5.625 degree
mean_snow<-aggregate(data=LOVECLIM_seasonal_snow_depth,.~lon+lat+age,FUN=mean, na.action = na.omit) #get the mean in each bin
err_snow<-aggregate(data=LOVECLIM_seasonal_snow_depth,.~lon+lat+age,FUN=sd, na.action = na.omit) #get the mean in each bin

identical(mean_snow[,c("lon","lat","age")],err_snow[,c("lon","lat","age")])
LOVECLIM_seasonal_snow_depth_binned<-cbind.data.frame(mean_snow,err_snow[,-c(1:3)])
colnames(LOVECLIM_seasonal_snow_depth_binned)<-c("lon","lat","age","snow_depth_1","snow_depth_2","snow_depth_3","snow_depth_4",
                                                 "err_snow_depth_1","err_snow_depth_2","err_snow_depth_3","err_snow_depth_4")

write.csv(LOVECLIM_seasonal_snow_depth_binned,paste(wd,"Output data/LOVECLIM/LOVECLIM_seasonal_snow_depth_binned_land.csv",sep=""))

##################################################################################################
##################################################################################################
####################                                                             #################
####################   5. Convert LOVECLIM snow depth to LOVECLIM snow cover     #################
####################                                                             #################
##################################################################################################
##################################################################################################

############################ modern_snow_depth

modern_snow_depth_raster<-brick(paste(wd,"Input data/Decomposition of albedo/monthly snow depth water equivalent 2005 0.1Deg.nc",sep=""))
modern_snow_depth_raster <- rasterToPoints(modern_snow_depth_raster, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
modern_snow_depth_raster <- spTransform(modern_snow_depth_raster, CRS(modern_snow_depth_raster@proj4string@projargs))# reproject sp object
modern_snow_depth.all <- data.frame(modern_snow_depth_raster@data, lon=coordinates(modern_snow_depth_raster)[,1],lat=coordinates(modern_snow_depth_raster)[,2])    # Assign coordinates to @data slot                     
colnames(modern_snow_depth.all)<-c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12","lon","lat")
modern_snow_depth.all$lon<-ifelse(modern_snow_depth.all$lon<=180,modern_snow_depth.all$lon,modern_snow_depth.all$lon-360)

modern_snow_depth_monthly<-modern_snow_depth.all
modern_snow_depth_monthly$lon<-ceiling(modern_snow_depth_monthly$lon/5.625) *5.625
modern_snow_depth_monthly$lat<-ceiling(modern_snow_depth_monthly$lat/5.625) *5.625
modern_snow_depth_monthly<-aggregate(. ~ lon + lat, data = modern_snow_depth_monthly, FUN = mean)

modern_snow_depth_seasonal<-cbind.data.frame(lon=modern_snow_depth_monthly[,"lon"],
                                             lat=modern_snow_depth_monthly[,"lat"],
                                             season1=rowSums(modern_snow_depth_monthly[,c("m12","m1","m2")])/3,
                                             season2=rowSums(modern_snow_depth_monthly[,c("m3","m4","m5")])/3,
                                             season3=rowSums(modern_snow_depth_monthly[,c("m6","m7","m8")])/3,
                                             season4=rowSums(modern_snow_depth_monthly[,c("m9","m10","m11")])/3)
modern_snow_depth_seasonal<-melt(modern_snow_depth_seasonal,id.vars = c("lon","lat"))
colnames(modern_snow_depth_seasonal)=c("lon","lat","season","snow_depth")

############################ modern_snow_cover

modern_snow_cover_raster<-brick(paste(wd,"Input data/Decomposition of albedo/monthly snow cover 2005 0.1Deg.nc",sep=""))
modern_snow_cover_raster <- rasterToPoints(modern_snow_cover_raster, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
modern_snow_cover_raster <- spTransform(modern_snow_cover_raster, CRS(modern_snow_cover_raster@proj4string@projargs))# reproject sp object
modern_snow_cover.all <- data.frame(modern_snow_cover_raster@data, lon=coordinates(modern_snow_cover_raster)[,1],lat=coordinates(modern_snow_cover_raster)[,2])    # Assign coordinates to @data slot                     
colnames(modern_snow_cover.all)<-c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12","lon","lat")
modern_snow_cover.all$lon<-ifelse(modern_snow_cover.all$lon<=180,modern_snow_cover.all$lon,modern_snow_cover.all$lon-360)
modern_snow_cover.all[,c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12")]<-modern_snow_cover.all[,c("m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12")]/100

modern_snow_cover_monthly<-modern_snow_cover.all
modern_snow_cover_monthly$lon<-ceiling(modern_snow_cover_monthly$lon/5.625) *5.625
modern_snow_cover_monthly$lat<-ceiling(modern_snow_cover_monthly$lat/5.625) *5.625
modern_snow_cover_monthly<-aggregate(. ~ lon + lat, data = modern_snow_cover_monthly, FUN = mean)

modern_snow_cover_seasonal<-cbind.data.frame(lon=modern_snow_cover_monthly[,"lon"],
                                             lat=modern_snow_cover_monthly[,"lat"],
                                             season1=rowSums(modern_snow_cover_monthly[,c("m12","m1","m2")])/3,
                                             season2=rowSums(modern_snow_cover_monthly[,c("m3","m4","m5")])/3,
                                             season3=rowSums(modern_snow_cover_monthly[,c("m6","m7","m8")])/3,
                                             season4=rowSums(modern_snow_cover_monthly[,c("m9","m10","m11")])/3)
modern_snow_cover_seasonal<-melt(modern_snow_cover_seasonal,id.vars = c("lon","lat"))
colnames(modern_snow_cover_seasonal)=c("lon","lat","season","snow_cover")

###################  Get the snow conversion relationship using modern data

identical(modern_snow_depth_seasonal[,c("lon","lat","season")],modern_snow_cover_seasonal[,c("lon","lat","season")])
df_modern_snow<-cbind.data.frame(modern_snow_depth_seasonal,snow_cover=modern_snow_cover_seasonal[,"snow_cover"]) 

fsnow <- function(snow_depth, b0,b1) {
  b0+(1-b0)*(1-exp(b1*snow_depth))
}
fit <- nls(snow_cover ~ fsnow(snow_depth,b0,b1), data=df_modern_snow,start = list(b0=0.5,b1=0.5))
summary(fit)
parameter_info_snow_conversion=summary(fit)[["coefficients"]][,c("Estimate","Std. Error")]
write.csv(parameter_info_snow_conversion,paste(wd,"Output data/Decomposition of albedo/parameter_info_snow_conversion_seasonal.csv",sep=""))

#plot
p<-ggplot(data=df_modern_snow,aes(snow_depth,snow_cover))+theme_bw()+geom_point(size=0.8)+ 
  stat_function(fun=function(snow_depth) parameter_info_snow_conversion["b0","Estimate"]+
                  (1-parameter_info_snow_conversion["b0","Estimate"])*
                  (1-exp(parameter_info_snow_conversion["b1","Estimate"]*snow_depth)),
                col="red")+labs(x="snow depth water equivalent (m)",y="snow cover")

ggsave(file=paste(wd,"Output data/Decomposition of albedo/Snow conversion between depth and cover.png",sep=""),
       p,width=7,height=4)

########################### Convert LOVECLIM snow depth into snow cover

LOVECLIM_seasonal_snow_cover_binned_land<-read.csv(paste(wd,"Output data/LOVECLIM/LOVECLIM_seasonal_snow_depth_binned_land.csv",sep=""),row.names=1)
para_snow<-read.csv(paste(wd,"Output data/Decomposition of albedo/parameter_info_snow_conversion_seasonal.csv",sep=""),row.names=1)

if(!require(Deriv)){ install.packages("Deriv");library(Deriv)}
fcover<-function(snow_depth,b0,b1){
  b0+(1-b0)*(1-exp(b1*snow_depth))
}
err_fcover<-function(snow_depth,b0,b1,err_snow_depth,err_b0,err_b1){
  D_snow_depth=Deriv(fcover,"snow_depth")
  D_b0=Deriv(fcover,"b0")
  D_b1=Deriv(fcover,"b1")
  sqrt( D_snow_depth(snow_depth,b0,b1)^2*err_snow_depth^2 +
          D_b0(snow_depth,b0,b1)^2*err_b0^2 +
          D_b1(snow_depth,b0,b1)^2*err_b1^2)
}

for(i in 1:4){
  
  snow_depth<-LOVECLIM_seasonal_snow_cover_binned_land[,paste("snow_depth_",i,sep="")]
  err_snow_depth<-LOVECLIM_seasonal_snow_cover_binned_land[,paste("err_snow_depth_",i,sep="")]
  
  LOVECLIM_seasonal_snow_cover_binned_land[,paste("snow_cover_",i,sep="")]<-fcover(snow_depth=snow_depth,
                                                                                   b0=para_snow["b0","Estimate"],
                                                                                   b1=para_snow["b1","Estimate"])
  LOVECLIM_seasonal_snow_cover_binned_land[,paste("err_snow_cover_",i,sep="")]<-err_fcover(snow_depth=snow_depth,
                                                                                           b0=para_snow["b0","Estimate"],
                                                                                           b1=para_snow["b1","Estimate"],
                                                                                           err_snow_depth=err_snow_depth,
                                                                                           err_b0=para_snow["b0","Std..Error"],
                                                                                           err_b1=para_snow["b1","Std..Error"])
  
}
LOVECLIM_seasonal_snow_cover_binned_land<-LOVECLIM_seasonal_snow_cover_binned_land[,c("lon","lat","age",paste("snow_cover_",1:4,sep=""),paste("err_snow_cover_",1:4,sep="") )]
write.csv(LOVECLIM_seasonal_snow_cover_binned_land,paste(wd,"Output data/LOVECLIM/LOVECLIM_seasonal_snow_cover_binned_land.csv",sep=""))



