########################################################################################
########################################################################################
#######################                                           ######################
#######################       1. Prepare the LOVECLIM prior       ######################
#######################                                           ######################
########################################################################################
########################################################################################

rm(list=ls())

wd<-"D:/PhD Project/Feedback paper_Albedo/Data and codes/"
if(!require(fxTWAPLS)){ install.packages("fxTWAPLS");library(fxTWAPLS)}
if(!require(Deriv)){ install.packages("Deriv");library(Deriv)}
source("D:/PhD Project/Feedback paper_Albedo/Data and codes/Feedback paper_Albedo_Functions.R", encoding = 'UTF-8')

###################

#Load albedo parameter
para_albedo<-read.csv(paste(wd,"Output data/Decomposition of albedo/parameter_info_albedo.csv",sep=""),row.names=1)
ks=para_albedo["ks","Estimate"]
err_ks=para_albedo["ks","Std..Error"]

#Load data
LOVECLIM_seasonal_snow_cover_binned_land<-read.csv(paste(wd,"Output data/LOVECLIM/LOVECLIM_seasonal_snow_cover_binned_land.csv",sep=""),row.names=1)
LOVECLIM_seasonal_albedo_binned_land<-read.csv(paste(wd,"Output data/LOVECLIM/LOVECLIM_seasonal_albedo_binned_land.csv",sep=""),row.names=1)
prior<-merge(LOVECLIM_seasonal_albedo_binned_land,LOVECLIM_seasonal_snow_cover_binned_land,by=c("lon","lat","age"))

#Use mean ssr of 50 -30 ka for each grid as the weight
LOVECLIM_seasonal_ssr<-read.csv(paste(wd,"Output data/LOVECLIM/LOVECLIM_seasonal_ssr_binned_land.csv",sep=""),row.names=1)
weight<-aggregate(data=LOVECLIM_seasonal_ssr[,c("lon","lat","ssr_1","ssr_2","ssr_3","ssr_4")],.~lon+lat,FUN=function(x) mean(x), na.action = na.omit)

#################### Decompose data

#get prior_whole
prior_whole<-prior
prior_whole<-get_ssr_weighted_albedo(data=prior_whole[,c("lon","lat","age",paste("albedo_",1:4,sep=""),paste("err_albedo_",1:4,sep=""))],weight)

#get prior_veg_constant
prior_veg_constant<-prior[,c("lon","lat","age")]
prior_veg_constant[,c("albedo_1","err_albedo_1")]<-hold_veg_constant(prior[,c("lon","lat","age","albedo_1","err_albedo_1","snow_cover_1","err_snow_cover_1")])
prior_veg_constant[,c("albedo_2","err_albedo_2")]<-hold_veg_constant(prior[,c("lon","lat","age","albedo_2","err_albedo_2","snow_cover_2","err_snow_cover_2")])
prior_veg_constant[,c("albedo_3","err_albedo_3")]<-hold_veg_constant(prior[,c("lon","lat","age","albedo_3","err_albedo_3","snow_cover_3","err_snow_cover_3")])
prior_veg_constant[,c("albedo_4","err_albedo_4")]<-hold_veg_constant(prior[,c("lon","lat","age","albedo_4","err_albedo_4","snow_cover_4","err_snow_cover_4")])
prior_veg_constant<-get_ssr_weighted_albedo(data=prior_veg_constant[,c("lon","lat","age",paste("albedo_",1:4,sep=""),paste("err_albedo_",1:4,sep=""))],weight)

#get prior_snow_constant
prior_snow_constant<-prior[,c("lon","lat","age")]
prior_snow_constant[,c("albedo_1","err_albedo_1")]<-hold_snow_constant(prior[,c("lon","lat","age","albedo_1","err_albedo_1","snow_cover_1","err_snow_cover_1")])
prior_snow_constant[,c("albedo_2","err_albedo_2")]<-hold_snow_constant(prior[,c("lon","lat","age","albedo_2","err_albedo_2","snow_cover_2","err_snow_cover_2")])
prior_snow_constant[,c("albedo_3","err_albedo_3")]<-hold_snow_constant(prior[,c("lon","lat","age","albedo_3","err_albedo_3","snow_cover_3","err_snow_cover_3")])
prior_snow_constant[,c("albedo_4","err_albedo_4")]<-hold_snow_constant(prior[,c("lon","lat","age","albedo_4","err_albedo_4","snow_cover_4","err_snow_cover_4")])
prior_snow_constant<-get_ssr_weighted_albedo(data=prior_snow_constant[,c("lon","lat","age",paste("albedo_",1:4,sep=""),paste("err_albedo_",1:4,sep=""))],weight)

#get prior_both_constant
mean_albedo<-aggregate(data=prior[,c("lon","lat","albedo_1","albedo_2","albedo_3","albedo_4")],.~lon+lat,FUN=function(x) mean(x), na.action = na.omit)
err_mean_albedo<-aggregate(data=prior[,c("lon","lat","err_albedo_1","err_albedo_2","err_albedo_3","err_albedo_4")],.~lon+lat,FUN=function(x) sqrt(sum(x^2))/length(x), na.action = na.omit)
prior_both_constant<-merge(mean_albedo,err_mean_albedo,by=c("lon","lat"))
prior_both_constant<-merge(prior[,c("lon","lat","age")],prior_both_constant,by=c("lon","lat"))
prior_both_constant<-get_ssr_weighted_albedo(data=prior_both_constant[,c("lon","lat","age",paste("albedo_",1:4,sep=""),paste("err_albedo_",1:4,sep=""))],weight)

############### Combine

identical(prior[,c("lon","lat","age")],prior_whole[,c("lon","lat","age")])
identical(prior[,c("lon","lat","age")],prior_veg_constant[,c("lon","lat","age")])
identical(prior[,c("lon","lat","age")],prior_snow_constant[,c("lon","lat","age")])
identical(prior[,c("lon","lat","age")],prior_both_constant[,c("lon","lat","age")])

prior_decomposition<-cbind.data.frame(prior[,c("lon","lat","age")],
                                      prior_whole[,c("albedo","err_albedo")],
                                      prior_veg_constant[,c("albedo","err_albedo")],
                                      prior_snow_constant[,c("albedo","err_albedo")],
                                      prior_both_constant[,c("albedo","err_albedo")])

colnames(prior_decomposition)=c("lon","lat","age",
                                "albedo","err_albedo",
                                "albedo_veg_constant","err_albedo_veg_constant",
                                "albedo_snow_constant","err_albedo_snow_constant",
                                "albedo_both_constant","err_albedo_both_constant")

write.csv(prior_decomposition,paste(wd,"Output data/3D var/3D var input/LOVECLIM_albedo_decomposition.csv",sep=""))


########################################################################################
########################################################################################
#######################                                           ######################
#######################     2. Prepare the observation input      ######################
#######################                                           ######################
########################################################################################
########################################################################################

rm(list=ls())

wd<-"D:/PhD Project/Feedback paper_Albedo/Data and codes/"
if(!require(fxTWAPLS)){ install.packages("fxTWAPLS");library(fxTWAPLS)}
if(!require(Deriv)){ install.packages("Deriv");library(Deriv)}
source("D:/PhD Project/Feedback paper_Albedo/Data and codes/Feedback paper_Albedo_Functions.R", encoding = 'UTF-8')

###################

recon_adjusted<-read.csv(paste(wd,"Output data/ACER reconstructions/Reconstruction/recon_adjusted.csv",sep=""),row.names=1)
LOVECLIM_seasonal_snow_cover<-read.csv(paste(wd,"Output data/LOVECLIM/LOVECLIM_seasonal_snow_cover_binned_land.csv",sep=""),row.names=1)
recon<-merge(recon_adjusted,LOVECLIM_seasonal_snow_cover,by=c("lon","lat","age"))

#Load parameters of the albedo decomposition model
para_albedo<-read.csv(paste(wd,"Output data/Decomposition of albedo/parameter_info_albedo.csv",sep=""),row.names=1)

#Use mean ssr of 50 -30 ka for each grid as the weight
LOVECLIM_seasonal_ssr<-read.csv(paste(wd,"Output data/LOVECLIM/LOVECLIM_seasonal_ssr_binned_land.csv",sep=""),row.names=1)
weight<-aggregate(data=LOVECLIM_seasonal_ssr[,c("lon","lat","ssr_1","ssr_2","ssr_3","ssr_4")],.~lon+lat,FUN=function(x) mean(x), na.action = na.omit)

#################### Decompose data

#whole
recon_whole<-recon
recon_whole<-compose_albedo(recon_whole,para_albedo,weight)
recon_whole<-get_ssr_weighted_albedo(data=recon_whole[,c("lon","lat","age",paste("albedo_",1:4,sep=""),paste("err_albedo_",1:4,sep=""))],weight)

#hold vegetation constant
mean_veg<-aggregate(data=recon[,c("lon","lat","fapar","hveg")],.~lon+lat,FUN=function(x) mean(x), na.action = na.omit)
err_veg<-aggregate(data=recon[,c("lon","lat","sse_fapar","sse_hveg")],.~lon+lat,FUN=function(x) sqrt(sum(x^2))/length(x), na.action = na.omit)
grid_veg<-merge(x=mean_veg,y=err_veg,by=c("lon","lat")) #use the mean over 50~30ka for each grid
recon_veg_constant<-merge(x=recon[,c("lon","lat","age","site_id",paste("snow_cover_",1:4,sep=""),paste("err_snow_cover_",1:4,sep=""))],
                          y=grid_veg,by=c("lon","lat"))
recon_veg_constant<-compose_albedo(recon_veg_constant,para_albedo,weight)
recon_veg_constant<-get_ssr_weighted_albedo(data=recon_veg_constant[,c("lon","lat","age",paste("albedo_",1:4,sep=""),paste("err_albedo_",1:4,sep=""))],weight)

#hold snow constant
mean_snow<-aggregate(data=recon[,c("lon","lat",paste("snow_cover_",1:4,sep=""))],.~lon+lat,FUN=function(x) mean(x), na.action = na.omit)
err_snow<-aggregate(data=recon[,c("lon","lat",paste("err_snow_cover_",1:4,sep=""))],.~lon+lat,FUN=function(x) sqrt(sum(x^2))/length(x), na.action = na.omit)
grid_snow<-merge(x=mean_snow,y=err_snow,by=c("lon","lat")) #use the mean over 50~30ka for each grid
recon_snow_constant<-merge(x=recon[,c("lon","lat","age","site_id","fapar","hveg","sse_fapar","sse_hveg")],
                           y=grid_snow,by=c("lon","lat"))
recon_snow_constant<-compose_albedo(recon_snow_constant,para_albedo,weight)
recon_snow_constant<-get_ssr_weighted_albedo(data=recon_snow_constant[,c("lon","lat","age",paste("albedo_",1:4,sep=""),paste("err_albedo_",1:4,sep=""))],weight)

#hold both constant
mean_both<-aggregate(data=recon[,c("lon","lat","fapar","hveg",paste("snow_cover_",1:4,sep=""))],.~lon+lat,FUN=function(x) mean(x), na.action = na.omit)
err_both<-aggregate(data=recon[,c("lon","lat","sse_fapar","sse_hveg",paste("err_snow_cover_",1:4,sep=""))],.~lon+lat,FUN=function(x) sqrt(sum(x^2))/length(x), na.action = na.omit)
grid_both<-merge(x=mean_both,y=err_both,by=c("lon","lat")) #use the mean over 50~30ka for each grid
recon_both_constant<-merge(x=recon[,c("lon","lat","age","site_id")],
                           y=grid_both,by=c("lon","lat"))
recon_both_constant<-compose_albedo(recon_both_constant,para_albedo,weight)
recon_both_constant<-get_ssr_weighted_albedo(data=recon_both_constant[,c("lon","lat","age",paste("albedo_",1:4,sep=""),paste("err_albedo_",1:4,sep=""))],weight)

##################### Combine

identical(recon[,c("lon","lat","age")],recon_whole[,c("lon","lat","age")])
identical(recon[,c("lon","lat","age")],recon_veg_constant[,c("lon","lat","age")])
identical(recon[,c("lon","lat","age")],recon_snow_constant[,c("lon","lat","age")])
identical(recon[,c("lon","lat","age")],recon_both_constant[,c("lon","lat","age")])

recon_decomposition<-cbind.data.frame(recon[,c("lon","lat","age","site_id","Tmean","sse_Tmean")],
                                      recon_whole[,c("albedo","err_albedo")],
                                      recon_veg_constant[,c("albedo","err_albedo")],
                                      recon_snow_constant[,c("albedo","err_albedo")],
                                      recon_both_constant[,c("albedo","err_albedo")])

colnames(recon_decomposition)=c("lon","lat","age","site_id","Tmean","sse_Tmean",
                                "albedo","err_albedo",
                                "albedo_veg_constant","err_albedo_veg_constant",
                                "albedo_snow_constant","err_albedo_snow_constant",
                                "albedo_both_constant","err_albedo_both_constant")
write.csv(recon_decomposition,paste(wd,"Output data/3D var/3D var input/recon_decomposition.csv",sep=""))


########################################################################################
########################################################################################
#######################                                           ######################
#######################     3. Apply 3D variational technique     ######################
#######################                                           ######################
########################################################################################
########################################################################################

rm(list=ls())

wd<-"D:/PhD Project/Feedback paper_Albedo/Data and codes/"
if(!require(fxTWAPLS)){ install.packages("fxTWAPLS");library(fxTWAPLS)}

source("D:/PhD Project/Feedback paper_Albedo/Data and codes/Feedback paper_Albedo_Functions.R", encoding = 'UTF-8')
source("D:/PhD Project/Feedback paper_Albedo/Data and codes/Feedback paper_Albedo_Functions_3D variational technique.R", encoding = 'UTF-8')

###########################################################
#################  Get the Ls and Lt  
###########################################################

DO_obs<-read.csv(paste(wd,"Output data/3D var/3D var input/recon_decomposition.csv",sep=""),row.names=1)

#Get Ls
site<-unique(DO_obs[,c("lon","lat")])
dist<-fxTWAPLS::get_distance(site)
mean(dist) 
Ls=mean(dist)/1000

#Get Lt
reso<-get_average_resolution(DO_obs[,c("site_id","age")])
Lt=mean(reso$avg_reso)

#We recommend run the 3D var code on HPC
interval=25

##################################################################################################
###############################   3D var temperature
##################################################################################################

DO_model<-read.csv(paste(wd,"Output data/3D var/3D var input/LOVECLIM_T_binned_land.csv",sep=""),row.names=1)
DO_obs<-read.csv(paste(wd,"Output data/3D var/3D var input/recon_decomposition.csv",sep=""),row.names=1)

DO_model<-DO_model[order(DO_model$age),]
DO_model<-na.omit(DO_model)

DO_obs<-DO_obs[,c("site_id","lon","lat","age","Tmean","sse_Tmean")]
DO_obs<-na.omit(DO_obs)

output_all_age<-data.frame()

for(t in unique(DO_model$age)){
  
  model<-DO_model[which(DO_model$age>=(t-interval) &DO_model$age<=(t+interval)),]
  obs<-DO_obs[which(DO_obs$age>=(t-interval) &DO_obs$age<=(t+interval)),]
  
  #extract useful variables
  model<-model[,c("lon","lat","age","T","sd_T")]
  obs<-obs[,c("lon","lat","age","Tmean","sse_Tmean")]
  colnames(model)=c("lon","lat","age","value","err")
  colnames(obs)=c("lon","lat","age","value","err")
  
  match<-merge(x=model,y=obs,by=c("lon","lat","age"))
  obs<-match[,c("lon","lat","age","value.y","err.y")]
  colnames(obs)=c("lon","lat","age","value","err")
  
  if(nrow(obs)>1){
    
    print(t)
    model$X<-seq(1:nrow(model))
    obs$X<-seq(1:nrow(obs))
    #plot(match$value.y~match$value.x);abline(a=0,b=1)
    
    output = get_xa(model=model,obs=obs,Ls=Ls,Lt=Lt)
    output_all_age<-rbind.data.frame(output_all_age,output)
  }
}

write.csv(output_all_age, paste(wd,"Output data/3D var/3D var output/3Dvar temperature.csv",sep=""))

##################################################################################################
###############################   3D var albedo_whole model
##################################################################################################

DO_model<-read.csv(paste(wd,"Output data/3D var/3D var input/LOVECLIM_albedo_decomposition.csv",sep=""),row.names=1)
DO_obs<-read.csv(paste(wd,"Output data/3D var/3D var input/recon_decomposition.csv",sep=""),row.names=1)

DO_obs<-na.omit(DO_obs)
DO_model<-na.omit(DO_model)

DO_obs<-DO_obs[order(DO_obs$age),]
DO_model<-DO_model[order(DO_model$age),]

output_all_age<-data.frame()

for(t in unique(DO_model$age)){
  
  model<-DO_model[which(DO_model$age>=(t-interval) &DO_model$age<=(t+interval)),]
  obs<-DO_obs[which(DO_obs$age>=(t-interval) &DO_obs$age<=(t+interval)),]
  
  #extract useful variables
  model<-model[,c("lon","lat","age","albedo","err_albedo")]
  obs<-obs[,c("lon","lat","age","albedo","err_albedo")]
  colnames(model)=c("lon","lat","age","value","err")
  colnames(obs)=c("lon","lat","age","value","err")
  
  match<-merge(x=model,y=obs,by=c("lon","lat","age"))
  obs<-match[,c("lon","lat","age","value.y","err.y")]
  colnames(obs)=c("lon","lat","age","value","err")
  
  if(nrow(obs)>1){
    
    print(t)
    model$X<-seq(1:nrow(model))
    obs$X<-seq(1:nrow(obs))
    #plot(match$value.y~match$value.x);abline(a=0,b=1)
    
    output = get_xa(model=model,obs=obs,Ls=Ls,Lt=Lt)
    output_all_age<-rbind.data.frame(output_all_age,output)
  }
}

write.csv(output_all_age, paste(wd,"Output data/3D var/3D var output/3Dvar albedo_whole model.csv",sep=""))


##################################################################################################
###############################   3D var albedo_veg constant model
##################################################################################################

DO_model<-read.csv(paste(wd,"Output data/3D var/3D var input/LOVECLIM_albedo_decomposition.csv",sep=""),row.names=1)
DO_obs<-read.csv(paste(wd,"Output data/3D var/3D var input/recon_decomposition.csv",sep=""),row.names=1)

DO_obs<-na.omit(DO_obs)
DO_model<-na.omit(DO_model)

DO_obs<-DO_obs[order(DO_obs$age),]
DO_model<-DO_model[order(DO_model$age),]

output_all_age<-data.frame()

for(t in unique(DO_model$age)){
  
  model<-DO_model[which(DO_model$age>=(t-interval) &DO_model$age<=(t+interval)),]
  obs<-DO_obs[which(DO_obs$age>=(t-interval) &DO_obs$age<=(t+interval)),]
  
  #extract useful variables
  model<-model[,c("lon","lat","age","albedo_veg_constant","err_albedo_veg_constant")]
  obs<-obs[,c("lon","lat","age","albedo_veg_constant","err_albedo_veg_constant")]
  colnames(model)=c("lon","lat","age","value","err")
  colnames(obs)=c("lon","lat","age","value","err")
  
  match<-merge(x=model,y=obs,by=c("lon","lat","age"))
  obs<-match[,c("lon","lat","age","value.y","err.y")]
  colnames(obs)=c("lon","lat","age","value","err")
  
  if(nrow(obs)>1){
    
    print(t)
    model$X<-seq(1:nrow(model))
    obs$X<-seq(1:nrow(obs))
    #plot(match$value.y~match$value.x);abline(a=0,b=1)
    
    output = get_xa(model=model,obs=obs,Ls=Ls,Lt=Lt)
    output_all_age<-rbind.data.frame(output_all_age,output)
  }
}

write.csv(output_all_age, paste(wd,"Output data/3D var/3D var output/3Dvar albedo_veg constant model.csv",sep=""))

##################################################################################################
###############################   3D var albedo_snow constant model
##################################################################################################

DO_model<-read.csv(paste(wd,"Output data/3D var/3D var input/LOVECLIM_albedo_decomposition.csv",sep=""),row.names=1)
DO_obs<-read.csv(paste(wd,"Output data/3D var/3D var input/recon_decomposition.csv",sep=""),row.names=1)

DO_obs<-na.omit(DO_obs)
DO_model<-na.omit(DO_model)

DO_obs<-DO_obs[order(DO_obs$age),]
DO_model<-DO_model[order(DO_model$age),]

output_all_age<-data.frame()

for(t in unique(DO_model$age)){
  
  model<-DO_model[which(DO_model$age>=(t-interval) &DO_model$age<=(t+interval)),]
  obs<-DO_obs[which(DO_obs$age>=(t-interval) &DO_obs$age<=(t+interval)),]
  
  #extract useful variables
  model<-model[,c("lon","lat","age","albedo_snow_constant","err_albedo_snow_constant")]
  obs<-obs[,c("lon","lat","age","albedo_snow_constant","err_albedo_snow_constant")]
  colnames(model)=c("lon","lat","age","value","err")
  colnames(obs)=c("lon","lat","age","value","err")
  
  match<-merge(x=model,y=obs,by=c("lon","lat","age"))
  obs<-match[,c("lon","lat","age","value.y","err.y")]
  colnames(obs)=c("lon","lat","age","value","err")
  
  if(nrow(obs)>1){
    
    print(t)
    model$X<-seq(1:nrow(model))
    obs$X<-seq(1:nrow(obs))
    #plot(match$value.y~match$value.x);abline(a=0,b=1)
    
    output = get_xa(model=model,obs=obs,Ls=Ls,Lt=Lt)
    output_all_age<-rbind.data.frame(output_all_age,output)
  }
}

write.csv(output_all_age, paste(wd,"Output data/3D var/3D var output/3Dvar albedo_snow constant model.csv",sep=""))

##################################################################################################
###############################   3D var albedo_both constant model
##################################################################################################

DO_model<-read.csv(paste(wd,"Output data/3D var/3D var input/LOVECLIM_albedo_decomposition.csv",sep=""),row.names=1)
DO_obs<-read.csv(paste(wd,"Output data/3D var/3D var input/recon_decomposition.csv",sep=""),row.names=1)

DO_obs<-na.omit(DO_obs)
DO_model<-na.omit(DO_model)

DO_obs<-DO_obs[order(DO_obs$age),]
DO_model<-DO_model[order(DO_model$age),]

output_all_age<-data.frame()

for(t in unique(DO_model$age)){
  
  model<-DO_model[which(DO_model$age>=(t-interval) &DO_model$age<=(t+interval)),]
  obs<-DO_obs[which(DO_obs$age>=(t-interval) &DO_obs$age<=(t+interval)),]
  
  #extract useful variables
  model<-model[,c("lon","lat","age","albedo_both_constant","err_albedo_both_constant")]
  obs<-obs[,c("lon","lat","age","albedo_both_constant","err_albedo_both_constant")]
  colnames(model)=c("lon","lat","age","value","err")
  colnames(obs)=c("lon","lat","age","value","err")
  
  match<-merge(x=model,y=obs,by=c("lon","lat","age"))
  obs<-match[,c("lon","lat","age","value.y","err.y")]
  colnames(obs)=c("lon","lat","age","value","err")
  
  if(nrow(obs)>1){
    
    print(t)
    model$X<-seq(1:nrow(model))
    obs$X<-seq(1:nrow(obs))
    #plot(match$value.y~match$value.x);abline(a=0,b=1)
    
    output = get_xa(model=model,obs=obs,Ls=Ls,Lt=Lt)
    output_all_age<-rbind.data.frame(output_all_age,output)
  }
}

write.csv(output_all_age, paste(wd,"Output data/3D var/3D var output/3Dvar albedo_both constant model.csv",sep=""))





