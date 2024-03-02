rm(list=ls())

wd<-"D:/PhD Project/Feedback paper_Albedo/Data and codes/"

if(!require(egg)){ install.packages("egg");library(egg)}
if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(ggmap)){ install.packages("ggmap");library(ggmap)}
if(!require(ggsn)){ install.packages("ggsn");library(ggsn)}
if(!require(maps)){ install.packages("maps");library(maps)}
if(!require(mapdata)){ install.packages("mapdata");library(mapdata)}
world <- map_data("world") 

source("D:/PhD Project/Feedback paper_Albedo/Data and codes/Feedback paper_Albedo_Functions.R", encoding = 'UTF-8')


##############################################################################################
##############################################################################################
####################                                                 #########################
####################      1. Global mean land temperature and assr   #########################
####################                                                 #########################
##############################################################################################
##############################################################################################

land_T<-read.csv(paste(wd,"Output data/3D var/3D var output/3Dvar temperature.csv",sep=""),row.names=1)
land_whole_model<-read.csv(paste(wd,"Output data/3D var/3D var output/3Dvar albedo_whole model.csv",sep=""),row.names=1)
land_veg_constant_model<-read.csv(paste(wd,"Output data/3D var/3D var output/3Dvar albedo_veg constant model.csv",sep=""),row.names=1)
land_snow_constant_model<-read.csv(paste(wd,"Output data/3D var/3D var output/3Dvar albedo_snow constant model.csv",sep=""),row.names=1)
land_both_constant_model<-read.csv(paste(wd,"Output data/3D var/3D var output/3Dvar albedo_both constant model.csv",sep=""),row.names=1)

################################## smooth over the age interval

smooth_over_interval<-function(data){
  mean<-aggregate(data[,c("lon","lat","age","xa")],.~lon+lat+age,FUN=function(x) mean(x), na.action = na.omit) #get the mean in each bin
  err<-aggregate(data[,c("lon","lat","age","err_xa")],.~lon+lat+age,FUN=function(x) sqrt(sum(x^2))/length(x), na.action = na.omit)#get the err in each bin
  output<-merge(mean,err,by=c("lon","lat","age"))
  return(output)
}

land_T <- smooth_over_interval(land_T)
land_whole_model <- smooth_over_interval(land_whole_model)
land_veg_constant_model <- smooth_over_interval(land_veg_constant_model)
land_snow_constant_model <- smooth_over_interval(land_snow_constant_model)
land_both_constant_model <- smooth_over_interval(land_both_constant_model)

############## Get assr using LOVECLIM ssr (mean value over 50 ~ 30 ka for each grid)

LOVECLIM_ssr<-read.csv(paste(wd,"Output data/LOVECLIM/LOVECLIM_seasonal_ssr_binned_land.csv",sep=""),row.names=1)
LOVECLIM_ssr$ssr<-rowMeans(LOVECLIM_ssr[,c("ssr_1","ssr_2","ssr_3","ssr_4")])
grid_ssr<-aggregate(data=LOVECLIM_ssr[,c("lon","lat","ssr")],.~lon+lat,FUN=function(x) mean(x), na.action = na.omit)

land_whole_model<-merge(land_whole_model,grid_ssr,by=c("lon","lat"))
land_whole_model$assr<-land_whole_model$ssr*(1-land_whole_model$xa)
land_whole_model$err_assr<-land_whole_model$ssr*land_whole_model$err_xa

land_veg_constant_model<-merge(land_veg_constant_model,grid_ssr,by=c("lon","lat"))
land_veg_constant_model$assr<-land_veg_constant_model$ssr*(1-land_veg_constant_model$xa)
land_veg_constant_model$err_assr<-land_veg_constant_model$ssr*land_veg_constant_model$err_xa

land_snow_constant_model<-merge(land_snow_constant_model,grid_ssr,by=c("lon","lat"))
land_snow_constant_model$assr<-land_snow_constant_model$ssr*(1-land_snow_constant_model$xa)
land_snow_constant_model$err_assr<-land_snow_constant_model$ssr*land_snow_constant_model$err_xa

land_both_constant_model<-merge(land_both_constant_model,grid_ssr,by=c("lon","lat"))
land_both_constant_model$assr<-land_both_constant_model$ssr*(1-land_both_constant_model$xa)
land_both_constant_model$err_assr<-land_both_constant_model$ssr*land_both_constant_model$err_xa

################## Use both constant model as the base to get individual effects
get_effect<-function(data,base){
  
  data<-merge(x=data,y=base,by=c("lon","lat","age"))
  
  data$xa<-data$xa.x-data$xa.y #xa is the 3D var estimated albedo
  data$err_xa<-sqrt(data$err_xa.x^2+data$err_xa.y^2) #err_xa is the error of 3D var estimated albedo
  data$assr<-data$assr.x-data$assr.y #assr is the 3D var estimated albedo
  data$err_assr<-sqrt(data$err_assr.x^2+data$err_assr.y^2) #err_assr is the error of 3D var estimated albedo
  
  output<-data[,c("lon","lat","age","xa","err_xa","assr","err_assr")]
  return(output)
}

land_whole_effect<-get_effect(land_whole_model,land_both_constant_model)
land_veg_effect<-get_effect(land_snow_constant_model,land_both_constant_model)
land_snow_effect<-get_effect(land_veg_constant_model,land_both_constant_model)


############## Get global mean land temperature and assr

#Use the self-defined get_global_mean () function to get the cos(lat) weighted mean
#the input should have such columns: c("lon","lat","age","x","err_x")
land_T_mean_globe<-get_global_mean(land_T[,c("lon","lat","age","xa","err_xa")])
land_whole_effect_mean_globe<-get_global_mean(land_whole_effect[,c("lon","lat","age","assr","err_assr")])
land_veg_effect_mean_globe<-get_global_mean(land_veg_effect[,c("lon","lat","age","assr","err_assr")])
land_snow_effect_mean_globe<-get_global_mean(land_snow_effect[,c("lon","lat","age","assr","err_assr")])


##############################################################################################
##############################################################################################
####################                                                 #########################
####################         2. Identify minimum and maximum         #########################
####################                                                 #########################
##############################################################################################
##############################################################################################

win=2000
min_left=100
min_right=100
duration=500

#Official DO start dates
DO_timing1 <- as.matrix(read.csv(paste(wd,"Input data/Timing/DO_timing1.csv",sep=""), row.names=1, stringsAsFactors = FALSE))
DO_label<-cbind.data.frame(k=5:12,age=DO_timing1[5:12])

###################### Identify minimum and maximum

store<-data.frame()

for(k in 5:12){
  
  t_DO=DO_label[which(DO_label$k==k),"age"]
  output_T<-get_event(data=land_T_mean_globe,t_DO, win,min_left,min_right,duration)
  output_whole_effect<-get_event(data=land_whole_effect_mean_globe,t_DO, win,min_left,min_right,duration)
  output_veg_effect<-get_event(data=land_veg_effect_mean_globe,t_DO, win,min_left,min_right,duration)
  output_snow_effect<-get_event(data=land_snow_effect_mean_globe,t_DO, win,min_left,min_right,duration)
  
  output<-cbind.data.frame(k,t_DO, output_T,output_whole_effect,output_veg_effect,output_snow_effect)
  colnames(output)<-c("k","t_DO",
                      paste(colnames(output_T),"_land_T",sep=""),
                      paste(colnames(output_whole_effect),"_land_whole_effect",sep=""),
                      paste(colnames(output_veg_effect),"_land_veg_effect",sep=""),
                      paste(colnames(output_snow_effect),"_land_snow_effect",sep=""))
  
  store<-rbind.data.frame(store,output)
}

write.csv(store,paste(wd,"Output data/3D var/3Dvar_start and end.csv",sep=""))

#Get the time needed to reach maximum from minimum
round(mean(store$t_start_land_T-store$t_end_land_T),digits=0)
round(mean(store$t_start_land_veg_effect-store$t_end_land_veg_effect),digits=0)
round(mean(store$t_start_land_snow_effect-store$t_end_land_snow_effect),digits=0)

##############################################################################################
##############################################################################################
####################                                                 #########################
####################          3. Plot global mean change             #########################
####################                                                 #########################
##############################################################################################
##############################################################################################

#Clip off the edge since the two ends are affected by the ±25 year smoothing for 3D var output

clip_edge<-function(data,left_edge,right_edge){
  data[which(data$age>=left_edge & data$age<=right_edge),]
}

left_edge<-DO_timing1[5]-1000;right_edge<-DO_timing1[12]+1000

land_T_mean_globe<-clip_edge(land_T_mean_globe,left_edge,right_edge)
land_whole_effect_mean_globe<-clip_edge(land_whole_effect_mean_globe,left_edge,right_edge)
land_veg_effect_mean_globe<-clip_edge(land_veg_effect_mean_globe,left_edge,right_edge)
land_snow_effect_mean_globe<-clip_edge(land_snow_effect_mean_globe,left_edge,right_edge)

#Plot global mean changes

p_whole_effect<-ggplot()+theme_bw()+
  geom_point(data=land_whole_effect_mean_globe,aes(age,value_w),size=0.3)+
  geom_line(data=land_whole_effect_mean_globe,aes(age,value_w),size=0.3)+
  geom_ribbon(data=land_whole_effect_mean_globe,aes(x=age,y=value_w, 
                                                                ymin=value_w-1.96*err_value_w,
                                                                ymax=value_w+1.96*err_value_w),alpha=0.3)+
  labs(y= bquote("assr"["whole effect"]~ "("~ W ~ m ^-2~")"),x="Age (yr BP)")+
  geom_vline(xintercept=DO_timing1[5:12],col="red",linetype="dashed",size=0.3)+
  geom_text(data=DO_label,aes(x=age, label=k, y=Inf), hjust=0.4, vjust=-0.6,size=4)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.text=element_text(size=11),axis.title=element_text(size=12))+  
  coord_cartesian(clip = "off")+theme(plot.margin = margin(20, 20, 2, 2, unit = "pt"))+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  geom_rect(data=store,aes(xmin=t_end_land_whole_effect, xmax=t_start_land_whole_effect, 
                                ymin=min_land_whole_effect, ymax=max_land_whole_effect),fill="red",alpha=0.4)


p_veg_effect<-ggplot()+theme_bw()+
  geom_point(data=land_veg_effect_mean_globe,aes(age,value_w),size=0.3)+
  geom_line(data=land_veg_effect_mean_globe,aes(age,value_w),size=0.3)+
  geom_ribbon(data=land_veg_effect_mean_globe,aes(x=age,y=value_w, 
                                                              ymin=value_w-1.96*err_value_w,
                                                              ymax=value_w+1.96*err_value_w),alpha=0.3)+
  labs(y= bquote("assr"["vegetation effect"]~ "("~ W ~ m ^-2~")"),x="Age (yr BP)")+
  geom_vline(xintercept=DO_timing1[5:12],col="red",linetype="dashed",size=0.3)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.text=element_text(size=11),axis.title=element_text(size=12))+
  geom_rect(data=store,aes(xmin=t_end_land_veg_effect, xmax=t_start_land_veg_effect, 
                                ymin=min_land_veg_effect, ymax=max_land_veg_effect),fill="red",alpha=0.4)

p_snow_effect<-ggplot()+theme_bw()+
  geom_point(data=land_snow_effect_mean_globe,aes(age,value_w),size=0.3)+
  geom_line(data=land_snow_effect_mean_globe,aes(age,value_w),size=0.3)+
  geom_ribbon(data=land_snow_effect_mean_globe,aes(x=age,y=value_w, 
                                                               ymin=value_w-1.96*err_value_w,
                                                               ymax=value_w+1.96*err_value_w),alpha=0.3)+
  labs(y= bquote("assr"["snow effect"]~ "("~ W ~ m ^-2~")"),x="Age (yr BP)")+
  geom_vline(xintercept=DO_timing1[5:12],col="red",linetype="dashed",size=0.3)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.text=element_text(size=11),axis.title=element_text(size=12))+
  geom_rect(data=store,aes(xmin=t_end_land_snow_effect, xmax=t_start_land_snow_effect, 
                                ymin=min_land_snow_effect, ymax=max_land_snow_effect),fill="red",alpha=0.4)

p_T<-ggplot()+theme_bw()+
  geom_point(data=land_T_mean_globe,aes(age,value_w),size=0.3)+
  geom_line(data=land_T_mean_globe,aes(age,value_w),size=0.3)+
  geom_ribbon(data=land_T_mean_globe,aes(x=age,y=value_w, 
                                                ymin=value_w-1.96*err_value_w,
                                                ymax=value_w+1.96*err_value_w),alpha=0.3)+
  labs(y= bquote(T["global mean land"]~(degree~C)),x="Age (yr BP)")+
  geom_vline(xintercept=DO_timing1[5:12],col="red",linetype="dashed",size=0.3)+
  theme(axis.text=element_text(size=11),axis.title=element_text(size=12))+
  geom_rect(data=store,aes(xmin=t_end_land_T, xmax=t_start_land_T, 
                           ymin=min_land_T, ymax=max_land_T),fill="red",alpha=0.4)

p<-ggarrange(p_whole_effect,p_veg_effect,p_snow_effect,p_T, ncol=1,labels=c("(a)","(b)","(c)","(d)"),
             label.args = list(gp = grid::gpar(font = 4,cex=1)))

ggsave(file=paste(wd,"Output data/3D var/3D var figure/Fig.2 3D var global mean land changes.png",sep=""),p,width=7,height=10)


##############################################################################################
##############################################################################################
####################                                                 #########################
####################           4. Plot global change map             #########################
####################                                                 #########################
##############################################################################################
##############################################################################################

##############################################################################################
##############################    Get map data for each event
##############################################################################################

######### land_T

map_T<-data.frame()
for(k in 5:12){
  each_T<-get_change_map(k, data=land_T[,c("lon","lat","age","xa","err_xa")],
                         t_start=store[which(store$k==k),"t_start_land_T"], 
                         t_end=store[which(store$k==k),"t_end_land_T"])
  map_T<-rbind.data.frame(map_T,each_T)
}

######### land_whole_effect

map_whole_effect<-data.frame()
for(k in 5:12){
  each_whole_effect<-get_change_map(k, data=land_whole_effect[,c("lon","lat","age","xa","err_xa")],
                         t_start=store[which(store$k==k),"t_start_land_whole_effect"], 
                         t_end=store[which(store$k==k),"t_end_land_whole_effect"])
  map_whole_effect<-rbind.data.frame(map_whole_effect,each_whole_effect)
}

######### land_veg_effect

map_veg_effect<-data.frame()
for(k in 5:12){
  each_veg_effect<-get_change_map(k, data=land_veg_effect[,c("lon","lat","age","xa","err_xa")],
                                    t_start=store[which(store$k==k),"t_start_land_veg_effect"], 
                                    t_end=store[which(store$k==k),"t_end_land_veg_effect"])
  map_veg_effect<-rbind.data.frame(map_veg_effect,each_veg_effect)
}

######### land_snow_effect

map_snow_effect<-data.frame()
for(k in 5:12){
  each_snow_effect<-get_change_map(k, data=land_snow_effect[,c("lon","lat","age","xa","err_xa")],
                                    t_start=store[which(store$k==k),"t_start_land_snow_effect"], 
                                    t_end=store[which(store$k==k),"t_end_land_snow_effect"])
  map_snow_effect<-rbind.data.frame(map_snow_effect,each_snow_effect)
}

##############################################################################################
########################    Plot change map for the mean over DO 5 ~ 12
##############################################################################################

mean_map_T<-aggregate(map_T[,c("k","lon","lat","change")],.~lon+lat,FUN=function(x) mean(x))
mean_map_whole_effect<-aggregate(map_whole_effect[,c("k","lon","lat","change")],.~lon+lat,FUN=function(x) mean(x))
mean_map_veg_effect<-aggregate(map_veg_effect[,c("k","lon","lat","change")],.~lon+lat,FUN=function(x) mean(x))
mean_map_snow_effect<-aggregate(map_snow_effect[,c("k","lon","lat","change")],.~lon+lat,FUN=function(x) mean(x))


library(RColorBrewer)

#### temperature
summary(mean_map_T[,"change"])
breaks_T = c(-2,-1,-0.2, 0.2,1,2,5,15)

mean_map_T$valuefactor <- cut(mean_map_T$change, breaks=breaks_T)
levels(mean_map_T$valuefactor)<-rewrite_label(levels(mean_map_T$valuefactor))
cols_T=rev(brewer.pal(n = 9, name = "RdBu")[c(1:7)])
names(cols_T)<-levels(mean_map_T$valuefactor)

p_T<-ggplot() + theme_dark()+
  geom_point(data = mean_map_T, aes(x = lon, y = lat, color = valuefactor), size=3, shape=15,position="identity")+
  geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0,color='grey40') +
  scale_color_manual(values=cols_T)+
  labs(color=expression(Delta~T~"(K)    "))+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0,0,0,0.7,"cm"))+
  theme(legend.position = "bottom",
        legend.title=element_text(size=17),
        legend.text=element_text(size=17))+
  guides(color = guide_legend(nrow = 4))


#### albedo
summary(mean_map_whole_effect[,"change"])
summary(mean_map_veg_effect[,"change"])
summary(mean_map_snow_effect[,"change"])
breaks_albedo = c(-0.30,-0.15,-0.09,-0.05,-0.02,-0.01,0.01,0.02,0.06)

mean_map_whole_effect$valuefactor <- cut(mean_map_whole_effect$change, breaks=breaks_albedo)
mean_map_veg_effect$valuefactor <- cut(mean_map_veg_effect$change, breaks=breaks_albedo)
mean_map_snow_effect$valuefactor <- cut(mean_map_snow_effect$change,breaks_albedo)

levels(mean_map_whole_effect$valuefactor)<-rewrite_label(levels(mean_map_whole_effect$valuefactor))
levels(mean_map_veg_effect$valuefactor)<-rewrite_label(levels(mean_map_veg_effect$valuefactor))
levels(mean_map_snow_effect$valuefactor)<-rewrite_label(levels(mean_map_snow_effect$valuefactor))

cols_albedo=brewer.pal(n = 11, name = "RdBu")[c(1:8)]
levels_albedo<-levels(mean_map_whole_effect$valuefactor)
names(cols_albedo)<-levels_albedo

## whole_effect
p_whole_effect<-ggplot() + theme_dark()+
  geom_point(data = mean_map_whole_effect, aes(x = lon, y = lat, color = valuefactor), size=3, shape=15,position="identity")+
  geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0,color='grey40') +
  scale_color_manual(values=cols_albedo,limits=levels_albedo)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0,0,0,0.7,"cm"))+
  theme(legend.position = "none")

## veg_effect
p_veg_effect<-ggplot() + theme_dark()+
  geom_point(data = mean_map_veg_effect, aes(x = lon, y = lat, color = valuefactor), size=3, shape=15,position="identity")+
  geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0,color='grey40') +
  scale_color_manual(values=cols_albedo,limits=levels_albedo)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0,0,0,0.7,"cm"))+
  theme(legend.position = "none")

## snow_effect
p_snow_effect<-ggplot() + theme_dark()+
  geom_point(data = mean_map_snow_effect, aes(x = lon, y = lat, color = valuefactor), size=3, shape=15,position="identity")+
  geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0,color='grey40') +
  scale_color_manual(values=cols_albedo,limits=levels_albedo)+
  labs(color=expression(Delta~Albedo~"    "))+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0,0,0,0.7,"cm"))+
  theme(legend.position = "bottom",
        legend.title=element_text(size=17),
        legend.text=element_text(size=17))+
  guides(color = guide_legend(nrow = 4))


## Put them together
p<-ggarrange(p_whole_effect,p_veg_effect,p_snow_effect,p_T,ncol=2,
             labels=c("(a)","(b)","(c)","(d)"),
             label.args = list(gp = grid::gpar(font = 4,cex=1.5)) )

ggsave(file=paste(wd,"Output data/3D var/3D var figure/Fig.3 3D var land change map_mean.png",sep=""),
       p,width=11,height=6.6)


##############################################################################################
###########################    Plot change map for each DO event
##############################################################################################

summary(map_T[,"change"])
summary(map_whole_effect[,"change"])
summary(map_veg_effect[,"change"])
summary(map_snow_effect[,"change"])

breaks_T = c(-20,-5,-2,-1,-0.2, 0.2,1,2,5,20)
breaks_albedo = c(-0.48,-0.21,-0.15,-0.09,-0.05,-0.01,0.01,0.05,0.09,0.15,0.21,0.48)

map_T$valuefactor <- cut(map_T$change, breaks=breaks_T)
map_whole_effect$valuefactor <- cut(map_whole_effect$change, breaks=breaks_albedo)
map_veg_effect$valuefactor <- cut(map_veg_effect$change, breaks=breaks_albedo)
map_snow_effect$valuefactor <- cut(map_snow_effect$change,breaks_albedo)

levels(map_T$valuefactor)<-rewrite_label(levels(map_T$valuefactor))
levels(map_whole_effect$valuefactor)<-rewrite_label(levels(map_whole_effect$valuefactor))
levels(map_veg_effect$valuefactor)<-rewrite_label(levels(map_veg_effect$valuefactor))
levels(map_snow_effect$valuefactor)<-rewrite_label(levels(map_snow_effect$valuefactor))

cols_T=rev(brewer.pal(n = 9, name = "RdBu"))
levels_T<-levels(map_T$valuefactor)
names(cols_T)<-levels_T

cols_albedo=brewer.pal(n = 11, name = "RdBu")
levels_albedo<-levels(map_whole_effect$valuefactor)
names(cols_albedo)<-levels_albedo


p_T<-plot_change_map(data=map_T,cols=cols_T,levels=levels_T,legend_name=expression(Delta~T~"(K)  "))
ggsave(file=paste(wd,"Output data/3D var/3D var figure/3D var land change map_land T.png",sep=""),p_T,width=9,height=12)

p_whole_effect<-plot_change_map(data=map_whole_effect,cols=cols_albedo,levels=levels_albedo,legend_name=expression(Delta~Albedo~"  "))
ggsave(file=paste(wd,"Output data/3D var/3D var figure/3D var land change map_land whole effect.png",sep=""),p_whole_effect,width=9,height=12)

p_veg_effect<-plot_change_map(data=map_veg_effect,cols=cols_albedo,levels=levels_albedo,legend_name=expression(Delta~Albedo~"  "))
ggsave(file=paste(wd,"Output data/3D var/3D var figure/3D var land change map_land veg effect.png",sep=""),p_veg_effect,width=9,height=12)

p_snow_effect<-plot_change_map(data=map_snow_effect,cols=cols_albedo,levels=levels_albedo,legend_name=expression(Delta~Albedo~"  "))
ggsave(file=paste(wd,"Output data/3D var/3D var figure/3D var land change map_land snow effect.png",sep=""),p_snow_effect,width=9,height=12)



##############################################################################################
##############################################################################################
####################                                                 #########################
####################             5. Feedback and gain                #########################
####################                                                 #########################
##############################################################################################
##############################################################################################
#rm(list=ls())

wd<-"D:/PhD Project/Feedback paper_Albedo/Data and codes/"

if(!require(egg)){ install.packages("egg");library(egg)}
if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}

source("D:/PhD Project/Feedback paper_Albedo/Data and codes/Feedback paper_Albedo_Functions.R", encoding = 'UTF-8')


##############################################################################################
############  Get the relationship between global temperature and global land temperature
##############################################################################################

LOVECLIM_T_binned<-read.csv(paste(wd,"Output data/LOVECLIM/LOVECLIM_T_binned.csv",sep=""),row.names=1)
LOVECLIM_T_binned_land<-read.csv(paste(wd,"Output data/LOVECLIM/LOVECLIM_T_binned_land.csv",sep=""),row.names=1)

#Use the self-defined get_global_mean () function to get the cos(lat) weighted mean
#the input should have such columns: c("lon","lat","age","x","err_x")
LOVECLIM_T_binned_mean_globe<-get_global_mean(LOVECLIM_T_binned)
LOVECLIM_T_binned_mean_globe_land<-get_global_mean(LOVECLIM_T_binned_land)

relationship_T<-merge(LOVECLIM_T_binned_mean_globe,LOVECLIM_T_binned_mean_globe_land,by="age")
colnames(relationship_T)<-c("age","T_global","err_T_global","T_global_land","err_T_global_land")

sum<-summary(lm(data=relationship_T,T_global~T_global_land));print(sum)
b<-sum[["coefficients"]]["T_global_land","Estimate"];err_b<-sum[["coefficients"]]["T_global_land","Std. Error"]
round(b,digits=2);round(err_b,digits=2)

p<-ggplot(data=relationship_T,aes(T_global_land,T_global))+theme_bw()+geom_point(size=0.5)+
  labs(x=expression(paste("LOVECLIM global mean land temperature "," ", (degree~C))),
       y=expression(paste("LOVECLIM global mean temperature "," ", (degree~C))))

ggsave(file=paste(wd,"Output data/LOVECLIM/LOVECLIM temperature relationship.png",sep=""),
       p,width=7,height=4)


##############################################################################################
###############################    Get feedback 
##############################################################################################

store<-read.csv(paste(wd,"Output data/3D var/3Dvar_start and end.csv",sep=""),row.names=1)

store$dT_land<-store$max_land_T-store$min_land_T
store$err_dT_land<-sqrt(store$err_max_land_T^2+store$err_min_land_T^2)

store$dT<-b*store$dT_land
store$err_dT<-sqrt(err_b^2*store$dT_land^2 + b^2*store$err_dT_land^2)

store$RF_land_whole_effect<-store$max_land_whole_effect-store$min_land_whole_effect
store$err_RF_land_whole_effect<-sqrt(store$err_max_land_whole_effect^2+store$err_min_land_whole_effect^2)

store$RF_land_veg_effect<-store$max_land_veg_effect-store$min_land_veg_effect
store$err_RF_land_veg_effect<-sqrt(store$err_max_land_veg_effect^2+store$err_min_land_veg_effect^2)

store$RF_land_snow_effect<-store$max_land_snow_effect-store$min_land_snow_effect
store$err_RF_land_snow_effect<-sqrt(store$err_max_land_snow_effect^2+store$err_min_land_snow_effect^2)

write.csv(store,paste(wd,"Output data/3D var/3Dvar_start and end_with dT and RF.csv",sep=""))

#store[which(store$k==12),c("RF_land_veg_effect","err_RF_land_veg_effect")]<-NA
#store[which(store$k==12),]<-NA

###### Deming regression
if(!require(deming)){ install.packages("deming");library(deming)}

fit<-deming(data=store,RF_land_whole_effect~dT-1,ystd=err_RF_land_whole_effect,xstd=err_dT);print(fit)
b_whole_effect<-fit[["coefficients"]][2];ci95_b_whole_effect<-fit[["ci"]][2,]
err_b_whole_effect<-(max(ci95_b_whole_effect)-min(ci95_b_whole_effect))/(2*1.96)

fit<-deming(data=store,RF_land_veg_effect~dT-1,ystd=err_RF_land_veg_effect,xstd=err_dT);print(fit)
b_veg_effect<-fit[["coefficients"]][2];ci95_b_veg_effect<-fit[["ci"]][2,]
err_b_veg_effect<-(max(ci95_b_veg_effect)-min(ci95_b_veg_effect))/(2*1.96)

fit<-deming(data=store,RF_land_snow_effect~dT-1,ystd=err_RF_land_snow_effect,xstd=err_dT);print(fit)
b_snow_effect<-fit[["coefficients"]][2];ci95_b_snow_effect<-fit[["ci"]][2,]
err_b_snow_effect<-(max(ci95_b_snow_effect)-min(ci95_b_snow_effect))/(2*1.96)

b_interact_effect<-b_whole_effect - b_veg_effect - b_snow_effect
err_b_interact_effect<-sqrt(err_b_whole_effect^2 + err_b_veg_effect^2 + err_b_snow_effect^2)

feedback<-cbind.data.frame(name=c("Whole effect","Vegetation effect","Snow effect","Interaction effect"),
                           b=c(b_whole_effect,b_veg_effect,b_snow_effect,b_interact_effect),
                           err_b=c(err_b_whole_effect,err_b_veg_effect,err_b_snow_effect,err_b_interact_effect))

feedback$b_CI95<-paste(round(feedback$b,digits=2),"±",round(feedback$err_b*1.96,digits=2))
write.csv(feedback,paste(wd,"Output data/Feedback/Albedo feedback.csv",sep=""))

##############################################################################################
###############################    Plot feedback
##############################################################################################

summary(store)

p1<-ggplot(data=store,aes(dT,RF_land_whole_effect,label=k))+
  geom_point()+geom_text(vjust=1.5,hjust=-0.2,size=4)+
  geom_abline(slope=b_whole_effect)+
  geom_abline(slope=ci95_b_whole_effect[1],linetype="dashed")+
  geom_abline(slope=ci95_b_whole_effect[2],linetype="dashed")+
  expand_limits(x=0,y=0)+
  theme_bw(base_size = 14)+labs(y= bquote('Whole effect ('~ W ~ m ^-2~')'), 
                  x = "Global mean temperature change (K)")+
  geom_pointrange(aes(xmin=dT-1.96*err_dT, xmax=dT+1.96*err_dT),size=0.3)+
  geom_pointrange(aes(ymin=RF_land_whole_effect-1.96*err_RF_land_whole_effect, 
                      ymax=RF_land_whole_effect+1.96*err_RF_land_whole_effect),size=0.3)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  theme(axis.text=element_text(size=12))


p2<-ggplot(data=store,aes(dT,RF_land_veg_effect,label=k))+
  geom_point()+geom_text(vjust=1.5,hjust=-0.2,size=4)+
  geom_abline(slope=b_veg_effect)+
  geom_abline(slope=ci95_b_veg_effect[1],linetype="dashed")+
  geom_abline(slope=ci95_b_veg_effect[2],linetype="dashed")+
  expand_limits(x=0,y=0)+
  theme_bw(base_size = 14)+labs(y= bquote('Vegetation effect ('~ W ~ m ^-2~')'), 
                  x = "Global mean temperature change (K)")+
  geom_pointrange(aes(xmin=dT-1.96*err_dT, xmax=dT+1.96*err_dT),size=0.3)+
  geom_pointrange(aes(ymin=RF_land_veg_effect-1.96*err_RF_land_veg_effect, 
                      ymax=RF_land_veg_effect+1.96*err_RF_land_veg_effect),size=0.3)+
  theme(axis.text=element_text(size=12))

p3<-ggplot(data=store,aes(dT,RF_land_snow_effect,label=k))+
  geom_point()+geom_text(vjust=1.5,hjust=-0.2,size=4)+
  geom_abline(slope=b_snow_effect)+
  geom_abline(slope=ci95_b_snow_effect[1],linetype="dashed")+
  geom_abline(slope=ci95_b_snow_effect[2],linetype="dashed")+
  expand_limits(x=0,y=0)+
  theme_bw(base_size = 14)+labs(y= bquote('Snow effect ('~ W ~ m ^-2~')'), 
                  x = "Global mean temperature change (K)")+
  geom_pointrange(aes(xmin=dT-1.96*err_dT, xmax=dT+1.96*err_dT),size=0.3)+
  geom_pointrange(aes(ymin=RF_land_snow_effect-1.96*err_RF_land_snow_effect, 
                      ymax=RF_land_snow_effect+1.96*err_RF_land_snow_effect),size=0.3)+
  theme(axis.text=element_text(size=12))

p<-ggarrange(p1,p2,p3,ncol = 2,labels=c("(a)","(b)","(c)"))
ggsave(file=paste(wd,"Output data/Feedback/Fig.4 Albedo feedback regression.jpeg",sep=""),p,width=9,height=8)


##############################################################################################
###############################    Get gain
##############################################################################################

################ IPCC AR6 lambda0
alpha_net<- -1.16 #net feedback parameter
err_alpha_net<- abs(-1.81-(-0.51))/(2*1.645) #90% CI
alpha_albedo<- 0.35 #surface albedo feedback parameter
err_alpha_albedo<- abs(0.10 - 0.60)/(2*1.645) #90% CI

alpha0<- alpha_net-alpha_albedo
err_alpha0<- sqrt(err_alpha_net^2+err_alpha_albedo^2)

lambda0<- -1/alpha0
err_lambda0<- abs(lambda0*err_alpha0/alpha0)

round(lambda0,digits=2)
round(err_lambda0,digits=2)


feedback$g<-feedback$b*lambda0 #the slope (b) is the feedback (c) in equation g=cλ_0
feedback$err_g<-sqrt(feedback$b^2*err_lambda0^2+feedback$err_b^2*lambda0^2)

feedback$g_CI95<-paste(round(feedback$g,digits=2),"±",round(feedback$err_g*1.96,digits=2))
write.csv(feedback,paste(wd,"Output data/Feedback/Albedo feedback.csv",sep=""))

