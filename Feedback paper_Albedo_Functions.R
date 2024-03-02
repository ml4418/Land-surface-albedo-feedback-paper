########################################################################################
###### Define a function to get average age resolution for each site
########################################################################################

get_average_resolution<-function(data){
  
  colnames(data)<-c("site_id","age")
  resolution<-data.frame(matrix(nrow=0,ncol=2))
  
  for(i in unique(data$site_id)){
    subdata<-data[which(data$site_id==i),]
    subdata<-subdata[order(subdata$age),]
    resolution_subdata<-cbind.data.frame(i,mean(subdata[-1,"age"]-subdata[-nrow(subdata),"age"]))
    resolution<-rbind.data.frame(resolution,resolution_subdata)
  }
  colnames(resolution)<-c("site_id","avg_reso")
  return(resolution)
  
}
#####################################################################################################
###### Define functions to decompose LOVECLIM albedo, which also includes the propagation of error
#####################################################################################################

logit<-function(x){log(x/(1-x))}
err_logit<-function(x,err_x){D<-Deriv(logit,"x");sqrt(D(x)^2 *err_x^2)}

rev_logit<-function(x){1/(1+exp(-x))}
err_rev_logit<-function(x,err_x){D<-Deriv(rev_logit,"x");sqrt(D(x)^2 *err_x^2)}

hold_veg_constant <- function(data) {
  
  colnames(data)<-c("lon","lat","age","albedo","err_albedo","snow_cover","err_snow_cover")
  
  data$logit_albedo<-logit(data$albedo)
  data$err_logit_albedo<-err_logit(data$albedo,data$err_albedo)
  
  data$logit_albedo_snow_removed<-data$logit_albedo-ks*data$snow_cover
  data$err_logit_albedo_snow_removed<-sqrt(data$err_logit_albedo^2 + ks^2*data$err_snow_cover^2+err_ks^2*data$snow_cover^2)
  
  mean_veg_albedo<-aggregate(data=data,logit_albedo_snow_removed~lon+lat,FUN=function(x) mean(x), na.action = na.omit)
  err_mean_veg_albedo<-aggregate(data=data,err_logit_albedo_snow_removed~lon+lat,FUN=function(x) sqrt(sum(x^2))/length(x), na.action = na.omit)
  grid<-merge(x=mean_veg_albedo,y=err_mean_veg_albedo,by=c("lon","lat"))
  
  data$logit_albedo_snow_removed<-NULL
  data$err_logit_albedo_snow_removed<-NULL
  data<-merge(data,grid,by=c("lon","lat"))
  
  data$logit_albedo_snow_added<-data$logit_albedo_snow_removed+ks*data$snow_cover
  data$err_logit_albedo_snow_added<-sqrt(data$err_logit_albedo_snow_removed^2 + ks^2*data$err_snow_cover^2+err_ks^2*data$snow_cover^2)
  
  data$albedo_veg_constant<-rev_logit(data$logit_albedo_snow_added)
  data$err_albedo_veg_constant<-err_rev_logit(data$logit_albedo_snow_added,data$err_logit_albedo_snow_added)
  
  #plot(data=data,albedo_veg_constant~albedo);abline(0,1)
  output<-data[,c("albedo_veg_constant","err_albedo_veg_constant")]
  return(output)
}

hold_snow_constant <- function(data) {
  
  colnames(data)<-c("lon","lat","age","albedo","err_albedo","snow_cover","err_snow_cover")
  
  data$logit_albedo<-logit(data$albedo)
  data$err_logit_albedo<-err_logit(data$albedo,data$err_albedo)
  
  data$logit_albedo_snow_removed<-data$logit_albedo-ks*data$snow_cover
  data$err_logit_albedo_snow_removed<-sqrt(data$err_logit_albedo^2 + ks^2*data$err_snow_cover^2+err_ks^2*data$snow_cover^2)
  
  mean_snow_cover<-aggregate(data=data,snow_cover~lon+lat,FUN=function(x) mean(x), na.action = na.omit)
  err_mean_snow_cover<-aggregate(data=data,err_snow_cover~lon+lat,FUN=function(x) sqrt(sum(x^2))/length(x), na.action = na.omit)
  grid<-merge(x=mean_snow_cover,y=err_mean_snow_cover,by=c("lon","lat"))
  
  data$snow_cover<-NULL
  data$err_snow_cover<-NULL
  data<-merge(data,grid,by=c("lon","lat"))
  
  data$logit_albedo_snow_added<-data$logit_albedo_snow_removed+ks*data$snow_cover
  data$err_logit_albedo_snow_added<-sqrt(data$err_logit_albedo_snow_removed^2 + ks^2*data$err_snow_cover^2+err_ks^2*data$snow_cover^2)
  
  data$albedo_snow_constant<-rev_logit(data$logit_albedo_snow_added)
  data$err_albedo_snow_constant<-err_rev_logit(data$logit_albedo_snow_added,data$err_logit_albedo_snow_added)
  
  #plot(data=data,albedo_snow_constant~albedo);abline(0,1)
  output<-data[,c("albedo_snow_constant","err_albedo_snow_constant")]
  return(output)
}

########################################################################################
###### Define a function to weight seasonal albedo by seasonal ssr
########################################################################################

get_ssr_weighted_albedo<-function(data,weight){
  data<-merge(data,weight,by=c("lon","lat"))
  data$albedo<-(data$albedo_1*data$ssr_1+
                  data$albedo_2*data$ssr_2+
                  data$albedo_3*data$ssr_3+
                  data$albedo_4*data$ssr_4)/(data$ssr_1+data$ssr_2+data$ssr_3+data$ssr_4)
  data$err_albedo<-sqrt(data$err_albedo_1^2*data$ssr_1^2+
                          data$err_albedo_2^2*data$ssr_2^2+
                          data$err_albedo_3^2*data$ssr_3^2+
                          data$err_albedo_4^2*data$ssr_4^2)/(data$ssr_1+data$ssr_2+data$ssr_3+data$ssr_4)
  return(data)
}

#Define functions to compose reconstructed albedo, which also includes the propagation of error
falbedo <- function(fapar,snow_cover,hveg,b,kf,ks,kh) {
  logit=b+kf*(1-fapar)+ks*snow_cover+kh*hveg
  1/(1+exp(-logit))
}
err_falbedo<-function(fapar,snow_cover,hveg,b,kf,ks,kh,
                      err_fapar,err_snow_cover,err_hveg,err_b,err_kf,err_ks,err_kh){
  D_fapar<-Deriv(falbedo,"fapar")
  D_snow_cover<-Deriv(falbedo,"snow_cover")
  D_hveg<-Deriv(falbedo,"hveg")
  D_b<-Deriv(falbedo,"b")
  D_kf<-Deriv(falbedo,"kf")
  D_ks<-Deriv(falbedo,"ks")
  D_kh<-Deriv(falbedo,"kh")
  sqrt( D_fapar(fapar,snow_cover,hveg, b,kf,ks,kh)^2 *err_fapar^2 +
          D_snow_cover(fapar,snow_cover,hveg, b,kf,ks,kh)^2 *err_snow_cover^2 +
          D_hveg(fapar,snow_cover,hveg, b,kf,ks,kh)^2 *err_hveg^2 +
          D_b(fapar,snow_cover,hveg, b,kf,ks,kh)^2 *err_b^2 +
          D_kf(fapar,snow_cover,hveg, b,kf,ks,kh)^2 *err_kf^2 +
          D_ks(fapar,snow_cover,hveg, b,kf,ks,kh)^2 *err_ks^2 +
          D_kh(fapar,snow_cover,hveg, b,kf,ks,kh)^2 *err_kh^2 )
}
compose_albedo<-function(data,para_albedo,weight){
  
  for(i in 1:4){
    snow_cover<-data[,paste("snow_cover_",i,sep="")]
    err_snow_cover<-data[,paste("err_snow_cover_",i,sep="")]
    data[,paste("albedo_",i,sep="")]<-falbedo(fapar=data$fapar,
                                              snow_cover=snow_cover,
                                              hveg=data$hveg, 
                                              b=para_albedo["b","Estimate"],
                                              kf=para_albedo["kf","Estimate"],
                                              ks=para_albedo["ks","Estimate"],
                                              kh=para_albedo["kh","Estimate"])
    data[,paste("err_albedo_",i,sep="")]<-err_falbedo(fapar=ifelse(data$fapar>0,data$fapar,10^(-9)), #there might be NAs when it's <=0
                                                      snow_cover=ifelse(snow_cover>0,snow_cover,10^(-9)), #there might be NAs when it's <=0
                                                      hveg=ifelse(data$hveg>0,data$hveg,10^(-9)), #there might be NAs when it's <=0
                                                      b=para_albedo["b","Estimate"],
                                                      kf=para_albedo["kf","Estimate"],
                                                      ks=para_albedo["ks","Estimate"],
                                                      kh=para_albedo["kh","Estimate"],
                                                      err_fapar=data$sse_fapar,
                                                      err_snow_cover=err_snow_cover,
                                                      err_hveg=data$sse_hveg, 
                                                      err_b=para_albedo["b","Std..Error"],
                                                      err_kf=para_albedo["kf","Std..Error"],
                                                      err_ks=para_albedo["ks","Std..Error"],
                                                      err_kh=para_albedo["kh","Std..Error"])
  }
  
  return(data)
}

########################################################################################
###### Define a function to get the global mean, which also includes the propagation of error
########################################################################################

get_global_mean<-function(Data){
  
  colnames(Data)<-c("lon","lat","age","value","err_value")
  
  mean<-data.frame()
  
  for(age in unique(Data$age)){
    
    data=Data[which(Data$age==age),]
    data$weight<-cos(data$lat*2*pi/360)
    value_w<-sum(data$value*data$weight)/sum(data$weight)
    err_value_w<-sqrt(sum(data$err_value^2*data$weight^2,na.rm=T))/sum(data$weight)
    each_mean<-cbind.data.frame(age,value_w,err_value_w)
    mean<-rbind.data.frame(mean,each_mean)
    
  }
  return(mean)
}

########################################################################################
###### Define a function to identify the minimum and maximum
########################################################################################

get_event<-function(data,t_DO,win,min_left,min_right,duration){
  
  colnames(data)<-c("Age","Value","Err")
  tnow<-t_DO
  
  event<-data[which(data$Age>=(tnow-win)&data$Age<=(tnow+win)),]
  colnames(event)<-c("Age","Value","Err");event<-na.omit(event)
  
  eventmin<-event[which(event$Age>=(tnow-min_left)&event$Age<=(tnow+min_right)),]
  min<-min(eventmin[,"Value"])
  err_min<-eventmin[which.min(eventmin[,"Value"]),"Err"]
  t_start<-eventmin[which.min(eventmin[,"Value"]),"Age"]
  
  eventmax<-event[which(event$Age>=(t_start-duration)&event$Age<=t_start),]
  max<-max(eventmax[,"Value"])
  err_max<-eventmax[which.max(eventmax[,"Value"]),"Err"]
  t_end<-eventmax[which.max(eventmax[,"Value"]),"Age"]
  
  outputdata<-cbind.data.frame(t_start,t_end,min,max,err_min,err_max)
  return(outputdata)
  
}

########################################################################################
###### Define a function to get the change map (end minus start)
########################################################################################

get_change_map<-function(k,data,t_start,t_end){
  
  colnames(data)<-c("lon","lat","age","value","err_value")
  map_start<-data[which(data$age==t_start),]
  map_end<-data[which(data$age==t_end),]
  
  map<-merge(x=map_start,y=map_end,by=c("lon","lat"))
  map$change<-map$value.y - map$value.x
  map$err_change<-sqrt(map$err_value.y^2 + map$err_value.x^2)
  map$k<-k
  
  output<-map[,c("k","lon","lat","change","err_change")]
  
  return(output)
  
}

########################################################################################
###### Define a function to get common legend
########################################################################################

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

########################################################################################
###### Define a function to plot change map
########################################################################################

plot_change_map<-function(data,cols,levels,legend_name){
  
  list<-list()
  
  for(k in 5:12){
    
    map_change<-data[which(data$k==k),]
    
    p_each<-ggplot() + theme_dark()+
      geom_point(data = map_change, aes(x = lon, y = lat, color = valuefactor), size = 3, shape=15,position="identity")+
      geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0,color='grey40') +
      scale_color_manual(values=cols,limits=levels)+
      labs(color=legend_name)+
      ggtitle(paste("DO",k))+
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())+
      theme(legend.position = "bottom",
            legend.title = element_text(size=13),legend.text = element_text(size=13))+
      guides(color = guide_legend(nrow = 2))
    
    legend<-get_legend(p_each)
    list[[k-4]]<-p_each
    gc()
  }
  
  p1<-list[[1]]+theme(legend.position = "none")
  p2<-list[[2]]+theme(legend.position = "none")
  p3<-list[[3]]+theme(legend.position = "none")
  p4<-list[[4]]+theme(legend.position = "none")
  p5<-list[[5]]+theme(legend.position = "none")
  p6<-list[[6]]+theme(legend.position = "none")
  p7<-list[[7]]+theme(legend.position = "none")
  p8<-list[[8]]+theme(legend.position = "none")
  
  DO_map<-grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,legend,ncol=2,
                       layout_matrix=rbind(cbind(1,2),cbind(1,2),cbind(3,4),cbind(3,4),
                                           cbind(5,6),cbind(5,6),cbind(7,8),cbind(7,8),9))
  return(DO_map)
}

########################################################################################
###### Define a function to rewrite labels
########################################################################################
rewrite_label<-function(value_factor){
  # Replace characters to convert to ~ notation
  value_factor <- gsub("\\(", "", value_factor)  # Replace ( with blank
  value_factor <- gsub("\\]", "", value_factor)  # Replace ] with blank
  value_factor <- gsub(",", "~", value_factor)    # Replace , with ~
  return(value_factor)
}
