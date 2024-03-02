rm(list=ls())

wd<-"D:/PhD Project/Feedback paper_Albedo/Data and codes/"

if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(dtw)){ install.packages("dtw");library(dtw)} #package for dynamic time warping

# function to add fake points
fake<-function(age,value,err_value,agebin){
  data<-cbind.data.frame(age,value,err_value)
  data_new<-data.frame(age=seq(round(min(age),digits=0),max(age),by=agebin))
  
  for(row in 1:nrow(data_new)){
    age_new=data_new[row,"age"]
    row_closet_in_data=which.min(abs(age_new-data[,"age"]))
    age_closet_in_data=data[row_closet_in_data,"age"]
    
    if(age_new==age_closet_in_data){
      value_new=value[row_closet_in_data,]
      err_value_new=err_value[row_closet_in_data,]
      
    }else if(age_new>age_closet_in_data){
      y1<-as.numeric(value[row_closet_in_data,])
      y2<-as.numeric(value[row_closet_in_data+1,])
      err_y1<-as.numeric(err_value[row_closet_in_data,])
      err_y2<-as.numeric(err_value[row_closet_in_data+1,])
      x1<-as.numeric(data[row_closet_in_data,"age"])
      x2<-as.numeric(data[row_closet_in_data+1,"age"])
      value_new<-y1+(y2-y1)*(age_new-x1)/(x2-x1) #y1*(1-(age_new-x1)/(x2-x1))+y2*(age_new-x1)/(x2-x1)
      err_value_new<-sqrt( err_y1^2*(1-(age_new-x1)/(x2-x1))^2 + err_y2^2*((age_new-x1)/(x2-x1))^2 )
      
    }else{
      y1<-as.numeric(value[row_closet_in_data-1,])
      y2<-as.numeric(value[row_closet_in_data,])
      err_y1<-as.numeric(err_value[row_closet_in_data-1,])
      err_y2<-as.numeric(err_value[row_closet_in_data,])
      x1<-as.numeric(data[row_closet_in_data-1,"age"])
      x2<-as.numeric(data[row_closet_in_data,"age"])
      value_new<-y1+(y2-y1)*(age_new-x1)/(x2-x1) #y1*(1-(age_new-x1)/(x2-x1))+y2*(age_new-x1)/(x2-x1)
      err_value_new<-sqrt( err_y1^2*(1-(age_new-x1)/(x2-x1))^2 + err_y2^2*((age_new-x1)/(x2-x1))^2 )
      
    }#end of the if function
    
    if(length(value_new)==0){
      data_new[row,colnames(data)]=NA
    }else{
      data_new[row,colnames(data)]=c(age_new,value_new,err_value_new)
    }   
  }#end of the for loop
  
  return(data_new)
}#overall end

########################################################################################
###################  Dynamic time warping
########################################################################################

#Load official start dates which are used to make plots
DO_timing1 <- as.matrix(read.csv(paste(wd,"Input data/Timing/DO_timing1.csv",sep=""), row.names=1, stringsAsFactors = FALSE))

#Load query
recon<-read.csv(paste(wd,"Output data/ACER reconstructions/Reconstruction/recon.csv",sep=""),row.names=1)
recon$lon<-ceiling(recon$lon/5.625)*5.625 #bin in 5.625 degree
recon$lat<-ceiling(recon$lat/5.625)*5.625 #bin in 5.625 degree
recon$age<-ceiling(recon$age/25)*25 #bin in 25 years
mean<-aggregate(data=recon[,c("lon","lat","age","site_id","Tmean","fapar","hveg")],.~lon+lat+age+site_id,FUN=function(x) mean(x)) #get the mean in each bin
err<-aggregate(data=recon[,c("lon","lat","age","site_id","sse_Tmean","sse_fapar","sse_hveg")],.~lon+lat+age+site_id,FUN=function(x) sqrt(sum(x^2))/length(x))#get the err in each bin
recon<-merge(mean,err,by=c("lon","lat","age","site_id"))

#Load reference
model<-read.csv(paste(wd,"Output data/LOVECLIM/LOVECLIM_T_binned_land.csv",sep=""),row.names=1)


#Define the normalization function
norm<-function(x){ (x-mean(x))/sd(x)}

#Adjust
recon_adjusted<-data.frame()

breaks_bin=(DO_timing1[-1]+DO_timing1[-length(DO_timing1)])/2 #use the middle of DO dates to divide the bin
#breaks_bin=c(30000,DO_timing1[5:12],50000)

for(i in unique(recon$site_id)){
  
  each<-recon[which(recon$site_id==i),]
  each<-each[order(each$age),]
  
  if(nrow(each)>2){
    tryCatch(
      {   
        each_adjusted<-data.frame()
        each_fake <- fake(age=each$age,
                          value=each[,c("Tmean","fapar","hveg")],
                          err_value=each[,c("sse_Tmean","sse_fapar","sse_hveg")],
                          agebin=25)

        lon<-unique(each$lon);lat<-unique(each$lat)
        each_model<-model[which(model$lon==lon&model$lat==lat&
                                  model$age>=min(each$age)& model$age<=max(each$age)),]
        
        #Divide
        each_fake$age_zone<-cut(each_fake$age,breaks=breaks_bin)
        each_model$age_zone<-cut(each_model$age,breaks=breaks_bin)
        
        for(t in unique(each_fake$age_zone)){
          sub_each_fake<-each_fake[which(each_fake$age_zone==t),]
          sub_each_model<-each_model[which(each_model$age_zone==t),]
          
          sub_each_fake$Tmean_norm<-norm(sub_each_fake$Tmean)
          sub_each_model$T_norm<-norm(sub_each_model$T)
          
          #Align
          query<-sub_each_fake$Tmean_norm
          reference<-sub_each_model$T_norm
          alignment<-dtw(query,reference,keep=T)
          #plot(alignment,type="twoway")
          
          wq<-warp(alignment,index.reference=FALSE)
          sub_each_adjusted<-sub_each_fake[wq,]
          sub_each_adjusted[,c("lon","lat","age")]<-sub_each_model[,c("lon","lat","age")]
          
          each_adjusted<-rbind.data.frame(each_adjusted,sub_each_adjusted)
          
        }
        
        each_adjusted$site_id<-unique(each$site_id)
        
        
        recon_adjusted<-rbind.data.frame(recon_adjusted,each_adjusted)
        
        p<-ggplot()+theme_bw()+xlim(30000,50000)+
          geom_point(data=each,aes(age,Tmean),size=0.5)+
          geom_line(data=each,aes(age,Tmean),size=0.5)+
          geom_ribbon(data=each,aes(x=age,y=Tmean,ymin=Tmean-sse_Tmean*1.96,ymax=Tmean+sse_Tmean*1.96),alpha=0.2)+
          geom_point(data=each_adjusted,aes(age,Tmean),col="dodgerblue2",size=0.3)+
          geom_line(data=each_adjusted,aes(age,Tmean),col="dodgerblue2",size=0.3)+
          geom_ribbon(data=each_adjusted,aes(x=age,y=Tmean,ymin=Tmean-sse_Tmean*1.96,ymax=Tmean+sse_Tmean*1.96),alpha=0.2,fill="dodgerblue2")+
          geom_point(data=each_model,aes(age,T),col="red",size=0.3)+
          geom_line(data=each_model,aes(age,T),col="red",size=0.3)+
          geom_ribbon(data=each_model,aes(x=age,y=T,ymin=T-sd_T*1.96,ymax=T+sd_T*1.96),alpha=0.2,fill="red")+
          geom_vline(xintercept=DO_timing1[5:12],linetype="dashed")+
          ggtitle(paste("site_id=",i,"; lon=",lon,"; lat=",lat))
        
        ggsave(paste(wd,"Output data/ACER reconstructions/Reconstruction/Reconstruction figure/Recon_Tmean_adjusted/site_id=",i,".png",sep=""),p,width=8,height=6)
        
      },error = function(e){})
  }
  
}

rownames(recon_adjusted)<-1:nrow(recon_adjusted)
write.csv(recon_adjusted,paste(wd,"Output data/ACER reconstructions/Reconstruction/recon_adjusted.csv",sep=""))
