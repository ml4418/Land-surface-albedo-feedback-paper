
if(!require(lbfgs)){ install.packages("lbfgs");library(lbfgs)} #load the package to do Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS)
if(!require(GPBayes)){ install.packages("GPBayes");library(GPBayes)} #load the package to get Bessel correlation
if(!require(pracma)){ install.packages("pracma");library(pracma)} #load the package to get Jacobian

########################################################################################
###################### Define the bessel correlation function to get CLs and CLt
########################################################################################

angle_between_points<-function(loc1,loc2){
  lon1<-loc1$lon;lon2<-loc2$lon
  lat1<-loc1$lat;lat2<-loc2$lat
  return(acos(sin(lat1) * sin(lat2)  + cos(lat1) * cos(lat2) * cos(lon2 - lon1)))
}
scaled_bessel_correlation<-function(theta,L,a){
  k=a/L
  scaled_angle=abs(sin(theta/2))
  bessel_input=k*scaled_angle
  
  if(isTRUE(scaled_angle==0)==TRUE){
    output=1
  }else{
    bessel=GPBayes::BesselK(nu=1,z=bessel_input)
    output=bessel_input*bessel
  }
  return(output)
}

##################################################################################################
################# Define function to get the spatial correlation Ls and temporal correlation Lt
#################################################################################################

get_CLs<-function(loc,Ls){
  loc<-loc*pi/180 #convert lon and lat to radians
  
  #Get the angle distance on a great circle of the Earth between locs
  angle_dist<-data.frame(matrix(nrow=nrow(loc),ncol=nrow(loc)))
  for(i in 1:nrow(angle_dist)){
    if(i%%100==0){print(i)}#show progress
    for(j in 1:i){
      angle_dist[i,j]=angle_between_points(loc[i,],loc[j,])
      angle_dist[j,i]=angle_dist[i,j]
    }
  }
  
  #Use the distance as input and calculate the Mat??rn covariance
  CLs<-data.frame(matrix(nrow=nrow(loc),ncol=nrow(loc)))
  as=6371 #radius of the Earth
  
  for(i in 1:nrow(CLs)){
    if(i%%100==0){print(i)}#show progress
    for(j in 1:i){
      CLs[i,j]=scaled_bessel_correlation(theta=angle_dist[i,j],L=Ls,a=as)
      CLs[j,i]=CLs[i,j]
    }
  }
  return(CLs)
}

get_Lt<-function(time,Lt){
  #Get the time distance between points
  time_dist<-data.frame(matrix(nrow=length(time),ncol=length(time)))
  for(i in 1:length(time)){
    if(i%%100==0){print(i)}#show progress
    for(j in 1:i){
      time_dist[i,j]= abs(time[i]-time[j])
      time_dist[j,i]=time_dist[i,j]
    }
  }
  #Use the time interval as input and calculate the Mat??rn covariance
  CLt<-data.frame(matrix(nrow=length(time),ncol=length(time)))
  at=(max(time)-min(time))/(2*pi)
  
  for(i in 1:nrow(CLt)){
    if(i%%100==0){print(i)}#show progress
    for(j in 1:i){
      CLt[i,j]=scaled_bessel_correlation(theta=time_dist[i,j],L=Lt,a=at)
      CLt[j,i]=CLt[i,j]
    }
  }
  return(CLt)
}


#################################################################################################
########################### Define the function to get estimation
#################################################################################################

get_xa<-function(model,obs,Ls,Lt){
  
  ########### Define function to get the correlation matrix C and covariance matrix B
  #Spatial correlation CLs
  loc<-unique(model[,c("lon","lat")])
  print(all.equal(rep(loc$lon,times=nrow(model)/nrow(loc)),model$lon)) #check whether the location works for all the ages
  print(all.equal(rep(loc$lat,times=nrow(model)/nrow(loc)),model$lat)) #check whether the location works for all the ages
  CLs = get_CLs(loc,Ls=Ls)
  
  
  #Temporal correlation CLt
  #Get the time information
  time<-seq(from=min(model$age),to=max(model$age),by=25)
  print(all.equal(rep(time,each=nrow(model)/length(time)),model$age)) #check whether the age works for all the locations
  CLt = get_Lt(time,Lt=Lt)
  
  
  #The prior error correlation matrix C
  C=kronecker(as.matrix(CLs), as.matrix(CLt), FUN = "*")#the Kronecker product of matrices
  C[is.na(C)] <- 0
  
  #The prior error covariance matrix B
  Sigma<-diag(sqrt(model$err))
  Sigma[is.na(Sigma)] <- 0
  
  B<-Sigma %*% as.matrix(C) %*% Sigma
  B_squared<-sqrt(B)
  
  ################ Write the cost function and get the best solution xa
  y=as.matrix(obs$value)
  R=diag(obs$err)
  R_inverse=solve(R)
  
  xb=as.matrix(model$value)
  
  match_model_obs=merge(x=model,y=obs,by=c("lon","lat","age"))
  match.model.seq=match_model_obs$X.x
  match.obs.seq=match_model_obs$X.y
  
  #define cost function J(w) and the gradient of cost function ∇J(w)
  w<-matrix(nrow=nrow(xb),ncol=1)
  h<-function(h_input){
    return(h_input[match.model.seq])
  }
  cost_w<-function(w){
    x=xb+B_squared%*%w
    diff=y-h(x)
    Jw=(1/2)*(t(w)%*%w)+(1/2)*( t(diff) %*% R_inverse %*% diff )
    return(as.numeric(Jw))
  }
  gradient_cost_w<-function(w){
    x=xb+B_squared%*%w
    diff=y-h(x)
    Hx=jacobian(h, x)
    gradient_Jw=w-B_squared %*% t(Hx) %*% R_inverse %*% diff
    return(as.numeric(gradient_Jw))
  }
  
  ################ Get w and xa
  output <- lbfgs(call_eval=cost_w,  call_grad=gradient_cost_w,  vars=rep(0,times=nrow(model)))
  model$w=output[["par"]]
  model$xb=model$value
  model$err_xb=model$err
  model$xa=as.vector(model$xb+B_squared %*% model$w)
  
  ################ Get the error of xa
  Hxb=jacobian(h, xb) 
  K=B %*% t(Hxb) %*% solve(Hxb %*% B %*% t(Hxb) + R) #analysis gain matrix
  I=diag(nrow(B)) #identity matrix
  A=(I-K %*% Hxb) %*% B
  model$err_xa=sqrt(diag(A))
  model<-model[,c("lon","lat","age","xb","err_xb","xa","err_xa")]
  
  return(model)
}

