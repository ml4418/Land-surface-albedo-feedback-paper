########################################################################################
########################################################################################
#######################                                           ######################
#######################   1. ACER temperature reconstructions     ######################
#######################                                           ######################
########################################################################################
########################################################################################

rm(list=ls())

wd<-"D:/PhD Project/Feedback paper_Albedo/Data and codes/"

if(!require(fxTWAPLS)){ install.packages("fxTWAPLS");library(fxTWAPLS)}
source("D:/PhD Project/Feedback paper_Albedo/Data and codes/Feedback paper_Albedo_Functions.R", encoding = 'UTF-8')


##############################################################################
####################### Pre-process the modern pollen record
##############################################################################

#Import dataset
smpdsv2_pollen_counts_amalgamated <- read.csv(paste(wd,"Input data/Pollen/Modern pollen/smpdsv2_pollen_counts_amalgamated.csv",sep=""),encoding = "UTF-8")
smpdsv2_metadata <- read.csv(paste(wd,"Input data/Pollen/Modern pollen/smpdsv2_metadata.csv",sep=""),encoding = "UTF-8")
colnames(smpdsv2_pollen_counts_amalgamated)[1]<-"ID_SAMPLE"
colnames(smpdsv2_metadata)[1]<-"ID_SITE"

#Merge climate data and pollen data by ID_sample
modern_data<-merge(x=smpdsv2_metadata,y=smpdsv2_pollen_counts_amalgamated,by="ID_SAMPLE",by.y="ID_SAMPLE")
modern_data<-na.omit(modern_data,cols = c("mtco", "mtwa","gdd0","mi")) #remove rows with NAs

taxaColMin <- which(colnames(modern_data) == "Abatia")
taxaColMax <- which(colnames(modern_data) == "Zygophyllaceae")
modern_taxa <- modern_data[, taxaColMin:taxaColMax]
modern_taxa<-modern_taxa/rowSums(modern_taxa)

modern_data_aggregate<-cbind.data.frame(modern_data[,c("longitude","latitude","elevation","mtco", "mtwa","gdd0","mi")],modern_taxa)
modern_data_aggregate<-aggregate(data=modern_data_aggregate,.~longitude+latitude+elevation,FUN=mean)#aggregate duplicate rows
write.csv(modern_data_aggregate,paste(wd,"Input data/Pollen/Modern pollen/modern_pollen_aggregated by lon lat elv.csv",sep=""))


##############################################################################
####################### Training
##############################################################################

modern_data<-read.csv(paste(wd,"Input data/Pollen/Modern pollen/modern_pollen_aggregated by lon lat elv.csv",sep=""),row.names=1)
taxaColMin <- which(colnames(modern_data) == "Abatia")
taxaColMax <- which(colnames(modern_data) == "Zygophyllaceae")
modern_taxa <- modern_data[, taxaColMin:taxaColMax]

#combine Quercus
modern_taxa$Quercus.combined<-modern_taxa$Quercus+modern_taxa$Quercus.deciduous+modern_taxa$Quercus.evergreen
modern_taxa[,c("Quercus","Quercus.deciduous","Quercus.evergreen")]<-NULL

#tidy modern data
modern_taxa <- modern_taxa[, which(colSums(modern_taxa>0)>=10)] #remove taxa less than 10 occurrences
modern_taxa <- modern_taxa[rowSums(modern_taxa)>0,] #remove rows with no taxa information
modern_taxa<-modern_taxa/rowSums(modern_taxa) #re-sum to 1
modern_taxa<-modern_taxa[,sort(colnames(modern_taxa))] #re-order column by names

#training
# MTCO
fit_mtco <- fxTWAPLS::TWAPLS.w2(modern_taxa, modern_data$mtco, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)
#MTWA
fit_mtwa <- fxTWAPLS::TWAPLS.w2(modern_taxa, modern_data$mtwa, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)

fxTWAPLS::plot_train(fit_mtco,3);fxTWAPLS::plot_residuals(fit_mtco,3)
fxTWAPLS::plot_train(fit_mtwa,3);fxTWAPLS::plot_residuals(fit_mtwa,3)

#############################################################################################
############################  Cross validation 
#############################################################################################

#This part takes too long so it's run on HPC
if(!require(foreach)){install.packages("foreach");library(foreach)}
if(!require(doParallel)){install.packages("doParallel");library(doParallel)}
`%>%` <- magrittr::`%>%`

CPUS<-detectCores()-1

point <- modern_data[, c("longitude", "latitude")]
dist <- fxTWAPLS::get_distance(point, cpus = CPUS)%>% fxTWAPLS::pb()
write.csv(dist, paste(wd,"Output data/ACER reconstructions/Cross validation/distance.csv",sep=""))

# Pseudo removed leave out cross validation
dist<-read.csv(paste(wd,"Output data/ACER reconstructions/Cross validation/distance.csv",sep=""),row.names=1)
gc()
options(future.globals.maxSize = 3 * 1024^3)

pseudo_mtco<-fxTWAPLS::get_pseudo(dist, modern_data$mtco, cpus = CPUS) %>% fxTWAPLS::pb()  
pseudo_mtwa<-fxTWAPLS::get_pseudo(dist, modern_data$mtwa, cpus = CPUS) %>% fxTWAPLS::pb()  

if(!require(rlist)){install.packages("rlist");library(rlist)} 
setwd(paste(wd,"Output data/ACER reconstructions/Cross validation",sep=""))
rlist::list.save(pseudo_mtco, 'pseudo_mtco.rdata')
rlist::list.save(pseudo_mtwa, 'pseudo_mtwa.rdata')

# Pseudo removed leave out cross validation
pseudo_mtco <- rlist::list.load('pseudo_mtco.rdata')
pseudo_mtwa <- rlist::list.load('pseudo_mtwa.rdata')

#tf2 pspline
if(!require(foreach)){install.packages("foreach");library(foreach)}

cv_mtco <- fxTWAPLS::cv.pr.w(modern_taxa,
                             modern_data$mtco,
                             nPLS = 10,
                             fxTWAPLS::TWAPLS.w2,
                             fxTWAPLS::TWAPLS.predict.w,
                             pseudo_mtco,
                             usefx = TRUE,
                             fx_method = "pspline",
                             bin = 0.02,
                             cpus = CPUS,
                             test_mode = FALSE)  %>% fxTWAPLS::pb()  

cv_mtwa <- fxTWAPLS::cv.pr.w(modern_taxa,
                             modern_data$mtwa,
                             nPLS = 10,
                             fxTWAPLS::TWAPLS.w2,
                             fxTWAPLS::TWAPLS.predict.w,
                             pseudo_mtwa,
                             usefx = TRUE,
                             fx_method = "pspline",
                             bin = 0.02,
                             cpus = CPUS,
                             test_mode = FALSE)   %>% fxTWAPLS::pb()  


######### Check the cross validation result
#random test
cv_mtco <- read.csv(paste(wd,"Output data/ACER reconstructions/Cross validation/cv_mtco with Quercus combined.csv",sep=""), row.names=1)
cv_mtwa <- read.csv(paste(wd,"Output data/ACER reconstructions/Cross validation/cv_mtwa with Quercus combined.csv",sep=""), row.names=1)

rand_mtco<-as.data.frame(fxTWAPLS::rand.t.test.w(cv_mtco))
rand_mtwa<-as.data.frame(fxTWAPLS::rand.t.test.w(cv_mtwa))

rand_mtco$name<-"mtco"
rand_mtwa$name<-"mtwa"

rand<-rbind.data.frame(rand_mtco[1:5,],rand_mtwa[1:5,])
write.csv(rand, paste(wd,"Output data/ACER reconstructions/Cross validation/random t test result for mtco and mtwa with Quercus combined.csv",sep=""))

##############################################################################
####################### Reconstruction using ACER fossil pollen
##############################################################################
nsig_mtco<-3
nsig_mtwa<-3

################ Load fossil pollen and get it formatted

setwd(paste(wd,"Input data/Pollen/Fossil pollen/",sep=""))
listcsv <- dir(pattern = "*.csv") # creates the list of all the csv files in the directory

allcore<-data.frame()
for(k in 1:length(listcsv)){
  
  #load
  fossil_pollen<-read.csv(listcsv[k])
  
  #tidy
  each_core0<-fossil_pollen[,11:ncol(fossil_pollen)] #get only taxa information
  each_core0[is.na(each_core0)]<-0 #replace NA with 0
  each_core<-each_core0
  
  #combine Quercus
  col_Quercus<-grep("Quercus",colnames(each_core))
  if(length(col_Quercus)==1){
    each_core$Quercus.combined<-each_core[,col_Quercus]
  }else{
    each_core$Quercus.combined<-rowSums(each_core[,col_Quercus])
  }
  each_core[,col_Quercus]<-NULL
  
  #get it formatted as modern pollen
  taxa_to_delete<-colnames(each_core)[!(colnames(each_core)%in% colnames(modern_taxa))]
  each_core[,taxa_to_delete]<-NULL
  taxa_to_add<-colnames(modern_taxa)[!(colnames(modern_taxa)%in% colnames(each_core))]
  each_core[,taxa_to_add]<-0
  each_core<-each_core[,order(colnames(each_core))]
  each_core<-each_core/rowSums(each_core)
  
  each_core<-cbind.data.frame(fossil_pollen[,1:10],each_core)
  
  allcore<-rbind.data.frame(allcore,each_core)
}
core<-allcore[,11:ncol(allcore)]

############################# Reconstruction

fossil_mtco<-fxTWAPLS::TWAPLS.predict.w(fit_mtco,core)
fossil_mtwa<-fxTWAPLS::TWAPLS.predict.w(fit_mtwa,core)

#Get the sample specific errors, use nboot=1000
`%>%` <- magrittr::`%>%`
sse_mtco<-fxTWAPLS::sse.sample(modern_taxa=modern_taxa,
                               modern_climate=modern_data$mtco,
                               fossil_taxa=core,
                               trainfun=fxTWAPLS::TWAPLS.w2,
                               predictfun=fxTWAPLS::TWAPLS.predict.w,
                               nboot=1000,
                               nPLS=5,
                               nsig=nsig_mtco,
                               usefx=TRUE,
                               fx_method = "pspline",
                               bin=0.02,
                               cpus = 6) %>% fxTWAPLS::pb()
sse_mtwa<-fxTWAPLS::sse.sample(modern_taxa=modern_taxa,
                               modern_climate=modern_data$mtwa,
                               fossil_taxa=core,
                               trainfun=fxTWAPLS::TWAPLS.w2,
                               predictfun=fxTWAPLS::TWAPLS.predict.w,
                               nboot=1000,
                               nPLS=5,
                               nsig=nsig_mtwa,
                               usefx=TRUE,
                               fx_method = "pspline",
                               bin=0.02,
                               cpus = 6) %>% fxTWAPLS::pb()

#Use the last significant number of components
core_sig<-cbind.data.frame(allcore[,1:10],
                           fossil_mtco[["fit"]][,nsig_mtco],sse_mtco,
                           fossil_mtwa[["fit"]][,nsig_mtwa],sse_mtwa)
colnames(core_sig)[11:ncol(core_sig)]<-c("mtco","sse_mtco","mtwa","sse_mtwa")

#Remove rows with wrong ages -9999
core_sig[core_sig==-9999]=NA 
core_sig<-na.omit(core_sig)

write.csv(core_sig,paste(wd,"Output data/ACER reconstructions/Reconstruction/core_sig_mtco_mtwa.csv",sep=""))


########################################################################################
########################################################################################
#######################                                           ######################
#######################   2. ACER vegetation reconstructions      ######################
#######################                                           ######################
########################################################################################
########################################################################################

rm(list=ls())

wd<-"D:/PhD Project/Feedback paper_Albedo/Data and codes/"

if(!require(fxTWAPLS)){ install.packages("fxTWAPLS");library(fxTWAPLS)}
source("D:/PhD Project/Feedback paper_Albedo/Data and codes/Feedback paper_Albedo_Functions.R", encoding = 'UTF-8')


##############################################################################
####################### Pre-process the modern pollen record
##############################################################################

####################### Load fapar and hveg
#Import annual maximum fapar data and get the mean over 1982~2016
load(paste(wd,"Input data/Decomposition of albedo/fapar95percentile 1982~2016 0.5Deg.rda",sep=""))
fapar_mean<-apply(fapar95percentile, c(1,2), FUN=function(x) mean(x))
fapar_mean<-as.data.frame(fapar_mean)
lon<-seq(-179.75,179.75,by=0.5) 
lat<- -seq(-89.75,89.75,by=0.5)
colnames(fapar_mean)<-lon
fapar_mean$lat<-lat

fapar<-melt(fapar_mean,id.vars = c("lat"))
colnames(fapar)<-c("lat","lon","fapar")
fapar$lon<-as.numeric(as.character(fapar$lon))

#Import vegetation height data at 2005
hveg_raster<-brick(paste(wd,"Input data/Decomposition of albedo/Tree height 2005 0.5Deg.nc",sep=""))
hveg_raster <- rasterToPoints(hveg_raster, spatial=TRUE)# Convert raster to SpatialPointsDataFrame
hveg_raster <- spTransform(hveg_raster, CRS(hveg_raster@proj4string@projargs))# reproject sp object
hveg <- data.frame(hveg_raster@data, lon=coordinates(hveg_raster)[,1],lat=coordinates(hveg_raster)[,2])    # Assign coordinates to @data slot                     
colnames(hveg)<-c("hveg","lon","lat")

################### Aggregate smpdsv2 by lon+lat, since fapar and hveg don't differentiate elv
#Import dataset
smpdsv2_pollen_counts_amalgamated <- read.csv(paste(wd,"Input data/Pollen/Modern pollen/smpdsv2_pollen_counts_amalgamated.csv",sep=""),encoding = "UTF-8")
smpdsv2_metadata <- read.csv(paste(wd,"Input data/Pollen/Modern pollen/smpdsv2_metadata.csv",sep=""),encoding = "UTF-8")
colnames(smpdsv2_pollen_counts_amalgamated)[1]<-"ID_SAMPLE"
colnames(smpdsv2_metadata)[1]<-"ID_SITE"

#Merge climate data and pollen data by ID_sample
modern_data<-merge(x=smpdsv2_metadata,y=smpdsv2_pollen_counts_amalgamated,by="ID_SAMPLE",by.y="ID_SAMPLE")

taxaColMin <- which(colnames(modern_data) == "Abatia")
taxaColMax <- which(colnames(modern_data) == "Zygophyllaceae")
modern_taxa <- modern_data[, taxaColMin:taxaColMax]
modern_taxa<-modern_taxa/rowSums(modern_taxa)

modern_data_aggregate<-cbind.data.frame(modern_data[,c("longitude","latitude","elevation","mtco", "mtwa","gdd0","mi")],modern_taxa)
modern_data_aggregate<-aggregate(data=modern_data_aggregate,.~longitude+latitude,FUN=mean)#aggregate duplicate rows

########################## Merge data
for(i in 1:nrow(modern_data_aggregate)){
  if(i%%100==0){print(i)}
  tryCatch({
    sitelon<-as.numeric(modern_data_aggregate[i,"longitude"])
    sitelat<-as.numeric(modern_data_aggregate[i,"latitude"])
    modern_data_aggregate[i,"hveg"]<-mean(hveg[which(hveg$lon>=(sitelon-0.5/2)&hveg$lon<=(sitelon+0.5/2)&
                                                       hveg$lat>=(sitelat-0.5/2)&hveg$lat<=(sitelat+0.5/2)),"hveg"],na.rm=TRUE)
    modern_data_aggregate[i,"fapar"]<-mean(fapar[which(fapar$lon>=(sitelon-0.5/2)&fapar$lon<=(sitelon+0.5/2)&
                                                         fapar$lat>=(sitelat-0.5/2)&fapar$lat<=(sitelat+0.5/2)),"fapar"],na.rm=TRUE)
    
  }, error=function(e){})  
}

modern_data_aggregate<-na.omit(modern_data_aggregate,cols = c("fapar", "hveg")) #remove rows with NAs

write.csv(modern_data_aggregate,paste(wd,"Input data/Pollen/Modern pollen/modern_pollen_aggregated by lon lat_with fapar and hveg.csv",sep=""))

##############################################################################
####################### Training
##############################################################################

modern_data<-read.csv(paste(wd,"Input data/Pollen/Modern pollen/modern_pollen_aggregated by lon lat_with fapar and hveg.csv",sep=""),row.names=1)
taxaColMin <- which(colnames(modern_data) == "Abatia")
taxaColMax <- which(colnames(modern_data) == "Zygophyllaceae")
modern_taxa <- modern_data[, taxaColMin:taxaColMax]

#combine Quercus
modern_taxa$Quercus.combined<-modern_taxa$Quercus+modern_taxa$Quercus.deciduous+modern_taxa$Quercus.evergreen
modern_taxa[,c("Quercus","Quercus.deciduous","Quercus.evergreen")]<-NULL

#tidy modern data
modern_taxa <- modern_taxa[, which(colSums(modern_taxa>0)>=10)] #remove taxa less than 10 occurences
modern_taxa <- modern_taxa[rowSums(modern_taxa)>0,] #remove rows with no taxa information
modern_taxa<-modern_taxa/rowSums(modern_taxa) #re-sum to 1
modern_taxa<-modern_taxa[,sort(colnames(modern_taxa))] #re-order column by names

#training
# fapar
fit_fapar <- fxTWAPLS::TWAPLS.w2(modern_taxa, modern_data$fapar, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.01)
#hveg
fit_hveg <- fxTWAPLS::TWAPLS.w2(modern_taxa, modern_data$hveg, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=1)

fxTWAPLS::plot_train(fit_fapar,2);fxTWAPLS::plot_residuals(fit_fapar,2)
fxTWAPLS::plot_train(fit_hveg,3);fxTWAPLS::plot_residuals(fit_hveg,3)

#############################################################################################
############################  Cross validation 
#############################################################################################

#This part takes too long so it's run on HPC
if(!require(foreach)){install.packages("foreach");library(foreach)}
if(!require(doParallel)){install.packages("doParallel");library(doParallel)}
`%>%` <- magrittr::`%>%`

CPUS<-6

point <- modern_data[, c("longitude", "latitude")]
dist <- fxTWAPLS::get_distance(point, cpus = CPUS)%>% fxTWAPLS::pb()
write.csv(dist, paste(wd,"Output data/ACER reconstructions/Cross validation/distance_fapar_hveg.csv",sep=""))

# Pseudo removed leave out cross validation
dist<-read.csv(paste(wd,"Output data/ACER reconstructions/Cross validation/distance_fapar_hveg.csv",sep=""),row.names=1)
gc()
options(future.globals.maxSize = 3 * 1024^3)

pseudo_fapar<-fxTWAPLS::get_pseudo(dist, modern_data$fapar, cpus = CPUS) %>% fxTWAPLS::pb()  
pseudo_hveg<-fxTWAPLS::get_pseudo(dist, modern_data$hveg, cpus = CPUS) %>% fxTWAPLS::pb()  

if(!require(rlist)){install.packages("rlist");library(rlist)} 
setwd(paste(wd,"Output data/ACER reconstructions/Cross validation",sep=""))
rlist::list.save(pseudo_fapar, 'pseudo_fapar.rdata')
rlist::list.save(pseudo_hveg, 'pseudo_hveg.rdata')

# Pseudo removed leave out cross validation
pseudo_fapar <- rlist::list.load('pseudo_fapar.rdata')
pseudo_hveg <- rlist::list.load('pseudo_hveg.rdata')

#tf2 pspline
if(!require(foreach)){install.packages("foreach");library(foreach)}

cv_fapar <- fxTWAPLS::cv.pr.w(modern_taxa,
                              modern_data$fapar,
                              nPLS = 10,
                              fxTWAPLS::TWAPLS.w2,
                              fxTWAPLS::TWAPLS.predict.w,
                              pseudo_fapar,
                              usefx = TRUE,
                              fx_method = "pspline",
                              bin = 0.01,
                              cpus = CPUS,
                              test_mode = FALSE)  %>% fxTWAPLS::pb()  

cv_hveg <- fxTWAPLS::cv.pr.w(modern_taxa,
                             modern_data$hveg,
                             nPLS = 10,
                             fxTWAPLS::TWAPLS.w2,
                             fxTWAPLS::TWAPLS.predict.w,
                             pseudo_hveg,
                             usefx = TRUE,
                             fx_method = "pspline",
                             bin = 1,
                             cpus = CPUS,
                             test_mode = FALSE)   %>% fxTWAPLS::pb()  


######### Check the cross validation result
#random test
cv_fapar <- read.csv(paste(wd,"Output data/ACER reconstructions/Cross validation/cv_fapar with Quercus combined.csv",sep=""), row.names=1)
cv_hveg <- read.csv(paste(wd,"Output data/ACER reconstructions/Cross validation/cv_hveg with Quercus combined.csv",sep=""), row.names=1)

rand_fapar<-as.data.frame(fxTWAPLS::rand.t.test.w(cv_fapar))
rand_hveg<-as.data.frame(fxTWAPLS::rand.t.test.w(cv_hveg))

rand_fapar$name<-"fapar"
rand_hveg$name<-"hveg"

rand<-rbind.data.frame(rand_fapar[1:5,],rand_hveg[1:5,])
write.csv(rand, paste(wd,"Output data/ACER reconstructions/Cross validation/random t test result for fapar and hveg with Quercus combined.csv",sep=""))

############################################################################################
####################### Reconstruction of fapar and hveg using ACER fossil pollen
############################################################################################
nsig_fapar<-2
nsig_hveg<-3

################ Load fossil pollen and get it formatted

setwd(paste(wd,"Input data/Pollen/Fossil pollen/",sep=""))
listcsv <- dir(pattern = "*.csv") # creates the list of all the csv files in the directory

allcore<-data.frame()
for(k in 1:length(listcsv)){
  
  #load
  fossil_pollen<-read.csv(listcsv[k])
  
  #tidy
  each_core0<-fossil_pollen[,11:ncol(fossil_pollen)] #get only taxa information
  each_core0[is.na(each_core0)]<-0 #replace NA with 0
  each_core<-each_core0
  
  #combine Quercus
  col_Quercus<-grep("Quercus",colnames(each_core))
  if(length(col_Quercus)==1){
    each_core$Quercus.combined<-each_core[,col_Quercus]
  }else{
    each_core$Quercus.combined<-rowSums(each_core[,col_Quercus])
  }
  each_core[,col_Quercus]<-NULL
  
  #get it formatted as modern pollen
  taxa_to_delete<-colnames(each_core)[!(colnames(each_core)%in% colnames(modern_taxa))]
  each_core[,taxa_to_delete]<-NULL
  taxa_to_add<-colnames(modern_taxa)[!(colnames(modern_taxa)%in% colnames(each_core))]
  each_core[,taxa_to_add]<-0
  each_core<-each_core[,order(colnames(each_core))]
  each_core<-each_core/rowSums(each_core)
  
  each_core<-cbind.data.frame(fossil_pollen[,1:10],each_core)
  
  allcore<-rbind.data.frame(allcore,each_core)
}
core<-allcore[,11:ncol(allcore)]

############################# Reconstruction

fossil_fapar<-fxTWAPLS::TWAPLS.predict.w(fit_fapar,core)
fossil_hveg<-fxTWAPLS::TWAPLS.predict.w(fit_hveg,core)

#Get the sample specific errors, use nboot=1000
`%>%` <- magrittr::`%>%`
sse_fapar<-fxTWAPLS::sse.sample(modern_taxa=modern_taxa,
                                modern_climate=modern_data$fapar,
                                fossil_taxa=core,
                                trainfun=fxTWAPLS::TWAPLS.w2,
                                predictfun=fxTWAPLS::TWAPLS.predict.w,
                                nboot=1000,
                                nPLS=5,
                                nsig=nsig_fapar,
                                usefx=TRUE,
                                fx_method = "pspline",
                                bin=0.01,
                                cpus = 6) %>% fxTWAPLS::pb()
sse_hveg<-fxTWAPLS::sse.sample(modern_taxa=modern_taxa,
                               modern_climate=modern_data$hveg,
                               fossil_taxa=core,
                               trainfun=fxTWAPLS::TWAPLS.w2,
                               predictfun=fxTWAPLS::TWAPLS.predict.w,
                               nboot=1000,
                               nPLS=5,
                               nsig=nsig_hveg,
                               usefx=TRUE,
                               fx_method = "pspline",
                               bin=1,
                               cpus = 6) %>% fxTWAPLS::pb()

#Use the last significant number of components
core_sig<-cbind.data.frame(allcore[,1:10],
                           fossil_fapar[["fit"]][,nsig_fapar],sse_fapar,
                           fossil_hveg[["fit"]][,nsig_hveg],sse_hveg)
colnames(core_sig)[11:ncol(core_sig)]<-c("fapar","sse_fapar","hveg","sse_hveg")

#Remove rows with wrong ages -9999
core_sig[core_sig==-9999]=NA 
core_sig<-na.omit(core_sig)

write.csv(core_sig,paste(wd,"Output data/ACER reconstructions/Reconstruction/core_sig_fapar_hveg.csv",sep=""))


########################################################################################
########################################################################################
#######################                                           ######################
#######################          3. Combine and plot              ######################
#######################                                           ######################
########################################################################################
########################################################################################
rm(list=ls())

wd<-"D:/PhD Project/Feedback paper_Albedo/Data and codes/"

if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(ggmap)){ install.packages("ggmap");library(ggmap)}
if(!require(ggsn)){ install.packages("ggsn");library(ggsn)}
if(!require(maps)){ install.packages("maps");library(maps)}
if(!require(mapdata)){ install.packages("mapdata");library(mapdata)}
world <- map_data("world") 

source("D:/PhD Project/Feedback paper_Albedo/Data and codes/Feedback paper_Albedo_Functions.R", encoding = 'UTF-8')

######################### Combine

core_sig_mtco_mtwa<-read.csv(paste(wd,"Output data/ACER reconstructions/Reconstruction/core_sig_mtco_mtwa.csv",sep=""),row.names = 1)
core_sig_fapar_hveg<-read.csv(paste(wd,"Output data/ACER reconstructions/Reconstruction/core_sig_fapar_hveg.csv",sep=""),row.names = 1)
identical(core_sig_mtco_mtwa[,1:10],core_sig_fapar_hveg[,1:10]) #check the sample information
recon<-cbind.data.frame(core_sig_mtco_mtwa,core_sig_fapar_hveg[,c("fapar","sse_fapar","hveg","sse_hveg")])

recon$Tmean<-(recon$mtwa+recon$mtco)/2
recon$sse_Tmean<-sqrt(recon$sse_mtwa^2+recon$sse_mtco^2)/2

################################ Add in site information

siteinfo<-read.csv(paste(wd,"Input data/Pollen/fossil site information.csv",sep=""),row.names = 1)
unique(recon$site_name)

#correct the spelling of names
recon[grep("Abric Roman",recon$site_name),"site_name"]<-"Abric Roman"
recon[which(recon$site_name=="Ca\xe7o"),"site_name"]<-"Cao"
recon[which(recon$site_name=="Col\xf4nia"),"site_name"]<-"Colnia"
recon[which(recon$site_name=="F\xfcramoos"),"site_name"]<-"Framoos"
recon[which(recon$site_name=="Navarr\xe9s"),"site_name"]<-"Navarrs"

#Check whether the site names are all in the file
unique(recon$site_name)%in%unique(siteinfo$site_name)

#Import the site information into the reconstruction file
for(i in 1:nrow(recon)){
  tryCatch({
    name<-recon[i,"site_name"]
    recon[i,c("lon","lat","elv")]<-siteinfo[which(siteinfo$site_name==name),c("long","lat","elev")]
  }, error=function(e){})
}

#Check the age
plot((recon$CLAM_max95+recon$CLAM_min95)/2~recon$CLAM_best);abline(a=0,b=1)
recon$age<-recon$CLAM_best

#Extract only 50~30ka
recon<-recon[which(recon$age>=30000&recon$age<=50000),]
write.csv(recon,paste(wd,"Output data/ACER reconstructions/Reconstruction/recon.csv",sep=""))

############################ plot ACER sites between 50 ka and 30 ka

xlab<-c(expression("90"~degree~W),expression("0"~degree~E),expression("90"~degree~E))
ylab<-c(expression("60"~degree~S),expression("0"~degree~N),expression("60"~degree~N))

p<-ggplot()+theme_bw()+
  geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0,color='black') +
  scale_x_continuous(breaks = c(-90,0,90),labels=xlab)+
  scale_y_continuous(breaks = c(-60,0,60),labels=ylab)+
  geom_point(data=recon,aes(lon,lat),color="red")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())

ggsave(paste(wd,"Output data/ACER reconstructions/ACER sites.png",sep=""),p,width=6,height=3.5)

############################### Check some information

#get the number of samples
nrow(recon)

#get the number of sites
length(unique(recon$site_id))

#get age uncertainty
summary((recon$CLAM_max95-recon$CLAM_min95)/(2*1.96))

#Get age resolution for each site
reso<-get_average_resolution(recon[,c("site_id","age")])
mean(reso$avg_reso)



