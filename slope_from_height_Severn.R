rm(list=ls())

library(raster)
library(rgdal)
library(tidyverse)
library(tidyr)
library(dplyr)




index_1k=readRDS('G:/My Drive/water_depth_data_Severn/processed inputs/1k_spatialindex.rds')
height1k=readRDS('G:/My Drive/water_depth_data_Severn/processed inputs/heights_1k.rds')
index_1k=data.frame(index_1k)%>%
  transmute(northing=V3,easting=V2)

folder="G:/My Drive/water_depth_data_Severn/processed inputs/"
index_90_files=list.files(folder,pattern='*coords.rds',full.names=FALSE)

slopematrix1k=matrix(data=NA,nrow=1,ncol=length(height1k))
slopematrix90=matrix(data=NA,nrow=1,ncol=length(height1k))
slopematrix1k_2=matrix(data=NA,nrow=length(height1k),ncol=length(index_90_files))
#slopematrix90_2=matrix(data=NA,nrow=length(height1k),ncol=length(index_90_files))




for (files in index_90_files){
  
  # files=index_90_files[5]
  heightfilename=paste(substr(files,1,(nchar(files)-13)),"_heights_90m",sep="")
  filename=substr(files,1,(nchar(files)-4)) 
  slopefilename=substr(files,1,(nchar(files)-13))
  index_90=data.frame(readRDS(paste(substr(folder,1,(nchar(folder))), filename,".rds",sep="")) ) 
  #height90=readRDS(paste(substr(folder,1,(nchar(folder))),heightfilename,".rds",sep=""))
  #a=isZero(height90)

  names(index_90)=c('easting','northing')
  
  
  distance_map=mapply(function(x,y) sqrt((x-index_90$easting)^2 +(y-index_90$northing)^2 ), x=index_1k$easting,y=index_1k$northing)
  mapping_index=apply(distance_map,1, function(x) which(x==min(x)))
  mapped_1k=height1k[mapping_index,]
  mapped_1k=as.matrix(mapped_1k)
  distance_map_1k=distance_map[mapping_index]
  
  
  
  j= which(index_90_files == files) 
  
  
  for (t in 1: dim(mapped_1k)[2]){
   
    slope1k=lm( mapped_1k[,t] ~  distance_map_1k)
    #slope90=lm( height90[,t] ~  distance_map_1k)
    slope1k=slope1k[["coefficients"]][2]
    #slope90=slope90[["coefficients"]][2]
    slopematrix1k[1,t]=slope1k 
    #slopematrix90[1,t]=slope90 
    
  }
  
  #slopematrix1k_2[,j]=t(slopematrix1k)
  #slopematrix90_2[,j]=t(slopematrix90)
  slopematrix1k_2=t(slopematrix1k)
  
  
  slopefilename=substr(files,1,(nchar(files)-13))
  saveRDS(slopematrix1k_2,file=paste("G:/My Drive/water_depth_data_Severn/Slopes_folder/",slopefilename,"_slope_1k.rds",sep=""))
  #saveRDS(slopematrix90_2,file=paste("G:/My Drive/water_depth_data_Severn/Slopes_folder/",slopefilename,"_slope_90m.rds",sep=""))
  
}
