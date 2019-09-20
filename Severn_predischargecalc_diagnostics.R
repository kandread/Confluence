#.rs.restartR()
rm(list=ls())

library(ncdf4)  
#library(raster)
#library(rgdal)
library(tidyverse)
library(tidyr)
library(dplyr)
library(schoolmath)
library(hydroTSM)
library(quantmod)
library(matrixStats)

#PEPSI_folder="G:/My Drive/water_depth_data_Severn/"
#PEPSI_in=nc_open(paste(PEPSI_folder,'Po.nc',sep=""))

#folder="G:/My Drive/water_depth_data_Severn/Netcdfs_folder/"
folder="G:/My Drive/Severn_river/Netcdfs_edited/"

riverreach_files=list.files(folder,pattern='*.nc',full.names=FALSE)

for (files in riverreach_files){
      #files=riverreach_files[1]
      Reach_in=nc_open(paste(folder,files,sep=""))
      QReach <- ncvar_get(Reach_in,"XS_Timseries/Q")
      #Reach_Q<- ncvar_get(Reach_in,"Reach_Timseries/Q")
      S1k <- t(ncvar_get(Reach_in,"Reach_Timeseries/S_1km"))
      S90 <- t(ncvar_get(Reach_in,"Reach_Timeseries/S_90m"))
      H1k <- t(ncvar_get(Reach_in,"XS_Timseries/H_1km"))
      H90 <- t(ncvar_get(Reach_in,"XS_Timseries/H_90m"))
      Wxs <- t(ncvar_get(Reach_in,"XS_Timseries/W"))
      ReachDS <- ncvar_get(Reach_in,"Reach_Timeseries/ID")
      Status <- as.array(seq(1,length(S1k),1)) #90 m
      

      
 ## Filtering the over bank flows using the 2-year flood criteria     
     
      H_Reach=rowMeans(H1k)
      H_Reach[is.negative(H_Reach)] <- 0
      
      # Create dates as a Date class object starting from 2016-01-01
      dates <- seq(as.Date("1990-01-01"), length = 9131, by = "days")
      # Use xts() to create smith
      H_Reach_mean <- xts(x = H_Reach,order.by=dates)
      H_Reach_ann <- to.yearly(H_Reach_mean)
      
      #H_Reach_ann <- H_Reach_ann[,2]
      
      WL_Fdc=fdc(H_Reach_ann)
      f <- splinefun(WL_Fdc[,2],H_Reach_ann[,2])
      OB_WL=f(c(.5))
      H_OBF_s= H_Reach > OB_WL
      Timestamp_OBF= which(H_OBF_s, arr.ind = FALSE, useNames = TRUE)

      H_Keep= H_Reach <= OB_WL
      Timestamp_OBF_Inverse=which(H_Keep, arr.ind = FALSE, useNames = TRUE)
      #H_OBF=H_Reach[H_Reach > OB_WL]
      #WL_Reach=H_Reach[H_Reach < OB_WL]
      H1k=H1k[H_Reach < OB_WL,]
      S1k=S1k[,Timestamp_OBF_Inverse]
      S90=S90[,Timestamp_OBF_Inverse]
      QReach=QReach[,Timestamp_OBF_Inverse]
      Wxs=Wxs[Timestamp_OBF_Inverse,]
      H90=H90[Timestamp_OBF_Inverse,]
      Status[Timestamp_OBF]= -999 #flagged for over-bank flow
      Status[Timestamp_OBF_Inverse]=1
      
      
   
      ## Filtering the high water surface slope  
      
      
      Slope_threshold_up <- 0.003
      Slope_threshold_down <- 0
      
      # Mountainous area check
      SReach_mountainous= S1k > Slope_threshold_up
      Timestamp_Slopeup= which(SReach_mountainous, arr.ind = FALSE, useNames = TRUE)
      
      SReach_Keep= S1k <= Slope_threshold_up
      Timestamp_Slopeup_Inverse=which(SReach_Keep, arr.ind = FALSE, useNames = TRUE)
      
      S1k=S1k[Timestamp_Slopeup_Inverse]
      S90=S90[Timestamp_Slopeup_Inverse]
      QReach=QReach[,Timestamp_Slopeup_Inverse]
      Wxs=Wxs[Timestamp_Slopeup_Inverse,]
      H1k=H1k[Timestamp_Slopeup_Inverse,]
      H90=H90[Timestamp_Slopeup_Inverse,]

      Status[Timestamp_Slopeup]= -99999 #flagged for slope threshold
      Status[Timestamp_Slopeup_Inverse]=1

      
      # Flowing upstream check
      
      
      SReach_flowupstream= S1k < Slope_threshold_down
      Timestamp_flowupstream= which(SReach_flowupstream, arr.ind = FALSE, useNames = TRUE)
      
      SReach_Keep= S1k >= Slope_threshold_down
      Timestamp_flowupstream_Inverse=which(SReach_Keep, arr.ind = FALSE, useNames = TRUE)
      
      
      S1k=S1k[Timestamp_flowupstream_Inverse]
      S90=S90[Timestamp_flowupstream_Inverse]
      QReach=QReach[,Timestamp_flowupstream_Inverse]
      #Qmean=colMeans(QReach)
      Wxs=Wxs[Timestamp_flowupstream_Inverse,]
      H1k=H1k[Timestamp_flowupstream_Inverse,]
      H90=H90[Timestamp_flowupstream_Inverse,]
      
      Status[Timestamp_flowupstream]= -99999 #flagged for slope threshold
      Status[Timestamp_flowupstream_Inverse]=1
      
      # Negative Discharge check
      QReachmean=colMeans(QReach)
      QReach_negative= QReachmean < 0
      Timestamp_Qneg= which(QReach_negative, arr.ind = FALSE, useNames = TRUE)
      
      QReach_Keep= QReachmean >= 0
      Timestamp_Qneg_Inverse=which(QReach_Keep, arr.ind = FALSE, useNames = TRUE)

      S1k=S1k[Timestamp_Qneg_Inverse]
      S90=S90[Timestamp_Qneg_Inverse]
      QReach=QReach[,Timestamp_Qneg_Inverse]
      Qmean=colMeans(QReach)
      Wxs=Wxs[Timestamp_Qneg_Inverse,]
      H1k=H1k[Timestamp_Qneg_Inverse,]
      H90=H90[Timestamp_Qneg_Inverse,]
      
      Status[Timestamp_Qneg]= -99999 #flagged for negative discharge
      Status[Timestamp_Qneg_Inverse]=1
      
 
      
      
      
      # for (i in 1:length(SReach)) {
      #   #i=2
      #   if (SReach[i] > Slope_threshold_up){
      #     
      #     Status[i] = "Slope too high - mountainous area" # flag if slope is larger then 300 cm/km (Frasson et al., 2019)
      #     Flag[i] <- 1
      #     # (will have to come up with an alternative way for orbit duration ~ 10 days)
      #     SReach[i][is.na(SReach[i])] <- 0 
      #     SReach[i] = mean(SReach[i+1],SReach[i-1]) # this works only for non-mountainous areas and daily measuremement 
      #     SReach[i][is.negative(SReach[i])] <- 0} # Data treatment?
      #   else if (SReach[i] < Slope_threshold_down) # flag if slope shows that water is flowing upstream
      #   { #print("Slope is flowing upstream")
      #     Status[i] = "Slope is flowing upstream"
      #     SReach[i] = 0} # Data treatment?
      #   else {
      #     #print("Proceed to discharge algorithm")
      #     Status[i] = "Proceed to discharge algorithm"# NO flag 
      #     SReach[i] = SReach[i]}
      #   
      # }

        
      
      
 
   
      

      nc_close(Reach_in)
      rm(Reach_in)
      rm(H_Reach)



      # rm(SReach)
      # 
      # ncpath <- "G:/My Drive/Severn_river/Netcdfs_MeanReachdischarge/"
      # ncname <- paste("SevernRiver_Reach_",reachname,sep="")  
      # ncfname <- paste(ncpath, ncname, ".nc", sep="")

      filename=substr(files,1,(nchar(files)-3)) 
      reachname=substr(files,19,(nchar(files)-3))




#-- writing the modified netcdf files -----------------
#---------- Creating the Netcdf file ------------------
ncpath <- "G:/My Drive/Severn_river/Prediagnostics_output2/"
ncname <- paste(filename,sep="")  
ncfname <- paste(ncpath, ncname, ".nc", sep="")

#---------- define variables-------------------------------

Xsec90m <- as.array(seq(1,dim(H90)[2],1)) #90 m
Xsec1k <- as.array(seq(1,dim(H90)[2],1)) # 1 km
time <- as.array(seq(1,length(S90),1))
Reach <- as.array(seq(1,1,1))

Xsec90mdim <- ncdim_def("XS_90m","orthogonals",as.integer(Xsec90m))
Xsec1kdim <- ncdim_def("Stageloc_1km","orthogonals",as.integer(Xsec1k))
XSTimestepsdim <- ncdim_def("Time steps","day",as.integer(time))
Reachdim <- ncdim_def("Reach","RiverReaches",as.integer(Reach))


#--------- Adding the priors---------------------------------

#Prior information includes prior distrbution of mean annual flow (Q),
#roughness coefficient (n), bankfull depth (z) and width (w)

fillvalue <- -9999
dlname <- "Width_LISFLOOD_derived"
XS_W <- ncvar_def("XS_Timseries/W","meters",list(Xsec90mdim,XSTimestepsdim),fillvalue,dlname,prec="double")
dlname <- "Water SUrface Elevation_90m"# 
XS_H_90m <- ncvar_def("XS_Timseries/H_90m","meters",list(Xsec90mdim,XSTimestepsdim),fillvalue,dlname,prec="double")
dlname <- "Water SUrface Elevation_1km"
XS_H_1km <- ncvar_def("XS_Timseries/H_1km","meters",list(Xsec1kdim,XSTimestepsdim),fillvalue,dlname,prec="double")
dlname <- "Discharge"
XS_Q <- ncvar_def("XS_Timseries/Q","cubic meters per second",list(Xsec1kdim,XSTimestepsdim),fillvalue,dlname,prec="double")
dlname <- "Mean Discharge"
Reach_Q <- ncvar_def("Reach_Timseries/Q","cubic meters per second",list(Reachdim,XSTimestepsdim),fillvalue,dlname,prec="double")
dlname <- "Slope_1km"
ReachSlope_1km <- ncvar_def("Reach_Timeseries/S_1km","m/m",list(Reachdim,XSTimestepsdim),fillvalue,dlname,prec="double")
dlname <- "Slope_90m"
ReachSlope_90m <- ncvar_def("Reach_Timeseries/S_90m","m/m",list(Reachdim,XSTimestepsdim),fillvalue,dlname,prec="double")
dlname <- "Downstream Reach ID"
Reach_DS <- ncvar_def("Reach_Timeseries/ID","dimensionless",list(Reachdim,Reachdim),fillvalue,dlname,prec="double")


#--------- create netCDF file and put arrays---------------------
ncout <- nc_create(ncfname,list(XS_W,XS_H_90m,XS_H_1km,XS_Q,Reach_Q,ReachSlope_90m,ReachSlope_1km,Reach_DS),force_v4=TRUE)

#---------------- Insert the variables 
ncvar_put(ncout,XS_W,Wxs)
ncvar_put(ncout,XS_H_90m,H90)
ncvar_put(ncout,XS_H_1km,H1k)
ncvar_put(ncout,XS_Q,QReach)
ncvar_put(ncout,Reach_Q,Qmean)
ncvar_put(ncout,ReachSlope_90m,S90)
ncvar_put(ncout,ReachSlope_1km,S1k)
ncvar_put(ncout,Reach_DS,ReachDS)


#---------------- add global attributes
title=paste("Severn SWOT-like data for reach", reachname)
ncatt_put(ncout,0,"title",title)

nc_close(ncout)

rm(QReach)
rm(S1k)
rm(S90)
rm(H1k)
rm(H90)
rm(Wxs)
rm(Status)
}
