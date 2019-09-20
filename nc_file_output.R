
rm(list=ls())


library(ncdf4)  
library(dplyr)
library(readxl)



#-------- Reading the inputs --------------
#index_1k=data.frame(readRDS('G:/My Drive/water_depth_data_Severn/processed inputs/1k_spatialindex.rds'))
index_1k=readRDS('G:/My Drive/water_depth_data_Severn/processed inputs/1k_spatialindex.rds')
index_1k=data.frame(index_1k)%>%
  transmute(northing=V3,easting=V2)
height1k=readRDS('G:/My Drive/water_depth_data_Severn/processed inputs/heights_1k.rds')
discharge1k=readRDS('G:/My Drive/water_depth_data_Severn/processed inputs/discharge_1k.rds')


folder="G:/My Drive/water_depth_data_Severn/processed inputs/"
index_90_files=list.files(folder,pattern='*coords.rds',full.names=FALSE)


routing_table <- read_excel('G:/My Drive/Severn_river/Routing_table.xlsx')


for (files in index_90_files){
  #files=index_90_files[9]
  
  
  
  heightfilename=paste(substr(files,1,(nchar(files)-13)),"_heights_90m",sep="")
  #height90=data.frame(readRDS(paste(substr(folder,1,(nchar(folder))),heightfilename,".rds",sep="")))
  height90=readRDS(paste(substr(folder,1,(nchar(folder))),heightfilename,".rds",sep=""))
  
  
  filename=substr(files,1,(nchar(files)-4)) 
  reachname=substr(files,1,(nchar(files)-13))
  a=routing_table[which(routing_table$GridID == reachname),c(1,2,3)]
  b=routing_table[which(routing_table$GridID == reachname),c(3)]
  #index_90=readRDS(paste(substr(folder,1,(nchar(folder))), filename,".rds",sep="")) 
  index_90=data.frame(readRDS(paste(substr(folder,1,(nchar(folder))), filename,".rds",sep="")) ) 
  names(index_90)=c('easting','northing')

  distance_map=mapply(function(x,y) sqrt((x-index_90$easting)^2 +(y-index_90$northing)^2 ), x=index_1k$easting,y=index_1k$northing)
  
  
  mapping_index=apply(distance_map,1, function(x) which(x==min(x)))
  mapped_1k_height=height1k[mapping_index,] #at stage points
  mapped_1k_height=as.matrix(mapped_1k_height)
  mapped_1k_discharge=discharge1k[mapping_index,] # at stage points
  mapped_1k_discharge=as.matrix(mapped_1k_discharge) 
  Mean_discharge=colMeans(mapped_1k_discharge,na.rm = TRUE,dims = 1)
  
  
  
  #Slope_measurements=readRDS('G:/My Drive/water_depth_data_Severn/Slopes_folder/7_3.rds')# one per reach
  slopefolder="G:/My Drive/water_depth_data_Severn/Slopes_folder/"
  Slope_measurements_1k=readRDS(paste(slopefolder,reachname,"_slope_1k.rds",sep=""))# one per reach
  Slope_measurements_90m=readRDS(paste(slopefolder,reachname,"_slope_90m.rds",sep=""))# one per reach
  widthsfolder="G:/My Drive/Severn_river/processed inputs/"                          
  width_measurements=readRDS(paste(widthsfolder,reachname,"_widths_90m.rds",sep="")) #at 90 m resolution
  width_measurements=as.matrix(width_measurements)
  height90=as.matrix(height90)
  
  
  
  #---------- Creating the Netcdf file ------------------
  ncpath <- "G:/My Drive/Severn_river/Netcdfs_edited/"
  ncname <- paste("SevernRiver_Reach_",reachname,sep="")  
  ncfname <- paste(ncpath, ncname, ".nc", sep="")
  
  #---------- define variables-------------------------------
  
  Xsec90m <- as.array(seq(1,dim(index_90)[1],1)) #90 m
  Xsec1k <- as.array(seq(1,dim(mapped_1k_height)[1],1)) # 1 km
  time <- as.array(seq(1,dim(mapped_1k_height)[2],1))
  Reach <- as.array(seq(1,1,1))
  
  Xsec90mdim <- ncdim_def("XS_90m","orthogonals",as.integer(Xsec90m))
  Xsec1kdim <- ncdim_def("Stageloc_1km","orthogonals",as.integer(Xsec1k))
  XSTimestepsdim <- ncdim_def("Time steps","day",as.integer(time))
  Reachdim <- ncdim_def("Reach","RiverReaches",as.integer(Reach))
  
  
  fillvalue <- -9999
  fillvalue2 <- "reachname"
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
  dlname <- "Reach ID"
  Reach_DS <- ncvar_def("Reach_Timeseries/ID","dimensionless",list(Reachdim,Reachdim),fillvalue,dlname,prec="double")
  
  
  #--------- create netCDF file and put arrays---------------------
  ncout <- nc_create(ncfname,list(XS_W,XS_H_90m,XS_H_1km,XS_Q,Reach_Q,ReachSlope_90m,ReachSlope_1km,Reach_DS),force_v4=TRUE)
  
  #---------------- Insert the variables 
  ncvar_put(ncout,XS_W,width_measurements)
  ncvar_put(ncout,XS_H_90m,height90)
  ncvar_put(ncout,XS_H_1km,mapped_1k_height)
  ncvar_put(ncout,XS_Q,mapped_1k_discharge)
  ncvar_put(ncout,Reach_Q,Mean_discharge)
  ncvar_put(ncout,ReachSlope_90m,Slope_measurements_90m)
  ncvar_put(ncout,ReachSlope_1km,Slope_measurements_1k)
  ncvar_put(ncout,Reach_DS,b)
  
  #---------------- add global attributes
  title=paste("Severn SWOT-like data for reach", reachname)
  ncatt_put(ncout,0,"title",title)
  
  nc_close(ncout)
}
