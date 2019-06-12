library(doParallel)
library(ncdf4)
library(bamr)
library(dplyr)
library(doParallel)
library(foreach)

testfunction <- function(filename) {
  library(doParallel)
  library(ncdf4)
  library(bamr)
  library(dplyr)
  library(doParallel)
  library(foreach)

  # helper function. 1
  calcdA_mat <- function(w, h) {
    stopifnot(all(dim(w) == dim(h)))
    dA <- w
    for (i in 1:nrow(dA)) {
      dA[i, ] <- calcdA_vec(w[i, ], h[i, ])
    }
    dA
  }

  # helper function 2
  calcdA_vec <- function(w, h) {
    words <- order(w)
    warr <- w[words]
    harr <- h[words]
    delh <- c(0, diff(harr))
    delA <- cumsum(warr * delh)
    dA <- 1:length(w)
    dA[words] <- delA
    dA
  }

  data_in <- nc_open(filename)
  W_obs <- ncvar_get(data_in, "XS_Timseries/W")
  Q_obs <- ncvar_get(data_in, "XS_Timseries/Q")
  H_obs <- ncvar_get(data_in, "XS_Timseries/H_1km")
  S_obs <- ncvar_get(data_in, "Reach_Timeseries/S_1km")

  # W_obs <- W_obs[1:2, , drop=FALSE]
  # H_obs <- H_obs[1:2, , drop=FALSE]

  W_obs= W_obs[, colSums(is.na(W_obs) ) < nrow(W_obs)]
  W_obs[W_obs==0]=1e-6
  # convert S_obs to matrix
  S_obs <- matrix(unlist(S_obs),ncol = length(S_obs), byrow = TRUE)
  S_obs[S_obs==0]=1e-6
  S_obs <- S_obs[rep(1:nrow(S_obs), times=nrow(W_obs)), ]

  dA=calcdA_mat(W_obs,H_obs)
  Q_prior <- abs(mean(Q_obs,na.rm=TRUE))

  bamdata=bamr::bam_data(w=W_obs,dA=dA,s=S_obs,Qhat=as.vector(Q_prior))
  priors=bam_priors(bamdata,logQ_sd = cv2sigma(1))
  #custom priors post Hagemann et al 2017
  priors$lowerbound_logQc=0.01
  priors$upperbound_logn=log(0.05)
  #lower A0 for small rivers
  if (max(dA)<2000){
    priors$lowerbound_A0= 1
    priors$upperbound_A0= 1e3
  }
  run_bam = bam_estimate(bamdata=bamdata,bampriors=priors, variant = 'manning_amhg', iter=10)
  required_data <- as.data.frame(run_bam)
  A_01 <- mean(required_data$`A0[1,1]`)
  #A_02 <- mean(required_data$`A0[1,2]`)

  # get the inverse of log n
  n <- (10^mean(required_data$`logn[1]`))
  current_result = list("A01"=c(A_01), "n"=c(n), "dA"=c(dA), "S"=c(S_obs), "W"=c(W_obs))
  return(current_result)
}

format_estimated_results <- function(results) {
  A0 <- unname(as.matrix(as.data.frame(results[1, ])))
  n_values <- unname(as.matrix(as.data.frame(results[2, ])))
  dA <- unname(as.matrix(as.data.frame(results[3, ])))
  Slope <- unname(as.matrix(as.data.frame(results[4, ])))
  Width <- unname(as.matrix(as.data.frame(results[5, ])))

  dim1_A <- ncdim_def( name = "Rows0", units = "NA", vals = as.double(1:nrow(A0)) )
  dim2_A <- ncdim_def( name = "Cols0", units = "NA", vals = as.double(1:ncol(A0)) )
  #========
  dim1_dA <- ncdim_def( name = "Rows1", units = "NA", vals = as.double(1:nrow(dA)) )
  dim2_dA <- ncdim_def( name = "Cols1", units = "NA", vals = as.double(1:ncol(dA)) )
  #=======
  dim1_S <- ncdim_def( name = "Rows2", units = "NA", vals = as.double(1:nrow(Slope)) )
  dim2_S <- ncdim_def( name = "Cols2", units = "NA", vals = as.double(1:ncol(Slope)) )
  # =====
  dim1_W <- ncdim_def( name = "Rows3", units = "NA", vals = as.double(1:nrow(Width)) )
  dim2_W <- ncdim_def( name = "Cols3", units = "NA", vals = as.double(1:ncol(Width)) )
  #====
  dim1_n <- ncdim_def( name = "Rows4", units = "NA", vals = as.double(1:nrow(n_values)) )
  dim2_n <- ncdim_def( name = "Cols4", units = "NA", vals = as.double(1:ncol(n_values)) )

  var_A <- ncvar_def(name = "A0", units =  "NA", dim = list(dim1_A, dim2_A), missval = -1, longname ="Area")
  var_dA <- ncvar_def(name = "dA", units =  "NA", dim = list(dim1_dA, dim2_dA), missval = -1, longname ="Cross Sectional Area")
  var_S <- ncvar_def(name = "Slope", units =  "NA", dim = list(dim1_S, dim2_S), missval = -1, longname ="Slope")
  var_W <- ncvar_def(name = "Width", units =  "NA", dim = list(dim1_W, dim2_W), missval = -1, longname ="Width")
  var_n <- ncvar_def(name = "n", units =  "NA", dim = list(dim1_n, dim2_n), missval = -1, longname ="n")

  allvars <- list(var_A, var_dA, var_S, var_W, var_n)
  nc_file <- nc_create(output_filename, vars = allvars)

  ncvar_put(nc = nc_file, varid = var_A, vals = A0)
  ncvar_put(nc = nc_file, varid = var_dA, vals = dA)
  ncvar_put(nc = nc_file, varid = var_S, vals = Slope)
  ncvar_put(nc = nc_file, varid = var_W, vals = Width)
  ncvar_put(nc = nc_file, varid = var_n, vals = n_values)
  nc_close(nc_file)
}

# cl <- makePSOCKcluster(detectCores())
cl <- makeCluster(future::availableCores())
#registerDoParallel(cl)
registerDoParallel()
doFuture::registerDoFuture()

filepath <- paste(getwd(), "/data", sep="")
output_filename <- paste(getwd(), "/bam_output.nc", sep="")
all_files <- list.files(filepath, patter="*nc", full.names=TRUE)
# all_files <- all_files[1]

start_time= Sys.time()
results <- foreach(filename=all_files, .combine = cbind,
                   .packages = c("ncdf4", "doParallel","foreach", "bamr", "dplyr", "parallel")) %dopar% testfunction(filename)
#results
# stopCluster(cl)
end_time = Sys.time()
elapsed_time=end_time-start_time
print(elapsed_time)
#results_df = as.data.frame(t(results))
#results_df
# requirements to create a netcdf file
format_estimated_results(results)

