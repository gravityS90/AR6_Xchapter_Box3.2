#library(ismev)
library(ncdf4)
library(extRemes)
library(evd)

# Turn off warnig messages
options(warn=-1)

# CDF of GEVD for observation	
pi_obs <- function(
                   file,                                                                          # Input file
                   in_variable,                                                                   # Input variable
                   unit,                                                                          # Input data unit
                   fmask,                                                                         # Fixed mask of OBS
                   obs.period,                                                                    # Input period
                   sel.period,                                                                    # Output period
                   fou){                                                                          # Output file
 obs_var  <- cdo("selvar",
                 args=in_variable,
                 input=file,
                 options="-s -O")
 obs_cal  <- cdo("setcalendar",
                 args="standard",
                 input=obs_var,
                 options="-s -O")
 obs_tim  <- cdo("setreftime",
                 args=paste0(obs.period[1],"-07-16,00:00:00"),
                 input=obs_cal,
                 options="-s -O")
 obs_prd  <- cdo("selyear",
                 args=paste0(sel.period[1],"/",sel.period[2]),
                 input=obs_tim,
                 options="-s -O")
 obs_in   <- cdo("mul",
                 input=c(obs_prd,fmask),
                 options="-s -O")

 # GEV fitting period (1953-2011)
 loc1 <- which(c(sel.period[1]:sel.period[2]) == 1953)
 loc2 <- which(c(sel.period[1]:sel.period[2]) == 2011)

 # Read input file
 ncid <- nc_open(obs_in)
 val  <- ncvar_get(ncid,in_variable)
 lon  <- ncvar_get(ncid,"lon")
 lat  <- ncvar_get(ncid,"lat")
 time <- ncvar_get(ncid,"time")
 nx   <- ncid$dim[['lon']]$len
 ny   <- ncid$dim[['lat']]$len
 nt   <- ncid$dim[['time']]$len
 nc_close(ncid)

 if(unit == "flux") val <- val*86400                                                      # Flux to mm

 # Estimate cumulative density function of GEVD
 cdf  <- array(NA,c(nx,ny,nt))
 for(iy in c(1:ny)){
  for(ix in c(1:nx)){
   if( all(is.na(val[ix,iy,c(loc1:loc2)])) |
       sum(val[ix,iy,c(loc1:loc2)] <= 1e-04 ,na.rm=T) > nt*0.35 ){                        # Exclude dry grid
   }else{
	tmp <- val[ix,iy,(loc1:loc2)]
#    gev_params <- ismev::gev.fit(tmp[!is.na(tmp)],show=F)$mle                            # GEV fitting
    gev_params <- extRemes::fevd(tmp[!is.na(tmp)])$results$par                            # GEV fitting
    tmp_cdf <- evd:::pgev(val[ix,iy,],
                          loc=gev_params[1],
                          scale=gev_params[2],
                          shape=gev_params[3])                                            # CDF of GEVD
    tmp_cdf[tmp_cdf == 1] <- NA
    cdf[ix,iy,] <- tmp_cdf
   } # if
  } # for ix
 } # for iy

 # Create new data ncfile
 xdim <- ncdim_def("lon","degrees_east",as.double(lon))
 ydim <- ncdim_def("lat","degrees_north",as.double(lat))
 tunits3 <- paste("years since ",
		  sel.period[1],
		  "-07-16 12:00:00.0-0:00", sep="")
 tdim    <- ncdim_def("time",tunits3,as.double(c(0:(nt-1))))
 # Set variables information
 fillvalue <- 1e32
 dlname    <- "Probability index"
 pi.def  <- ncvar_def("cdf","none",
                      list(xdim,ydim,tdim),
                      fillvalue,
                      dlname,
                      prec="single")
 # Open nc file which is writed
 ncfname <- paste(fou,sep="")
 ncout   <- nc_create(ncfname,list(pi.def))
 # Assign variables
 ncvar_put(ncout,pi.def,cdf)
 # Put additional attributes into dimension and data variables
 ncatt_put(ncout,"lon","axis","X")
 ncatt_put(ncout,"lat","axis","Y")
 ncatt_put(ncout,"time","axis","T")
 nc_close(ncout)
 unlink(c(obs_var,obs_cal,obs_tim,obs_prd,obs_in))
} # function



# CDF of GEVD for splitted model data
pi_model_splitted <- function(
                              file1,                                                                    # Input file1
                              file2,                                                                    # Input file2
                              variable,                                                                 # Input variable
                              unit,                                                                     # Input data unit
                              end.yr,                                                                   # Last year of file1
                              sel.period,                                                               # Output period
                              fou){                                                                     # Output file
 # For first file
 file1_copy <- cdo("copy",
                   input=file1,
                   options="-s -O -f nc")
 file1_in   <- cdo("selyear",
                   args=paste0(sel.period[1],"/",end.yr),
                   input=file1_copy,
                   options="-s -O")
 # For second file
 file2_copy <- cdo("copy",
                   input=file2,
                   options="-s -O -f nc")
 file2_in   <- cdo("selyear",
                   args=paste0((end.yr+1),"/",sel.period[2]),
                   input=file2_copy,
                   options="-s -O")
 # Merge files (file1 + file2)
 merged_in <- cdo("mergetime",
                  input=c(file1_in,file2_in),
                  options="-s -O")

 # GEV fitting period (1953-2011)
 loc1 <- which(c(sel.period[1]:sel.period[2]) == 1953)
 loc2 <- which(c(sel.period[1]:sel.period[2]) == 2011)

 # Read input file
 ncid <- nc_open(merged_in)
 val  <- ncvar_get(ncid,variable)
 lon  <- ncvar_get(ncid,"lon")
 lat  <- ncvar_get(ncid,"lat")
 time <- ncvar_get(ncid,"time")
 nx   <- ncid$dim[['lon']]$len
 ny   <- ncid$dim[['lat']]$len
 nt   <- ncid$dim[['time']]$len
 nc_close(ncid)

 if(unit == "flux") val <- val*86400                                                      # Flux to mm

 # Estimate cumulative density function of GEVD
 cdf  <- array(NA,c(nx,ny,nt))
 for(iy in c(1:ny)){
  for(ix in c(1:nx)){
   if( all(is.na(val[ix,iy,c(loc1:loc2)])) |
       sum(val[ix,iy,c(loc1:loc2)] <= 1e-04 ,na.rm=T) > nt*0.35 ){                        # Exclude dry grid
   }else{
	tmp <- val[ix,iy,(loc1:loc2)]
#    gev_params <- ismev::gev.fit(tmp[!is.na(tmp)],show=F)$mle                            # GEV fitting
    gev_params <- extRemes::fevd(tmp[!is.na(tmp)])$results$par                            # GEV fitting
    tmp_cdf <- evd:::pgev(val[ix,iy,],
                          loc=gev_params[1],
                          scale=gev_params[2],
                          shape=gev_params[3])                                            # CDF of GEVD
    tmp_cdf[tmp_cdf == 1] <- NA
    cdf[ix,iy,] <- tmp_cdf
   } # if
  } # for ix
 } # for iy

 # Create new data ncfile
 xdim <- ncdim_def("lon","degrees_east",as.double(lon))
 ydim <- ncdim_def("lat","degrees_north",as.double(lat))
 tunits3 <- paste("years since ",
                  sel.period[1],
                  "-07-16 12:00:00.0-0:00", sep="")
 tdim    <- ncdim_def("time",tunits3,as.double(c(0:(nt-1))))
 # Set variables information
 fillvalue <- 1e32
 dlname    <- "Probability index"
 pi.def  <- ncvar_def("cdf","none",
                      list(xdim,ydim,tdim),
                      fillvalue,
                      dlname,prec="single")
 # Open nc file which is writed
 ncfname <- paste(fou,sep="")
 ncout   <- nc_create(ncfname,list(pi.def))
 # Assign variables
 ncvar_put(ncout,pi.def,cdf)
 # Put additional attributes into dimension and data variables
 ncatt_put(ncout,"lon","axis","X")
 ncatt_put(ncout,"lat","axis","Y")
 ncatt_put(ncout,"time","axis","T")
 nc_close(ncout)
 unlink(c(file1_copy,file1_in,file2_copy,file2_in,merged_in))
} # function


# CDF of GEVD for non-splitted model data
pi_model_nsplitted <- function(
                               file,                                                                    # Input file
                               variable,                                                                # Input variable
                               unit,                                                                    # Input data unit
                               sel.period,                                                              # Output period
                               fou){                                                                    # Output file
 # Preprocessing
 model_copy <- cdo("copy",
                   input=file,
                   options="-s -O -f nc")
 model_in   <- cdo("selyear",
                   args=paste0(sel.period[1],"/",sel.period[2]),
                   input=model_copy,
                   options="-s -O")

 # GEV fitting period (1953-2011)
 loc1 <- which(c(sel.period[1]:sel.period[2]) == 1953)
 loc2 <- which(c(sel.period[1]:sel.period[2]) == 2011)

 # Read input file
 ncid <- nc_open(model_in)
 val  <- ncvar_get(ncid,variable)
 lon  <- ncvar_get(ncid,"lon")
 lat  <- ncvar_get(ncid,"lat")
 time <- ncvar_get(ncid,"time")
 nx   <- ncid$dim[['lon']]$len
 ny   <- ncid$dim[['lat']]$len
 nt   <- ncid$dim[['time']]$len
 nc_close(ncid)

 if(unit == "flux") val <- val*86400                                                      # Flux to mm

 # Estimate cumulative density function of GEVD
 cdf  <- array(NA,c(nx,ny,nt))
 for(iy in c(1:ny)){
  for(ix in c(1:nx)){
   if( all(is.na(val[ix,iy,c(loc1:loc2)])) |
       sum(val[ix,iy,c(loc1:loc2)] <= 1e-04 ,na.rm=T) > nt*0.35 ){                        # Exclude dry grid
   }else{
	tmp <- val[ix,iy,(loc1:loc2)]
#    gev_params <- ismev::gev.fit(tmp[!is.na(tmp)],show=F)$mle                            # GEV fitting
    gev_params <- extRemes::fevd(tmp[!is.na(tmp)])$results$par                            # GEV fitting
    tmp_cdf <- evd:::pgev(val[ix,iy,],
                          loc=gev_params[1],
                          scale=gev_params[2],
                          shape=gev_params[3])                                            # CDF of GEVD
    tmp_cdf[tmp_cdf == 1] <- NA
    cdf[ix,iy,] <- tmp_cdf
   } # if
  } # for ix
 } # for iy

 # Create new data ncfile
 xdim <- ncdim_def("lon","degrees_east",as.double(lon))
 ydim <- ncdim_def("lat","degrees_north",as.double(lat))
 tunits3 <- paste("years since ",
                  sel.period[1],
                  "-07-16 12:00:00.0-0:00", sep="")
 tdim    <- ncdim_def("time",tunits3,as.double(c(0:(nt-1))))
 # Set variables information
 fillvalue <- 1e32
 dlname    <- "Probability index"
 pi.def  <- ncvar_def("cdf","none",
                      list(xdim,ydim,tdim),
                      fillvalue,
                      dlname,prec="single")
 # Open nc file which is writed
 ncfname <- paste(fou,sep="")
 ncout   <- nc_create(ncfname,list(pi.def))
 # Assign variables
 ncvar_put(ncout,pi.def,cdf)
 # Put additional attributes into dimension and data variables
 ncatt_put(ncout,"lon","axis","X")
 ncatt_put(ncout,"lat","axis","Y")
 ncatt_put(ncout,"time","axis","T")
 nc_close(ncout)
 unlink(c(model_copy,model_in))
} # function

# Turn on warnig messages
options(warn=0)
