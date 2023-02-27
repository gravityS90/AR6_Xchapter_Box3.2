library(ncdf4)

# Make fixed observational mask
fixed_mask <- function(
                       fin,                                                               # Input file
                       in_variable,                                                       # Input variable
                       out_variable,                                                      # Output variable
                       raw_period,                                                        # Input period
                       sel_period,                                                        # Output period
                       fou){                                                              # Output file
 # Input period	
 raw.syr <- raw_period[1]
 raw.eyr <- raw_period[2]
 raw.year <- c(raw.syr:raw.eyr)
 # Extract period
 syr <- sel_period[1]
 eyr <- sel_period[2]
 nyr <- eyr-syr+1
 loc.syr <- which(raw.year == syr)
 loc.eyr <- which(raw.year == eyr)
 
 # Read observation
 ncid <- nc_open(fin)
 val  <- ncvar_get(ncid,in_variable)
 lon  <- ncvar_get(ncid,"longitude")
 lat  <- ncvar_get(ncid,"latitude")
 time <- ncvar_get(ncid,"time")
 nx   <- ncid$dim[['longitude']]$len
 ny   <- ncid$dim[['latitude']]$len
 nt   <- ncid$dim[['time']]$len
 nc_close(ncid)

 # Fixed mask (more than 70% & existing at least 3 years out of last 5 years)
 mask <- matrix(NA,nx,ny)
 for(iy in c(1:ny)){
  for(ix in c(1:nx)){
   if( out_variable == "TXx" | out_variable == "txx" ){
    if(sum(!is.na(val[ix,iy,c(loc.syr:loc.eyr)])) >= ceiling(nyr*0.7) & 
       sum(!is.na(val[ix,iy,c((loc.eyr-4):loc.eyr)])) >= 3 ){
      mask[ix,iy] <- 1
	}
   }else if( out_variable == "Rx1day" | out_variable == "rx1day" ){
    if(sum(!is.na(val[ix,iy,c(loc.syr:loc.eyr)])) >= ceiling(nyr*0.7) & 
       sum(!is.na(val[ix,iy,c((loc.eyr-3):loc.eyr)])) >= 2 ){
      mask[ix,iy] <- 1
    } 
   } # if out_variable
  } # for ix
 } # for iy

 #--- Save fixed mask file 
 # Create new data ncfile
 xdim <- ncdim_def("lon","degrees_east",as.double(lon))
 ydim <- ncdim_def("lat","degrees_north",as.double(lat))
 # Set variables information
 fillvalue <- 1e32
 dlname    <- "Fixed mask"
 mask.def  <- ncvar_def("mask","none",list(xdim,ydim),fillvalue,dlname,prec="single")
 # Open nc file which is writed
 ncfname <- paste(fou,sep="")
 ncout   <- nc_create(ncfname,list(mask.def))
 # Assign variables
 ncvar_put(ncout,mask.def,mask)
 # Put additional attributes into dimension and data variables
 ncatt_put(ncout,"lon","axis","X")
 ncatt_put(ncout,"lat","axis","Y")
 nc_close(ncout)
}
