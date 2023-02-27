library(ncdf4)

# Preprocess for observational data
preproc_obs <- function(file,                                                 # Input file
						in_variable,                                                      # Input variable
						out_variable,                                                     # Output variable
						fmask,                                                            # Fixed mask of OBS
						obs.period,                                                       # Input period
						sel.period){                                                      # Output period
 # Masking
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
 obs_name <- cdo("chname",
		         args=paste0(in_variable,",",out_variable),
				     input=obs_prd,
				     options="-s -O")
 obs_in   <- cdo("mul",
		         input=c(obs_name,fmask),
				     options="-s -O")
 # Make 5-yr averaged time series of anomaly
 obs_rprd <- cdo("selyear",
		         args=paste0(ref.period[1],"/",ref.period[2]),
				     input=obs_in,
				     options="-s -O")
 obs_clim <- cdo("timmean",
		         input=obs_rprd,
				     options="-s -O")
 obs_ano  <- cdo("sub",
		         input=c(obs_in,obs_clim),
				     options="-s -O")
 obs_fld  <- cdo("fldmean",
		         input=obs_ano,
				     options="-s -O")
 obs_ave  <- cdo("timselmean",
		         args=5,
				     input=obs_fld,
				     options="-s -O")
 input_obs <- ncvar_get(nc_open(obs_ave),out_variable)
 unlink(c(obs_var,
          obs_cal,
          obs_tim,
          obs_prd,
          obs_name,
          obs_in,
          obs_rprd,
          obs_clim,
          obs_ano,
          obs_fld,
          obs_ave))
 return(input_obs)
}


# Preprocess for model data (Splitted data)
preproc_model_splitted <- function(
                       file1,                                                 # Input file1
			     				     file2,                                                 # Input file2
            				   variable,                                              # Input variable
            	   		   fmask,                                                 # Fixed mask of OBS
								       end.yr,                                                # Last year of file1
 						           sel.period,                                            # Output period
						           ref.period){                                           # Reference period
 # Preprocessing for first file
 file1_copy <- cdo("copy",
                   input=file1,
                   options="-s -O -f nc")
 file1_in   <- cdo("selyear",
                   args=paste0(sel.period[1],"/",end.yr),
                   input=file1_copy,
                   options="-s -O")
 # Preprocessing for second file
 file2_copy <- cdo("copy",
                   input=file2,
                   options="-s -O -f nc")
 file2_in   <- cdo("selyear",
                   args=paste0((end.yr+1),"/",sel.period[2]),
                   input=file2_copy,
                   options="-s -O")
 # Merge files (file1 + file2)
 file_merged <- cdo("mergetime",
                    input=c(file1_in,file2_in),
                    options="-s -O")
 merged_mask  <- cdo("remapbil",
                     args=fmask,
                     input=file_merged,
                     options="-s -O")
 merged_in    <- cdo("mul",
                     input=c(merged_mask,fmask),
                     options="-s -O")
 # Make 5-yr averaged time series of anomaly 
 merged_ref  <- cdo("selyear",
                    args=paste0(ref.period[1],"/",ref.period[2]),
                    input=merged_in,
                    options="-s -O")
 merged_clim <- cdo("timmean",
                    input=merged_ref,
                    options="-s -O")
 merged_ano  <- cdo("sub",
                    input=c(merged_in,merged_clim),
                    options="-s -O")
 merged_fld  <- cdo("fldmean",
                    input=merged_ano,
                    options="-s -O")
 merged_ave  <- cdo("timselmean,5", 
                    input=merged_fld,
                    options="-s -O")
 input_all   <- ncvar_get(nc_open(merged_ave),variable)
 unlink(c(file1_copy,file1_in,
          file2_copy,file2_in,
          file_merged,merged_mask,merged_in,
          merged_ref,merged_clim,merged_ano,merged_fld,merged_ave))
 return(input_all)
}



# Preprocess for model data (Non-splitted data)
preproc_model_nsplitted <- function(
                  file,                                                 # Input file
                  variable,                                             # Input variable
									fmask,                                                # Fixed mask of OBS
									sel.period,                                           # Output period
									ref.period){                                          # Reference period
 # Masking 
 file_copy <- cdo("copy",
                  input=file,
                  options="-s -O -f nc")
 file_prd  <- cdo("selyear",
                  args=paste0(sel.period[1],"/",sel.period[2]),
                  input=file_copy,
                  options="-s -O")
 file_mask <- cdo("remapbil",
                  args=fmask,
                  input=file_prd,
                  options="-s -O")
 file_in   <- cdo("mul",
                  input=c(file_mask,fmask),
                  options="-s -O")
 # Make 5-yr averaged time series of anomaly 
 file_ref  <- cdo("selyear",
                  args=paste0(ref.period[1],"/",ref.period[2]),
                  input=file_in,
                  options="-s -O")
 file_clim <- cdo("timmean",
                  input=file_ref,
                  options="-s -O")
 file_ano  <- cdo("sub",
                  input=c(file_in,file_clim),
                  options="-s -O")
 file_fld  <- cdo("fldmean",
                  input=file_ano,
                  options="-s -O")
 file_ave  <- cdo("timselmean,5", 
                  input=file_fld,
                  options="-s -O")
 input_file <- ncvar_get(nc_open(file_ave),variable)
 unlink(c(file_copy,file_prd,file_mask,file_in,file_ref,
		  file_clim,file_ano,file_fld,file_ave))
 return(input_file)
}
