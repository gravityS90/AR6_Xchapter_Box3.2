library(ncdf4)
library(yaml)
library(extRemes)
library(ggplot2)
library(gridExtra)
source("./function/func_make_fixed_observational_mask.R")
source("./function/func_cdo.R")
source("./function/func_preproc.R")
source("./function/func_probability_index.R")


# Observation data
obs <- "HadEX3"

for( variable in c("TXx","Rx1day") ){
 # Set period
 obs.period <- c(1901,2018)
 ref.period <- c(1961,1990)                                                               # Reference period
 if(variable == "TXx")    sel.period.obs <- c(1953,2017)
 if(variable == "Rx1day") sel.period.obs <- c(1953,2011)
 sel.period     <- c(1953,2017)
 sel.period.nat <- c(1953,2012)                                                           # For only CMIP5

 # File path
 input_dir  <- paste0("./input/",obs,"/")
 output_dir <- "./output/"

 # Raw observation data
 fobs <- paste0(input_dir,
           obs,
           "_",
           variable,
           "_",
           obs.period[1],
           "-",
           obs.period[2],
 		       "_ADW_61-90_1.25x1.875deg.nc")
 
 # Make fixed mask
 fmask <- paste0(output_dir,
                 "Fixed_mask70_",
 				variable,
 				"_",
 				obs,
 				"_",
 				sel.period.obs[1],
 				"-",
 				sel.period.obs[2],
 				".nc")
 fixed_mask(fin=fobs,
            in_variable="Ann",
            out_variable=variable,
            raw_period=obs.period,
            sel_period=sel.period.obs,
            fou=fmask)
 
 # Read observation
 if( variable == "TXx" ){
  obs_value <- preproc_obs(file=fobs,
                           in_variable="Ann",
                           out_variable=variable,
                           fmask=fmask,
                           obs.period=obs.period,
                           sel.period=sel.period.obs)
 }else if( variable == "Rx1day" ){
  # CDF of GEVD for observation (for only Rx1day)
  fobs.pi <- paste0(input_dir,
           "cdf_",
 				   tolower(variable),
 				   "ETCCDI_",
 				   obs,
 				   "_",
 				   sel.period.obs[1],
 				   "-",
 				   sel.period.obs[2],
 				   ".nc")
  pi_obs(file=fobs,
         in_variable="Ann",
         unit="mm/day",                                                                   # Flux -> mutiply 86400
         fmask=fmask,
         obs.period=obs.period,
         sel.period=sel.period.obs,
         fou=fobs.pi)
  obs_value <- preproc_obs(file=fobs.pi,
                           in_variable="cdf",
                           out_variable="cdf",
                           fmask=fmask,
                           obs.period=obs.period,
                           sel.period=sel.period.obs)
  obs_value <- obs_value*100                                                              # Probability
 } # if variable
 

 # CMIP5 
 list0 <- yaml::read_yaml(paste0("./",variable,"_cmip5.yml"))
 cmip5_name <- unname(sapply(list0$datasets, "[[", "dataset"))
 cmip5_ensemble <- unname(sapply(list0$datasets, "[[", "ensemble"))
 cmip5_start_year <- unname(sapply(list0$datasets, "[[", "start_year"))
 cmip5_end_year <- unname(sapply(list0$datasets, "[[", "end_year"))
 cmip5_experiment <- unname(sapply(list0$datasets, "[[", "exp"))

 # CMIP5 file list (climdex data)
 template <- paste0(
     tolower(variable),
     "ETCCDI_yr_",
     cmip5_name,
     "_",
     cmip5_experiment,
     "_",
     cmip5_ensemble,
     "_",
     cmip5_start_year,
     "-",
     cmip5_end_year,
     ".nc")

 # CMIP5 file list (CDF of GEVD)
 if( variable == "Rx1day" ){
  cmip5.cdf.list <- function(x){
 	 switch(x,
          all=paste(
                    "cdf_",
                   	tolower(variable),
                    "ETCCDI_yr_",
                    cmip5_name[cmip5_experiment == "historical"],
                    "_",
                    cmip5_experiment[cmip5_experiment == "historical"],
                    "_",
                    cmip5_ensemble[cmip5_experiment == "historical"],
                    "_",
                    sel.period[1],
                    "-",
                    sel.period[2],
                    ".nc",
                    sep = ""),
         	nat=paste(
                  	"cdf_",
                  	tolower(variable),
                    "ETCCDI_yr_",
                    cmip5_name[cmip5_experiment == "historicalNat"],
                    "_",
                    cmip5_experiment[cmip5_experiment == "historicalNat"],
                    "_",
                    cmip5_ensemble[cmip5_experiment == "historicalNat"],
                    "_",
                    sel.period.nat[1],
                    "-",
                    sel.period.nat[2],
                    ".nc",
                    sep = ""))
  } # function
 } # if variable

 # CMIP5 file list
 work_dir <- "./input/CMIP5/"
 historical_file <- paste0(work_dir,template[grep("_historical_",template)])
 rcp45_file <- paste0(work_dir,template[grep("_rcp45_",template)])
 historicalNat_file <-  paste0(work_dir,template[grep("_historicalNat_",template)])

 # Fitting GEV for only Rx1day
 if( variable == "Rx1day" ){
  for( scenario in c("historical","historicalNat") ){
   climdex.input <- get(paste0(scenario,"_file"))
   cdf.climdex <- switch(scenario,
            		         historical=cmip5.cdf.list("all"),
             						 historicalNat=cmip5.cdf.list("nat"))
   for( model_idx in c(1:length(climdex.input)) ){
    if( !file.exists(paste0(work_dir,cdf.climdex)[model_idx]) ){                          # CDF file exist?
     switch(scenario,
    		    historical=pi_model_splitted(file1=climdex.input[model_idx],
 							                           file2=rcp45_file[model_idx],
                							           variable="rx1dayETCCDI",
                                         unit="mm/day",
                                         end.yr=2005,
 							                           sel.period=sel.period,
 							                           fou=paste0(work_dir,cdf.climdex[model_idx])),
            historicalNat=pi_model_nsplitted(
                                         file=climdex.input[model_idx],
             								             variable="rx1dayETCCDI",
                                         unit="mm/day",
 						             		             sel.period=sel.period.nat,
 								                         fou=paste0(work_dir,cdf.climdex[model_idx])))
    } # if 
   } # for model_idx
  } # for scenario
  historical_file <- paste0(work_dir,
                            switch("historical",
                                   historical=cmip5.cdf.list("all"),
                                   historicalNat=cmip5.cdf.list("nat")))
  historicalNat_file <- paste0(work_dir,
                               switch("historicalNat",
                                      historical=cmip5.cdf.list("all"),
                                      historicalNat=cmip5.cdf.list("nat")))
 } # if variable

 # 5-yr averaged time series of ALL (CMIP5)
 cmip5_each_all <- matrix(NA,length(historical_file),(sel.period[2]-sel.period[1]+1)/5)
 if( variable == "TXx" ){
  for(model_idx in c(1:length(historical_file))){
   cmip5_each_all[model_idx,] <- preproc_model_splitted(
                                   file1=historical_file[model_idx],
 		                               file2=rcp45_file[model_idx],
										               variable=paste0(tolower(variable),"ETCCDI"),
										               fmask=fmask,
                                   end.yr=2005,
										               sel.period=sel.period,
										               ref.period=ref.period)
  } # for model_idx
 }else if( variable == "Rx1day" ){
  for(model_idx in c(1:length(historical_file))){
   cmip5_each_all[model_idx,] <- preproc_model_nsplitted(
                                    file=historical_file[model_idx],
 		                                variable="cdf",
										                fmask=fmask,
                 								    sel.period=sel.period,
										                ref.period=ref.period)
  } # for model_idx
  cmip5_each_all <- cmip5_each_all*100                                                    # Probability
 } # if variable
 # Each model average
 cmip5_name_all <- unique(cmip5_name[grep("_historical_",template)])
 each_model_cmip5_each_all <- matrix(NA,length(cmip5_name_all),
  		                               (sel.period[2]-sel.period[1]+1)/5)
 for(model_idx in c(1:length(cmip5_name_all))){
  each_model_value <- cmip5_each_all[ grep(cmip5_name_all[model_idx], 
                                      grep("_historical_",template,value=T)) ,]
  if(!is.null(nrow(each_model_value))){
   each_model_cmip5_each_all[model_idx,] <- apply(each_model_value,2,mean)
  }else{
   each_model_cmip5_each_all[model_idx,] <- each_model_value
  } # if
 } # for model_idx

 # 5-yr averaged time series of NAT (CMIP5)
 cmip5_each_nat <- matrix(NA,length(historicalNat_file),
     		                  (sel.period.nat[2]-sel.period.nat[1]+1)/5)
 if( variable == "TXx" ){
  for(model_idx in c(1:length(historicalNat_file))){
   cmip5_each_nat[model_idx,] <- preproc_model_nsplitted(
                                    file=historicalNat_file[model_idx],
 		                                variable=paste0(tolower(variable),"ETCCDI"),
 										                fmask=fmask,
                 								    sel.period=sel.period.nat,
										                ref.period=ref.period)
  } # for model_idx
 }else if( variable == "Rx1day" ){
  for(model_idx in c(1:length(historicalNat_file))){
   cmip5_each_nat[model_idx,] <- preproc_model_nsplitted(
                                    file=historicalNat_file[model_idx],
 		                                variable="cdf",
										                fmask=fmask,
                 								    sel.period=sel.period.nat,
										                ref.period=ref.period)
  } # for model_dix
  cmip5_each_nat <- cmip5_each_nat*100                                                    # Probability
 } # if variable
 # Each model average
 cmip5_name_nat <- unique(cmip5_name[grep("_historicalNat_",template)])
 each_model_cmip5_each_nat <- matrix(NA,length(cmip5_name_nat),
 		                                 (sel.period.nat[2]-sel.period.nat[1]+1)/5)
 for(model_idx in c(1:length(cmip5_name_nat))){
  each_model_value <- cmip5_each_nat[ grep(cmip5_name_nat[model_idx],
 		                                  grep("_historicalNat_",template,value=T)) ,]
  if(!is.null(nrow(each_model_value))){
   each_model_cmip5_each_nat[model_idx,] <- apply(each_model_value,2,mean)
  }else{
   each_model_cmip5_each_nat[model_idx,] <- each_model_value
  } # if
 } # for model_idx

 # Multi model ensemble for CMIP5
 cmip5_mme_all <- apply(each_model_cmip5_each_all,2,mean)
 cmip5_mme_nat <- apply(each_model_cmip5_each_nat,2,mean)

 # Minumum and maximum range of each models
 cmip5_all_lb <- apply(each_model_cmip5_each_all,2,min)
 cmip5_all_ub <- apply(each_model_cmip5_each_all,2,max)
 cmip5_nat_lb <- apply(each_model_cmip5_each_nat,2,min)
 cmip5_nat_ub <- apply(each_model_cmip5_each_nat,2,max)


 # CMIP6 
 list0 <- yaml::read_yaml(paste0("./",variable,"_cmip6.yml"))
 cmip6_name <- unname(sapply(list0$datasets, "[[", "dataset"))
 cmip6_ensemble <- unname(sapply(list0$datasets, "[[", "ensemble"))
 cmip6_start_year <- unname(sapply(list0$datasets, "[[", "start_year"))
 cmip6_end_year <- unname(sapply(list0$datasets, "[[", "end_year"))
 cmip6_experiment <- unname(sapply(list0$datasets, "[[", "exp"))
 cmip6_grid <- unname(sapply(list0$datasets, "[[", "grid"))

 # CMIP6 file list (climdex data)
 template <- paste0(
  	tolower(variable),
    "ETCCDI_yr_",
    cmip6_name,
    "_",
    cmip6_experiment,
    "_",
    cmip6_ensemble,
    "_",
	  cmip6_grid,
 	"_",
    cmip6_start_year,
    "-",
    cmip6_end_year,
    ".nc")

 # CMIP6 file list (CDF of GEVD)
 cmip6.cdf.list <- function(x){
    switch(x,
           all=paste0(
                      "cdf_",
                      tolower(variable),
                      "ETCCDI_yr_",
                      cmip6_name[cmip6_experiment == "historical"],
                      "_",
                      cmip6_experiment[cmip6_experiment == "historical"],
                      "_",
                      cmip6_ensemble[cmip6_experiment == "historical"],
                      "_",
                      cmip6_grid[cmip6_experiment == "historical"],
                      "_",
                      sel.period[1],
                      "-",
                      sel.period[2],
                      ".nc"),
           nat=paste0(
                      "cdf_",
                      tolower(variable),
                      "ETCCDI_yr_",
                      cmip6_name[cmip6_experiment == "hist-nat"],
                      "_",
                      cmip6_experiment[cmip6_experiment == "hist-nat"],
                      "_",
                      cmip6_ensemble[cmip6_experiment == "hist-nat"],
                      "_",
                     	cmip6_grid[cmip6_experiment == "hist-nat"],
                    	"_",
                      sel.period.nat[1],
                      "-",
                      sel.period.nat[2],
                      ".nc"))
 } # function

 # CMIP6 file list
 work_dir <- "./input/CMIP6/"
 historical_file <- paste0(work_dir,template[grep("_historical_",template)])
 ssp245_file <- paste0(work_dir,template[grep("_ssp245_",template)])
historicalNat_file <-  paste0(work_dir,template[grep("_hist-nat_",template)])
 
 # Fitting GEV for only Rx1day
 if( variable == "Rx1day" ){
  for( scenario in c("historical","historicalNat") ){
   climdex.input <- get(paste0(scenario,"_file"))
   cdf.climdex <- switch(scenario,
                         historical=cmip6.cdf.list("all"),
                         historicalNat=cmip6.cdf.list("nat"))
   for( model_idx in c(1:length(climdex.input)) ){
    if( !file.exists(paste0(work_dir,cdf.climdex)[model_idx]) ){                          # CDF file exist?
     switch(scenario,
            historical=pi_model_splitted(file1=climdex.input[model_idx],
                                         file2=ssp245_file[model_idx],
                                         variable="rx1dayETCCDI",
                                         unit="mm/day",                                   # Flux -> multiply 86400
                                         end.yr=2014,
                                         sel.period=sel.period,
                                         fou=paste0(work_dir,cdf.climdex[model_idx])),
            historicalNat=pi_model_nsplitted(file=climdex.input[model_idx],
                                             variable="rx1dayETCCDI",
                                             unit="mm/day",                               # Flux -> multiply 86400
                                             sel.period=sel.period,
                                             fou=paste0(work_dir,cdf.climdex[model_idx])))

    } # if
   } # for model_idx
  } # for scenario
  historical_file <- paste0(work_dir,switch("historical",
                                             historical=cmip6.cdf.list("all"),
                                             historicalNat=cmip6.cdf.list("nat")))
  historicalNat_file <- paste0(work_dir,switch("historicalNat",
                                               historical=cmip6.cdf.list("all"),
                                               historicalNat=cmip6.cdf.list("nat")))
 } # if variable

 # 5-yr averaged time series of ALL (CMIP6)
 cmip6_each_all <- matrix(NA,length(historical_file),(sel.period[2]-sel.period[1]+1)/5)
 if( variable == "TXx" ){
  for(model_idx in c(1:length(historical_file))){
   cmip6_each_all[model_idx,] <- preproc_model_splitted(file1=historical_file[model_idx],
                                                        file2=ssp245_file[model_idx],
                                                        variable=paste0(tolower(variable),"ETCCDI"),
                                                        fmask=fmask,
                                                        end.yr=2014,
                                                        sel.period=sel.period,
                                                        ref.period=ref.period)
  } # for model_idx
 }else if( variable == "Rx1day" ){
  for(model_idx in c(1:length(historical_file))){
   cmip6_each_all[model_idx,] <- preproc_model_nsplitted(file=historical_file[model_idx],
                                                         variable="cdf",
                                                         fmask=fmask,
                                                         sel.period=sel.period,
                                                         ref.period=ref.period)
  } # for model_idx
  cmip6_each_all <- cmip6_each_all*100                                                     # Probability
 } # if variable
 # Each model average
 cmip6_name_all <- unique(cmip6_name[grep("_historical_",template)])
 each_model_cmip6_each_all <- matrix(NA,length(cmip6_name_all),
      		                           (sel.period[2]-sel.period[1]+1)/5)
 for(model_idx in c(1:length(cmip6_name_all))){
  each_model_value <- cmip6_each_all[ grep(cmip6_name_all[model_idx],
 		                                  grep("_historical_",template,value=T)) ,]
  if(!is.null(nrow(each_model_value))){
   each_model_cmip6_each_all[model_idx,] <- apply(each_model_value,2,mean)
  }else{
   each_model_cmip6_each_all[model_idx,] <- each_model_value
  } # if
 } # for model_idx

 # 5-yr averaged time series of NAT (CMIP6)
 cmip6_each_nat <- matrix(NA,length(historicalNat_file),
                          (sel.period[2]-sel.period[1]+1)/5)
 if( variable == "TXx" ){
  for(model_idx in c(1:length(historicalNat_file))){
   cmip6_each_nat[model_idx,] <- preproc_model_nsplitted(historicalNat_file[model_idx],
                                                         paste0(tolower(variable),"ETCCDI"),
                                                         fmask,
                                                         sel.period,
                                                         ref.period)
  } # for model_idx
 }else if( variable == "Rx1day" ){
  for(model_idx in c(1:length(historicalNat_file))){
   cmip6_each_nat[model_idx,] <- preproc_model_nsplitted(historicalNat_file[model_idx],
                                                         "cdf",
                                                         fmask,
                                                         sel.period,
                                                         ref.period)
  } # for model_idx
  cmip6_each_nat <- cmip6_each_nat*100                                                    # Probability
 } # if
 # Each model average
 cmip6_name_nat <- unique(cmip6_name[grep("_hist-nat_",template)])
 each_model_cmip6_each_nat <- matrix(NA,length(cmip6_name_nat),
 		                                 (sel.period[2]-sel.period[1]+1)/5)
 for(model_idx in c(1:length(cmip6_name_nat))){
  each_model_value <- cmip6_each_nat[ grep(cmip6_name_nat[model_idx],
 		                                  grep("_hist-nat_",template,value=T)) ,]
  if(!is.null(nrow(each_model_value))){
   each_model_cmip6_each_nat[model_idx,] <- apply(each_model_value,2,mean)
  }else{
   each_model_cmip6_each_nat[model_idx,] <- each_model_value
  } # if
 } # for model_idx

 # Multi model ensemble for cmip6
 cmip6_mme_all <- apply(each_model_cmip6_each_all,2,mean)
 cmip6_mme_nat <- apply(each_model_cmip6_each_nat,2,mean)

 if(variable == "TXx"){
  # Data frame for ALL
  dummy_all_list <- paste0("A",formatC(c(1:length(cmip6_name_all)),flag="0",width=2))
  cmip5_all_txx_rng  <- data.frame(yr=seq(1955,2015,5),
                                   lb=cmip5_all_lb,
                                   ub=cmip5_all_ub)
  cmip6_all_txx_each <- data.frame(model=rep(dummy_all_list,each=(sel.period[2]-sel.period[1]+1)/5),
 		                               yr=seq(1955,2015,5),
                                   y=as.numeric(t(each_model_cmip6_each_all)))
  txx_obs_and_mme_all <- data.frame(project=rep(LETTERS[1:3],each=(sel.period[2]-sel.period[1]+1)/5), 
 		                                yr=seq(1955,2015,5),
 									                  y=c(cmip6_mme_all,cmip5_mme_all,obs_value))
  # Data frame for NAT
  dummy_nat_list <- paste0("A",formatC(c(1:length(cmip6_name_nat)),flag="0",width=2))
  cmip5_nat_txx_rng  <- data.frame(yr=seq(1955,2015,5),
                                   lb=c(cmip5_nat_lb,NA),
                                   ub=c(cmip5_nat_ub,NA))
  cmip6_nat_txx_each <- data.frame(model=rep(dummy_nat_list,each=(sel.period[2]-sel.period[1]+1)/5),
 		                               yr=seq(1955,2015,5),
                                   y=as.numeric(t(each_model_cmip6_each_nat)))
  txx_obs_and_mme_nat <- data.frame(project=rep(LETTERS[1:3],each=(sel.period[2]-sel.period[1]+1)/5), 
		                                yr=seq(1955,2015,5),
  									                y=c(cmip6_mme_nat,c(cmip5_mme_nat,NA),obs_value))
  txx_ncmip5_all <- length(cmip5_name_all)
  txx_ncmip5_nat <- length(cmip5_name_nat)
  txx_ncmip6_all <- length(cmip6_name_all)
  txx_ncmip6_nat <- length(cmip6_name_nat)
 }else if(variable == "Rx1day"){
  # Data frame for ALL
  dummy_all_list <- paste0("A",formatC(c(1:length(cmip6_name_all)),flag="0",width=2))
  cmip5_all_rx1day_rng  <- data.frame(yr=seq(1955,2015,5),
                                      lb=cmip5_all_lb,
                                      ub=cmip5_all_ub)
  cmip6_all_rx1day_each <- data.frame(model=rep(dummy_all_list,each=(sel.period[2]-sel.period[1]+1)/5),
 		                                  yr=seq(1955,2015,5),
                                      y=as.numeric(t(each_model_cmip6_each_all)))
  rx1day_obs_and_mme_all <- data.frame(project=rep(LETTERS[1:3],each=(sel.period[2]-sel.period[1]+1)/5), 
 		                                   yr=seq(1955,2015,5),
 									                     y=c(cmip6_mme_all,cmip5_mme_all,c(obs_value,NA)))
  # Data frame for NAT
  dummy_nat_list <- paste0("A",formatC(c(1:length(cmip6_name_nat)),flag="0",width=2))
  cmip5_nat_rx1day_rng  <- data.frame(yr=seq(1955,2015,5),
                                      lb=c(cmip5_nat_lb,NA),
                                      ub=c(cmip5_nat_ub,NA))
  cmip6_nat_rx1day_each <- data.frame(model=rep(dummy_nat_list,each=(sel.period[2]-sel.period[1]+1)/5),
 		                                  yr=seq(1955,2015,5),
                                      y=as.numeric(t(each_model_cmip6_each_nat)))
  rx1day_obs_and_mme_nat <- data.frame(project=rep(LETTERS[1:3],each=(sel.period[2]-sel.period[1]+1)/5), 
 		                                   yr=seq(1955,2015,5),
 									                     y=c(cmip6_mme_nat,c(cmip5_mme_nat,NA),c(obs_value,NA)))
  rx1day_ncmip5_all <- length(cmip5_name_all)
  rx1day_ncmip5_nat <- length(cmip5_name_nat)
  rx1day_ncmip6_all <- length(cmip6_name_all)
  rx1day_ncmip6_nat <- length(cmip6_name_nat)
 } # if variable
} # for variable

# Plots
common <- theme_bw() +
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank()) +
          theme(plot.margin=unit(c(1.5,1,1.0,0.5),"cm")) +
          theme(axis.title.y=element_text(margin=margin(0,20,0,0)),
                axis.title.x=element_text(margin=margin(20,0,0,0)),
                axis.title=element_text(size=18, face="bold"),
                axis.text=element_text(size=18, face="bold")) +
          theme(legend.position=c(0.18,0.82),
                legend.text=element_text(size=15),
                legend.key.height=unit(1.5,"line"),
                legend.key.width=unit(3,"line"),
                legend.background=element_rect(colour="transparent",fill="transparent"),
                legend.key=element_rect(fill="transparent",colour="transparent"),
                legend.direction="vertical",
                legend.title=element_blank()) 

p_txx_all <- ggplot() + 
             geom_ribbon(data=cmip5_all_txx_rng,aes(x=yr,ymin=lb,ymax=ub),
                         fill=rgb(254/255,229/255,217/255)) +
             geom_hline(yintercept=0,linetype=2,colour="grey50",size=1) +
             geom_line(data=cmip6_all_txx_each, aes(x=yr,y=y,group=model),
					   colour=rgb(255/255,191/255,128/255),size=0.5) +
#             geom_line(data=txx_obs_and_mme_all,aes(x=yr,y=y,colour=project),size=1.5,na.rm=T) +
             geom_line(data=txx_obs_and_mme_all[txx_obs_and_mme_all$project == "B",],
					   aes(x=yr,y=y,colour="B"),size=1.5,na.rm=T) +
             geom_line(data=txx_obs_and_mme_all[txx_obs_and_mme_all$project == "A",],
					   aes(x=yr,y=y,colour="A"),size=1.5,na.rm=T) +
             geom_line(data=txx_obs_and_mme_all[txx_obs_and_mme_all$project == "C",],
					   aes(x=yr,y=y,colour="C"),size=1.5,na.rm=T) +
             scale_color_manual(values=c(rgb(180/255,90/255,0/255),rgb(222/255,45/255,38/255),"black"),
                                labels=c(paste("CMIP6 (",txx_ncmip6_all,")",sep=""),
                                         paste("CMIP5 (",txx_ncmip5_all,")",sep=""),
                                         "Observation"),guide=guide_legend(reverse = TRUE)) +
             scale_x_continuous(limits=c(1950,2020),breaks=c(seq(1950,2020,10))) +
             scale_y_continuous(limit=c(-1.2,2.4),breaks=seq(-1,2,1)) +
             xlab("Year") +
             ylab(expression(bold("TXx ("*~degree*C*")"))) + 
             common + 
     			   annotate("text",x=1969,y=2.4,size=7,label="Natural and Human forcing")

p_txx_nat <- ggplot() + 
             geom_ribbon(data=cmip5_nat_txx_rng,aes(x=yr,ymin=lb,ymax=ub),
                         fill=rgb(254/255,229/255,217/255)) +
             geom_hline(yintercept=0,linetype=2,colour="grey50",size=1) +
             geom_line(data=cmip6_nat_txx_each, aes(x=yr,y=y,group=model),
                       colour=rgb(255/255,191/255,128/255),size=0.5) +
#             geom_line(data=txx_obs_and_mme_nat,aes(x=yr,y=y,colour=project),size=1.5,na.rm=T) +
             geom_line(data=txx_obs_and_mme_nat[txx_obs_and_mme_nat$project == "B",],
           					   aes(x=yr,y=y,colour="B"),size=1.5,na.rm=T) +
             geom_line(data=txx_obs_and_mme_nat[txx_obs_and_mme_nat$project == "A",],
					             aes(x=yr,y=y,colour="A"),size=1.5,na.rm=T) +
             geom_line(data=txx_obs_and_mme_nat[txx_obs_and_mme_nat$project == "C",],
            				   aes(x=yr,y=y,colour="C"),size=1.5,na.rm=T) +
             scale_color_manual(values=c(rgb(180/255,90/255,0/255),rgb(222/255,45/255,38/255),"black"),
                                labels=c(paste("CMIP6 (",txx_ncmip6_nat,")",sep=""),
                                         paste("CMIP5 (",txx_ncmip5_nat,")",sep=""),
                                         "Observation"),guide=guide_legend(reverse = TRUE)) +
             scale_x_continuous(limits=c(1950,2020),breaks=c(seq(1950,2020,10))) +
             scale_y_continuous(limit=c(-1.2,2.4),breaks=seq(-1,2,1)) +
             xlab("Year") +
             ylab(expression(bold("TXx ("*~degree*C*")"))) + 
             common + 
			       annotate("text",x=1961,y=2.4,size=7,label="Natural forcing")

p_rx1day_all <- ggplot() + 
                geom_ribbon(data=cmip5_all_rx1day_rng,aes(x=yr,ymin=lb,ymax=ub),
                            fill=rgb(224/255,236/255,244/255)) +
                geom_hline(yintercept=0,linetype=2,colour="grey50",size=1) +
                geom_line(data=cmip6_all_rx1day_each, aes(x=yr,y=y,group=model),
					                colour=rgb(146/255,236/255,148/255),size=0.5) +
#                geom_line(data=rx1day_obs_and_mme_all,aes(x=yr,y=y,colour=project),size=1.5,na.rm=T) +
                geom_line(data=rx1day_obs_and_mme_all[rx1day_obs_and_mme_all$project == "B",],
	    			              aes(x=yr,y=y,colour="B"),size=1.5,na.rm=T) +
                geom_line(data=rx1day_obs_and_mme_all[rx1day_obs_and_mme_all$project == "A",],
		    		              aes(x=yr,y=y,colour="A"),size=1.5,na.rm=T) +
                geom_line(data=rx1day_obs_and_mme_all[rx1day_obs_and_mme_all$project == "C",],
			    	              aes(x=yr,y=y,colour="C"),size=1.5,na.rm=T) +
                scale_color_manual(values=c(rgb(116/255,196/255,118/255),rgb(40/255,130/255,189/255),"black"),
                                labels=c(paste("CMIP6 (",rx1day_ncmip6_all,")",sep=""),
                                         paste("CMIP5 (",rx1day_ncmip5_all,")",sep=""),
                                         "Observation"),guide=guide_legend(reverse = TRUE)) +
                scale_x_continuous(limits=c(1950,2020),breaks=c(seq(1950,2020,10))) +
                scale_y_continuous(limit=c(-5,10),breaks=seq(-5,10,5)) +
                xlab("Year") +
                ylab("Rx1day (%)") +
             		common + 
                annaotate("text",x=1969,y=10,size=7,label="Natural and Human forcing")

p_rx1day_nat <- ggplot() + 
                geom_ribbon(data=cmip5_nat_rx1day_rng,aes(x=yr,ymin=lb,ymax=ub),
                            fill=rgb(224/255,236/255,244/255)) +
                geom_hline(yintercept=0,linetype=2,colour="grey50",size=1) +
                geom_line(data=cmip6_nat_rx1day_each, aes(x=yr,y=y,group=model),
   					              colour=rgb(146/255,236/255,148/255),size=0.5) +
#                geom_line(data=rx1day_obs_and_mme_nat,aes(x=yr,y=y,colour=project),size=1.5,na.rm=T) +
                geom_line(data=rx1day_obs_and_mme_nat[rx1day_obs_and_mme_nat$project == "B",],
   			            		  aes(x=yr,y=y,colour="B"),size=1.5,na.rm=T) +
                geom_line(data=rx1day_obs_and_mme_nat[rx1day_obs_and_mme_nat$project == "A",],
   					              aes(x=yr,y=y,colour="A"),size=1.5,na.rm=T) +
                geom_line(data=rx1day_obs_and_mme_nat[rx1day_obs_and_mme_nat$project == "C",],
   					              aes(x=yr,y=y,colour="C"),size=1.5,na.rm=T) +
                scale_color_manual(values=c(rgb(116/255,196/255,118/255),rgb(40/255,130/255,189/255),"black"),
                                   labels=c(paste("CMIP6 (",rx1day_ncmip6_nat,")",sep=""),
                                            paste("CMIP5 (",rx1day_ncmip5_nat,")",sep=""),
                                            "Observation"),guide=guide_legend(reverse = TRUE)) +
                scale_x_continuous(limits=c(1950,2020),breaks=c(seq(1950,2020,10))) +
                scale_y_continuous(limit=c(-5,10),breaks=seq(-5,10,5)) +
                xlab("Year") +
                ylab("Rx1day (%)") +
         			  common + 
                annotate("text",x=1961,y=10,size=7,label="Natural forcing")

plot_dir <- paste0("./Figures/")
png(paste0(plot_dir,"Fig1_Xchapter_box3.2.png"),width=1200,height=1000)
grid.arrange(p_txx_all,p_rx1day_all,p_txx_nat,p_rx1day_nat,ncol=2)
dev.off()
