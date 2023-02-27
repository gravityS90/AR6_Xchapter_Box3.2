library(yaml)

# Select variable (TXx or Rx1day)
variable   <- "TXx" 

#=========================================================================================
#                                  1) Copy CMIP6 models
#=========================================================================================
# Obtain CMIP6 model list from yml
list0 <- yaml::read_yaml(paste0("./",variable,"_cmip6.yml"))
cmip6_name <- unname(sapply(list0$datasets, "[[", "dataset"))                             # Model 
cmip6_ensemble <- unname(sapply(list0$datasets, "[[", "ensemble"))                        # Ensemble
cmip6_start_year <- unname(sapply(list0$datasets, "[[", "start_year"))                    # Start year
cmip6_end_year <- unname(sapply(list0$datasets, "[[", "end_year"))                        # End year
cmip6_experiment <- unname(sapply(list0$datasets, "[[", "exp"))                           # Experiment
cmip6_grid <- unname(sapply(list0$datasets, "[[", "grid"))                                # Grid

# CMIP6 file list
template <- paste(
    paste0(tolower(variable),"ETCCDI"),
    "_yr_",
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
    ".nc",
    sep = "")

pou <- "/data2/mgseong/2020/00_AR6_Figs_code/_r/Final/input/CMIP6/"
for(sce in c("historical","ssp245","hist-nat")){
 flist <- template[grep(paste0("_",sce,"_"),template)]
 nfile <- length(flist)

 for(ifile in c(1:nfile)){
  modl <- strsplit(flist[ifile],"_")[[1]][3]
  ense <- strsplit(flist[ifile],"_")[[1]][5]
  grd  <- strsplit(flist[ifile],"_")[[1]][6]
  pin  <- paste0("/data2/mgseong/CMIP6/ann/TXx/",sce,"/",modl,"/")
  system(paste0("find ",pin," -name '",tolower(variable),
                "ETCCDI_yr*",ense,"*",grd,"*' -exec cp {} ", pou," \\;"),intern=F)
 } # for ifile
} # for sce

#=========================================================================================
#                                  2) Copy CMIP5 models
#=========================================================================================
# Obtain CMIP5 model list from yml
list0 <- yaml::read_yaml("./Rx1day_cmip5.yml")
cmip5_name <- unname(sapply(list0$datasets, "[[", "dataset"))
cmip5_ensemble <- unname(sapply(list0$datasets, "[[", "ensemble"))
cmip5_start_year <- unname(sapply(list0$datasets, "[[", "start_year"))
cmip5_end_year <- unname(sapply(list0$datasets, "[[", "end_year"))
cmip5_experiment <- unname(sapply(list0$datasets, "[[", "exp"))

# CMIP5 file list
template <- paste(
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
    ".nc",
    sep = "")

pou <- "/data2/mgseong/2020/00_AR6_Figs_code/_r/Final/input/CMIP5/rx1day/"
for(sce in c("historical","rcp45","historicalNat")){
 flist <- template[grep(paste0("_",sce,"_"),template)]
 nfile <- length(flist)

 for(ifile in c(1:nfile)){
  modl <- strsplit(flist[ifile],"_")[[1]][3]
  ense <- strsplit(flist[ifile],"_")[[1]][5]
  pin  <- paste0("/disk4/climdex/CMIP5_new/",sce,"/",modl,"/",ense,"/")
  system(paste0("find ",pin," -name '",tolower(variable),
                "ETCCDI_yr*' -exec cp {} ", pou," \\;"),intern=F)
 } # for ifile
} # for sce
