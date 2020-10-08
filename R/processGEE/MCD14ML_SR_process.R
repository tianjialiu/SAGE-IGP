# --------------------------------------------
# Pre-process MCD14ML SR files from EE
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# --------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')

# backfill empty files
setwd(file.path(geeHome,"MCD14ML_SR"))

inState <- "Punjab"
infiles <- do.call(c,strsplit(dir(".",inState),".csv"))

inDates <- paste0(sort(rep(2003:2018,4)),"_",dig2(9:12))

sr_files <- c(paste0("MCD14ML_SR_Terra_",inState,"_",inDates),
              paste0("MCD14ML_SR_Aqua_",inState,"_",inDates))
sr_area_files <- c(paste0("MCD14ML_SR_Area_Terra_",inState,"_",inDates),
                  paste0("MCD14ML_SR_Area_Aqua_",inState,"_",inDates))
allfiles <- c(sr_files,sr_area_files)

print(allfiles[!(allfiles %in% infiles)])

blank_df <- matrix(NA,1,13)
colnames(blank_df) <- c("YYYYMMDD","HHMM","conf","sat","FRP",
              "sur_refl_b01","sur_refl_b02","sur_refl_b03","sur_refl_b04","sur_refl_b05",
              "sur_refl_b06","sur_refl_b07",".geo")
  
for (i in seq_along(sr_files)) {
  inFileName <- paste0(sr_files[i],".csv")
  dim_file <- dim(try(read.csv(inFileName)))[1]

  if (is.null(dim_file)) {write.table(blank_df,inFileName,sep=",",row.names=F)}
}
