# --------------------------------------------
# Pre-process VIIRS FRP boost files from EE
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# --------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')

xYears <- 2012:2018
inState <- "Punjab"

setwd(file.path(geeHome,"VNP14IMGML_FRPboost_Grid"))
vnp14_files <- dir(".",inState)
vnp14 <- read.csv(dir(".",vnp14_files[1]))
vnp14_dates <- as.Date(as.character(vnp14$YYYYMMDD),"%Y%m%d")

gridIDs <- sort(as.numeric(substr(vnp14_files,22+nchar(inState),27+nchar(inState))))

FRPyrAll <- blankDates(1,12,xYears,paste0("G",gridIDs))
FRPyrAll <- FRPyrAll[FRPyrAll$Month %in% c(9:12),]
inDates <- as.Date(paste0(FRPyrAll$Year,"-",toDbl(FRPyrAll$Month),"-",toDbl(FRPyrAll$Day)))

for (id in seq_along(gridIDs)) {
  vnp14_id <- read.csv(paste0("VNP14IMGML_FRP_",inState,"_Grid_",gridIDs[id],".csv"))
  FRPyrAll[inDates %in% vnp14_dates, paste0("G",gridIDs[id])] <- vnp14_id$FRPboost
}

setwd(file.path(tabHome,"VIIRS_FRPboost_Grid"))
write.table(FRPyrAll,paste0("VNP14IMGML_FRPboost_Grid_",inState,".csv"),sep=",",row.names=F)
