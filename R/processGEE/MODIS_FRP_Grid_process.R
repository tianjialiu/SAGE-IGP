# --------------------------------------------
# Pre-process MxD14A1 FRP files from EE
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# --------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')

xYears <- 2003:2018
inState <- "Punjab"

setwd(file.path(geeHome,"MxD14A1_Grid"))
mxd14a1_files <- dir(".",paste0("MxD14A1_FRP_",inState))
mxd14a1 <- read.csv(dir(".",mxd14a1_files[1]))
mxd14a1_dates <- as.Date(as.character(mxd14a1$YYYYMMDD),"%Y%m%d")

gridIDs <- sort(as.numeric(substr(mxd14a1_files,19+nchar(inState),24+nchar(inState))))

FRPyrAll <- blankDates(1,12,xYears,paste0("G",gridIDs))
FRPyrAll_Aqua <- blankDates(1,12,xYears,paste0("G",gridIDs))
inDates <- as.Date(paste0(FRPyrAll$Year,"-",toDbl(FRPyrAll$Month),"-",toDbl(FRPyrAll$Day)))

for (id in seq_along(gridIDs)) {
  setwd(file.path(geeHome,"MxD14A1_Grid"))
  mxd14a1_id <- read.csv(paste0("MxD14A1_FRP_",inState,"_Grid_",gridIDs[id],".csv"))
  FRPyrAll[inDates %in% mxd14a1_dates, paste0("G",gridIDs[id])] <- mxd14a1_id$MxD14A1
  FRPyrAll_Aqua[inDates %in% mxd14a1_dates, paste0("G",gridIDs[id])] <- mxd14a1_id$MYD14A1
}

FRPyrAll[is.na(FRPyrAll)] <- 0
FRPyrAll_Aqua[is.na(FRPyrAll_Aqua)] <- 0

setwd(file.path(tabHome,"MODIS_FRP_Grid"))
write.table(FRPyrAll,paste0("MODIS_FRP_Grid_",inState,".csv"),sep=",",row.names=F)
write.table(FRPyrAll_Aqua,paste0("MODIS_Aqua_FRP_Grid_",inState,".csv"),sep=",",row.names=F)
