# -------------------------------------------------------
# Grid state-level dry matter at 0.25°x0.25° resolution
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# -------------------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')

xStates <- c("Punjab","Haryana","Uttar_Pradesh","Bihar","Rajasthan")

saveDMras <- function(inState,doVIIRS) {
  if (doVIIRS == T) {
    xYears <- 2012:2018
  } else {
    xYears <- 2003:2018
  }
  
  setwd(tabHome)
  adjDM <- read.csv(paste0("MODISadjFRP_DM/MODISadjFRP_DM_",inState,".csv"))
  FRPyrAll <- read.csv(paste0("MODIS_FRP_Grid/MODIS_FRP_Grid_",inState,".csv"))
  FRPyrAll_Aqua <- read.csv(paste0("MODIS_FRP_Grid/MODIS_Aqua_FRP_Grid_",inState,".csv"))
  
  if (doVIIRS == T) {
    FRPyrAll_VIIRS <- read.csv(paste0("VIIRS_FRPboost_Grid/VNP14IMGML_FRPboost_Grid_",inState,".csv"))
  }
  
  gridIDs <- as.numeric(substr(colnames(FRPyrAll)[-c(1:4)],2,7))
  gridRas <- raster(res=c(0.25,0.25))
  inMonths <- c(9:12)
  nDays <- 122
  
  vec2Ras <- function(x) {
    valVec <- rep(NA,length(gridRas))
    valVec[gridIDs] <- x
    valRas <- setValues(gridRas,valVec)
    valRas <- crop(valRas,regionExtent)
    return(valRas)
  }
  
  timestamp(prefix=paste(inState,": ##------ "))
  
  for (iYear in seq_along(xYears)) {
    
    setwd(tabHome)
    
    # adjusted DM
    adjDMyr <- adjDM[,paste0("Y",xYears[iYear])] 
    fireMat <- read.csv(paste0("MODISadjFRP/MODISadjFRP_",inState,"_",xYears[iYear],".csv"))
    adjFactor <- fireMat$adjFactor
    
    # MODIS FRP, daily, by grid cell
    FRPyr <- FRPyrAll[FRPyrAll$Year == xYears[iYear] & FRPyrAll$Month %in% inMonths, -c(1:4)]
    FRPyr <- sweep(FRPyr,1,adjFactor,"*")
    
    if (xYears[iYear] >= 2012 & doVIIRS == T) {
      FRPyr_VIIRS <- FRPyrAll_VIIRS[FRPyrAll_VIIRS$Year == xYears[iYear], -c(1:4)]
      FRPyr_VIIRS <- sweep(FRPyr_VIIRS,1,adjFactor,"*")
      
      FRPyr_Aqua <- FRPyrAll_Aqua[FRPyrAll_Aqua$Year == xYears[iYear] & FRPyrAll_Aqua$Month %in% inMonths, -c(1:4)]
      FRPyr_Aqua <- sweep(FRPyr_Aqua,1,adjFactor,"*")
      
      FRPyr_Terra <- FRPyr - FRPyr_Aqua
      VIIRS_Aqua_ratio <- FRPyr_VIIRS/FRPyr_Aqua
      VIIRS_Aqua_ratio[is.na(VIIRS_Aqua_ratio)] <- 0
      VIIRS_Aqua_ratio[VIIRS_Aqua_ratio==Inf] <- 1
      
      FRPyr <- FRPyr + FRPyr_VIIRS + FRPyr_Terra*VIIRS_Aqua_ratio
    }
    
    FRPyrTot <- colSums(FRPyr) # grid cell seasonal total
    FRPyrTotAll <- sum(FRPyrTot) # state seasonal total
    
    # Step 1: Disaggregate total DM spatially
    adjDMgrid <- FRPyrTot/FRPyrTotAll*sum(adjDMyr)
    
    # Step 2: Adjust daily timeseries by Gaussian distribution
    grid_xBeta <- apply(FRPyr,2,function(x) {
      tapply(xBeta,xBeta,function(iBeta) {getBrkPtdate(x,iBeta,smooth=F)})
    })
    
    soloDays <- which((grid_xBeta[3,]-grid_xBeta[1,])==0)
    grid_xBeta[1,soloDays] <- grid_xBeta[1,soloDays] - 1
    grid_xBeta[3,soloDays] <- grid_xBeta[3,soloDays] + 1
    grid_gauss <- apply(grid_xBeta,2,function(betaPts) getGaussShape(1:nDays,betaPts))
    
    adjDMgridTS <- matrix(0,nDays,length(adjDMgrid))
    for (id in seq(adjDMgrid)) {
      adjDMgridTS[,id] <- grid_gauss[,id] * adjDMgrid[id]
    }
    
    # Step 3: Conserve total grid DM to original state DM and FRP spatial distribution
    scaleGrid <- adjDMgrid/colSums(adjDMgridTS,na.rm=T)
    adjDMgridTS <- sweep(adjDMgridTS,2,scaleGrid,"*")
    
    tempGrid <- adjDMyr/rowSums(adjDMgridTS,na.rm=T)
    adjDMgridTS <- sweep(adjDMgridTS,1,tempGrid,"*")
    
    # to convergence
    xIter <- 0
    while(mean(abs(1 - scaleGrid[scaleGrid < Inf & !is.na(scaleGrid)])) > 0.05 & xIter < 1e3) {
      xIter <- xIter + 1
      scaleGrid <- adjDMgrid/colSums(adjDMgridTS,na.rm=T)
      adjDMgridTS <- sweep(adjDMgridTS,2,scaleGrid,"*")
      
      tempGrid <- adjDMyr/rowSums(adjDMgridTS,na.rm=T)
      adjDMgridTS <- sweep(adjDMgridTS,1,tempGrid,"*")
    }
    
    adjDMgridTS[is.na(adjDMgridTS)] <- 0
    
    # write to raster
    DMyr_ras <- vec2Ras(colSums(adjDMgridTS,na.rm=T))
    adjDMgridTSras <- stack(apply(adjDMgridTS,1,vec2Ras))
    
    setwd(file.path(rasHome,"adjFRP_DMstate"))
    if (doVIIRS == F) {
      writeRaster(DMyr_ras,paste0("adjFRP_DM_",xYears[iYear],"_",inState,".tif"),format="GTiff",overwrite=T)
      writeRaster(adjDMgridTSras,paste0("adjFRP_DM_daily_",xYears[iYear],"_",inState,".tif"),format="GTiff",overwrite=T)
    } else {
      writeRaster(DMyr_ras,paste0("adjFRP_DM_VIIRS_",xYears[iYear],"_",inState,".tif"),format="GTiff",overwrite=T)
      writeRaster(adjDMgridTSras,paste0("adjFRP_DM_dailyVIIRS_",xYears[iYear],"_",inState,".tif"),format="GTiff",overwrite=T)
    }

    timestamp(prefix=paste("Year",xYears[iYear],": ##------ "))
  }
}

for (doVIIRS in c(F,T)) {
  for (inState in xStates) {
    saveDMras(inState,doVIIRS)
  }
}
