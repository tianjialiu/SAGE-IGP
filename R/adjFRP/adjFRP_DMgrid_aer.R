# --------------------------------------------------
# Combine gridded state-level dry matter, aerosols
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# --------------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')

doVIIRS <- F
inSpecies <- "DMaer"
timeScale <- "daily"

if (doVIIRS == T) {
  xYears <- 2012:2018
} else {
  xYears <- 2003:2018
}

xStates <- c("Punjab","Haryana","Uttar_Pradesh","Bihar")

pBurnFracSt <- rep(NA,length(xStates))
for (iState in seq_along(xStates)) {
  fireMat <- read.csv(paste0("tables/MODISadjFRP/MODISadjFRP_",xStates[iState],"_",xYears[1],".csv"))
  pBurnFracSt[iState] <- sum(fireMat$pBurnBoost)/sum(fireMat[,c("satFRP","cloudBoost","pBurnBoost")])
}

if (timeScale == "seasonal") {
  for (inYear in xYears) {
    
    setwd(file.path(rasHome,"adjFRP_DMstate"))
    
    file_suffix <- ""
    if (doVIIRS == T) {file_suffix <- "_VIIRS"}
    
    readRaster <- function(stateName) {
      inRas <- raster(paste0("adjFRP_DM",file_suffix,"_",inYear,"_",stateName,".tif"))
      return(inRas)
    }
    
    haryana <- readRaster("Haryana")
    punjab <- readRaster("Punjab")
    up <- readRaster("Uttar_Pradesh")
    bihar <- readRaster("Bihar")
    rajasthan <- readRaster("Rajasthan")
    
    yrEmiSum <- sum(stack(punjab*(pBurnFracSt[1]*pBurnAdj + 1-pBurnFracSt[1]),
                          haryana*(pBurnFracSt[2]*pBurnAdj + 1-pBurnFracSt[2]),
                          up*(pBurnFracSt[3]*pBurnAdj + 1-pBurnFracSt[3]),
                          bihar*(pBurnFracSt[4]*pBurnAdj + 1-pBurnFracSt[4]),
                          rajasthan*(pBurnFracSt[2]*pBurnAdj + 1-pBurnFracSt[2])),na.rm=T)
    
    setwd(file.path(rasHome,"adjFRP_Inv"))
    writeRaster(yrEmiSum,paste0(inSpecies,"/adjFRP_",inSpecies,file_suffix,"_",inYear,".tif"),
                format="GTiff",overwrite=T)
    
    timestamp(prefix=paste("Year",inYear,": ##------ "))
  }
}

if (timeScale == "daily") {
  for (inYear in xYears) {
    
    setwd(file.path(rasHome,"adjFRP_DMstate"))
    
    file_suffix <- ""
    if (doVIIRS == T) {file_suffix <- "VIIRS"; file_suffix2 <- paste0("_",file_suffix)}

    readStack <- function(stateName) {
      stStack <- stack(paste0("adjFRP_DM_daily",file_suffix,"_",inYear,"_",stateName,".tif"))
      return(stStack)
    }
    
    haryana <- readStack("Haryana")
    punjab <- readStack("Punjab")
    up <- readStack("Uttar_Pradesh")
    bihar <- readStack("Bihar")
    rajasthan <- readStack("Rajasthan")
    
    yrEmi <- tapply(1:122,1:122,function(iDay) {
      dayEmi <- sum(stack(punjab[[iDay]]*(pBurnFracSt[1]*pBurnAdj + 1-pBurnFracSt[1]),
                          haryana[[iDay]]*(pBurnFracSt[2]*pBurnAdj + 1-pBurnFracSt[2]),
                          up[[iDay]]*(pBurnFracSt[3]*pBurnAdj + 1-pBurnFracSt[3]),
                          bihar[[iDay]]*(pBurnFracSt[4]*pBurnAdj + 1-pBurnFracSt[4]),
                          rajasthan[[iDay]]*(pBurnFracSt[2]*pBurnAdj + 1-pBurnFracSt[2])),na.rm=T)
      return(dayEmi)
    })
    
    yrEmi <- stack(do.call(list,yrEmi))
    
    setwd(file.path(rasHome,"adjFRP_Inv"))
    writeRaster(yrEmi,paste0(inSpecies,"/adjFRP_",inSpecies,"_daily",file_suffix,"_",inYear,".tif"),
                format="GTiff",overwrite=T)
    
    timestamp(prefix=paste("Year",inYear,": ##------ "))
  }
}
