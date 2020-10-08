# --------------------------------------------
# Combine gridded state-level dry matter
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# --------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')

doVIIRS <- F
inSpecies <- "DM"
timeScale <- "daily"

if (doVIIRS == T) {
  xYears <- 2012:2018
} else {
  xYears <- 2003:2018
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
    
    yrEmiSum <- sum(stack(haryana,punjab,up,bihar,rajasthan),na.rm=T)
    
    setwd(file.path(rasHome,"adjFRP_Inv"))
    writeRaster(yrEmiSum,paste0(inSpecies,"/adjFRP_DM",file_suffix,"_",inYear,".tif"),
                format="GTiff",overwrite=T)
    
    timestamp(prefix=paste("Year",inYear,": ##------ "))
  }
}

if (timeScale == "daily") {
  for (inYear in xYears) {
    setwd(file.path(rasHome,"adjFRP_DMstate"))
    
    file_suffix <- ""
    if (doVIIRS == T) {file_suffix <- "VIIRS"}
    
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
      dayEmi <- sum(stack(haryana[[iDay]],punjab[[iDay]],
                          up[[iDay]],bihar[[iDay]],rajasthan[[iDay]]),na.rm=T)
      return(dayEmi)
    })
    
    yrEmi <- stack(do.call(list,yrEmi))
    
    setwd(file.path(rasHome,"adjFRP_Inv"))
    writeRaster(yrEmi,paste0(inSpecies,"/adjFRP_",inSpecies,"_daily",file_suffix,"_",inYear,".tif"),
                format="GTiff",overwrite=T)
    
    timestamp(prefix=paste("Year",inYear,": ##------ "))
  }
}
