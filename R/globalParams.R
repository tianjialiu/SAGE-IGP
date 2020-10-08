# --------------------------------------------
# Global parameters for SAGE-IGP
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# --------------------------------------------
library("raster"); library("rgdal"); library("rgeos"); library("fields");
library("plyr"); library("pracma"); library("leaps"); library("usdm"); 
library("zoo"); library("drc"); library("questionr"); library("SDMTools")
library("ggplot2"); library("tidyverse"); library("ncdf4")

dirHome <- "/Users/TLiu/Google Drive/India/IndiaSurvey"
figHome <- file.path(dirHome,"figures")
tabHome <- file.path(dirHome,"tables")
geeHome <- file.path(dirHome,"gee")
rasHome <- file.path(dirHome,"rasters")

alpha <- 0.41 # FRE to DM, kg/MJ
t_fire <- 60*60*0.5 # duration of fire, sec
pBurnAdj <- 1.92 # Lasko and Vadrevu (2018)
xBeta <- c(0.1,0.5,0.9) # fire season breakpoints
efs_andreae <- data.frame(OC=4.9,BC=0.42,CO=76,CO2=1430,DM=1000)
regionExtent <- extent(c(72,89,23,33))
  
setwd(dirHome)
getSurvey <- function(surveyYr) {
  return(read.csv(paste0("tables/IGP_survey_",surveyYr,".csv")))
}

colorPalette <- c("indianred","skyblue","#fdbf6f",1)
colorPalette2 <- brewer.pal(4,"Pastel1")

getBrkPtdate <- function(x,beta,smooth=T) {
  if (sum(x,na.rm=T) > 0) {
    if (smooth == T) {npartial_sums <- get_nls(x)} else {npartial_sums <- cumsum(x)/sum(x)}
    npartial_sums_diff <- npartial_sums-beta
    beta_date <- which(npartial_sums_diff %in% npartial_sums_diff[npartial_sums_diff>=0])[1]
  } else {beta_date <- NA}
  return(beta_date)
}

getWtdPeak <- function(y,peakVal=F) {
  x <- seq_along(y)
  
  peakDate <- round(weighted.mean(x,y))
  peakMag <- y[which(x==peakDate)]
  if (peakVal==T) {return(c(peakDate,peakMag))} else {return(peakDate)}
}

getGaussPeak <- function(y,peak=F,usable_frac=NA) {
  x <- seq_along(y)
  peakDate <- getWtdPeak(y)
  total_y <- sum(y)
  y <- y/total_y
  
  if (length(usable_frac)==1) {usable_frac <- rep(1,length(x))}
  
  gauss_fit <- function(x,par) {
    m <- par[1]
    sd <- par[2]
    k <- par[3]
    rhat <- k * exp(-0.5 * ((x - m)/sd)^2)
    return(rhat)
  }
  
  min_fn <- function(par) { 
    rhat <- gauss_fit(x,par)
    lsq <- (y - rhat)^2
    
    return(weighted.mean(lsq,usable_frac)*length(lsq))
  }
  
  gauss_peak <- optim(c(peakDate, 7, 1), min_fn, method="BFGS", control=list(reltol=1e-9))
  gauss_y <- gauss_fit(x,gauss_peak$par)
  
  if (peak==T) {return(round(gauss_peak$par[1]))} else {return(gauss_y)}
  
}

getWtdBounds <- function(y,smooth=T) {
  x <- seq_along(y)
  totY <- sum(y)
  
  if (smooth == T) {npartial_sums <- get_nls(y)} else {npartial_sums <- cumsum(y)/sum(y)}
  y <- c(0,diff(npartial_sums))*totY
  
  x_distMid <- wtd.mean(x,y)
  x_distMin <- wtd.mean(x,y)-sqrt(wtd.var(x,y))*1.25
  x_distMax <- wtd.mean(x,y)+sqrt(wtd.var(x,y))*1.25
  
  return(round(c(x_distMin,x_distMid,x_distMax)))
}

getHistogram <- function(x) {
  histFreq <- do.call(c,lapply(lapply(strsplit(strsplit(as.character(x),"\\[")[[1]][-c(1:2)], ", "),
                                      function(x) {strsplit(x[2],"\\]")}),function(x) {x[[1]][1]}))
  return(as.numeric(histFreq))
}

get_nls <- function(y) {
  y <- cumsum(y)/sum(y)
  x <- seq_along(y)
  nls_fit <- nls(y ~ 1/(1 + exp(a + b*x)), start=list(a=0, b=0), control=nls.control(maxiter = 100))
  npartial_sums <- predict(nls_fit)/max(predict(nls_fit))
  return(npartial_sums)
}

get_d_nls <- function(y) {
  npartial_sums <- get_nls(y)
  d_nls <- c(0,npartial_sums[-1]-npartial_sums[-length(npartial_sums)])
  return(d_nls)
}

getGaussShape <- function(x,betaPts) {
  betaPts <- as.numeric(betaPts)
  gaussTS <- exp(-0.5 * ((x - betaPts[2])/((betaPts[3]-betaPts[1])/2.5))^2)
  return(gaussTS)
}

perChange <- function(x,y=NULL,doRound=T) {
  x <- as.numeric(x)
  
  if (is.null(y) == T) {
    perChange <- (x[length(x)]-x[1])/x[1]*100
  } else {
    y <- as.numeric(y)
    perChange <- (y-x)/x*100
  }
  
  if (doRound == T) {perChange <- round(perChange)}
  
  return(perChange)
}

toDbl <- function(x) {
  return(sprintf("%02d",x))
}

blankDates <- function(sMonth,eMonth,inYears,NAcols=F,cutMonths=F) {
  sYear <- inYears[1]
  eYear <- inYears[length(inYears)]
  
  if (inYears[length(inYears)] %% 4 == 0) {
    nDays <- c(31,29,31,30,31,30,31,31,30,31,30,31)
  } else {nDays <- c(31,28,31,30,31,30,31,31,30,31,30,31)}
  
  sDate <- paste0(sYear,"-",sprintf("%02d",sMonth),"-01")
  eDate <- paste0(eYear,"-",sprintf("%02d",eMonth),"-",sprintf("%02d",nDays[eMonth]))
  datesAll <- seq(as.Date(sDate),as.Date(eDate),"day")
    
  YEARS <- as.numeric(format(datesAll,"%Y"))
  MONTHS <- as.numeric(format(datesAll,"%m"))
  DAYS <- as.numeric(format(datesAll,"%d"))
  JDAYS <- as.numeric(format(datesAll,"%j"))
  
  datesTable <- data.frame(cbind(YEARS,MONTHS,DAYS,JDAYS))
  colnames(datesTable) <- c("Year","Month","Day","Julian")
  if (cutMonths == T) {datesTable <- datesTable[datesTable$Month %in% sMonth:eMonth,]}
  
  if (is.character(NAcols) == T | is.numeric(NAcols) == T ) {
    datesTable <- cbind(datesTable,matrix(NA,dim(datesTable)[1],length(NAcols)))
    colnames(datesTable) <- c(colnames(datesTable)[1:4],NAcols)
  }
  
  return(datesTable)
}

datefromYMD <- function(blankDatesMat) {
  ymdDate <- as.Date(paste0(blankDatesMat$Year,"-",
                     dig2(blankDatesMat$Month),"-",
                     dig2(blankDatesMat$Day)))
  return(ymdDate)
}

