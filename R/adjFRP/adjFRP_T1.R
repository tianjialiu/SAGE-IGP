# ---------------------------------------------------
# Adjusted FRP for Tier 1 states - Punjab, Haryana
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# ---------------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')

xMonths <- 9:12
inState <- "Punjab"
xYears <- 2003:2018

nWindow <- 7 # size of window to gap fill for clouds/haze
neighborSensAdj <- 0.75 # sensitivity of cloud-gap fill procedure

timestamp(prefix=paste(inState,": ##------ "))

# -------
# Survey
# -------
xStates <- c("Haryana","Punjab","Uttar_Pradesh","Bihar")
iState <- which(xStates==inState)

survey <- getSurvey(2017)
surveySt <- survey[survey$State==iState,]
surveySt_area <- surveySt$Operated_Landholding_ha

# method of burning: partial or complete burn
surveySt_methodBurn <- surveySt$Burn_entire_field_stack_center_burn
methodBurn <- as.numeric(tapply(surveySt_area,surveySt_methodBurn,sum))
methodBurn_per <- methodBurn/sum(methodBurn)
partialBurnRatio <- methodBurn_per[2]
wholeBurnRatio <- methodBurn_per[1]

nChoices <- 3
surveySt_diurnal <- surveySt[,c(paste0("Time_of_burning_",1:nChoices),"Operated_Landholding_ha")]
row_idx <- rowSums(!is.na(surveySt_diurnal[,1:nChoices]))

# diurnal cycle
diurnal <- rep(0,4)
for (iTime in 1:nChoices) {
  surveySt_diurnalSub <- surveySt_diurnal[row_idx==iTime,]
  
  for (iChoice in 1:iTime) {
    surveySt_diurnalTab <- tapply(surveySt_diurnalSub[,nChoices+1],surveySt_diurnalSub[,iTime],sum)
    diurnal_idx <- as.numeric(names(surveySt_diurnalTab))
    diurnal[diurnal_idx] <- as.numeric(surveySt_diurnalTab)*(1/iTime) + diurnal[diurnal_idx]
  }
}

diurnalPer <- diurnal/sum(diurnal)
diurnalSat <- diurnalPer[2]

# ----------
# Satellite
# ----------
# cloud-gap fill detection rate cutoffs
cutoff_st <- 0.145
cutoff_end <-  0.375

# viirs-aqua ratio
viirsAquaRatio <- c(0.88,0.85,1.64,1.43)[iState]

# clouds/haze rate of detection approximation
getObsFrac <- function(SRred) {
  
  if (SRred >= cutoff_end) {obsFrac <- 0}
  if (SRred <= cutoff_st) {obsFrac <- 1}
  if (SRred < cutoff_end & SRred > cutoff_st) {
    obsFrac <- 1-(SRred-cutoff_st)/(cutoff_end-cutoff_st)
  }
  return(obsFrac)
}

histSteps <- seq(0.005,0.995,0.01)
xObsFrac <- tapply(histSteps,seq_along(histSteps),function(x) {getObsFrac(x)})

windowBlock1 <- 1:nWindow
windowBlock2 <- (nWindow+2):(nWindow*2+1)
midWindow <- nWindow + 1

for (iYear in xYears) {
  
  fireMat <- blankDates(xMonths[1],xMonths[length(xMonths)],iYear,
                        NAcols=c("aquaObsFrac","terraObsFrac","satObsFrac","terraFRP","aquaFRP","viirsFRP","satFRP",
                                 "viirsFRPboost","cloudBoost","diurnalBoost","pBurnBoost","adjFRP"))
  nObs <- length(fireMat$adjFRP)
  
  # Step 1: Add VIIRS small fire boost
  
  for (inSat in c("Aqua","Terra")) {
    if (inSat == "Terra") {inMODISsat <- "MOD14A1"}
    if (inSat == "Aqua") {inMODISsat <- "MYD14A1"}
    # read MODIS FRP
    modisFRPgrid <- read.csv(paste0("gee/MxD14A1/MxD14A1_FRP_",inState,".csv"))
    
    for (iMonth in xMonths) {
      monthday <- which(fireMat$Month==iMonth)
      
      if (iYear >= 2012) {
        # read VIIRS FRP
        viirs_frp <- read.csv(paste0("gee/VNP14IMGML_FRPboost/VNP14IMGML_FRP_",inState,"_",
                                     iYear,"_",sprintf("%02d",iMonth),".csv"))
        fireMat$viirsFRP[monthday] <- viirs_frp$FRP
        fireMat$viirsFRPboost[monthday] <- viirs_frp$FRPboost
      }
      
      # read cloud/haze gap fill tables
      fire_sr_area_name <- paste0("gee/MCD14ML_SR/MCD14ML_SR_Area_",inSat,"_",inState,
                                  "_",iYear,"_",sprintf("%02d",iMonth),".csv")
      fire_sr_area <- read.csv(fire_sr_area_name)
      
      for (iDay in seq_along(monthday)) {
        frpDay <- modisFRPgrid[which(as.Date(as.character(modisFRPgrid$YYYYMMDD),"%Y%m%d")==
                                       as.Date(paste0(iYear,"-",sprintf("%02d",iMonth),"-",sprintf("%02d",iDay)))),
                               inMODISsat]
        dayIdx <- which(fireMat$Month==iMonth & fireMat$Day==iDay)
        if (length(frpDay) == 0) {frpDay <- 0}
        
        if (!nchar(as.character(fire_sr_area$histogram[iDay])) %in% c(0,NA)) {
          hist_wsteps_day <- getHistogram(fire_sr_area$histogram[iDay])
          obsFracDay <- weighted.mean(xObsFrac,hist_wsteps_day,na.rm=T)
        } else {obsFracDay <- 1}
        
        fireMat[dayIdx,paste0(tolower(inSat),"FRP")] <- frpDay
        fireMat[dayIdx,paste0(tolower(inSat),"ObsFrac")] <- obsFracDay
      }
    }
    
    # add boost from overall lack of detection from clouds/haze
    obsFrac <- fireMat[,paste0(tolower(inSat),"ObsFrac")]
    if (inSat == "Aqua") {
      # for years prior to 2012, approximate that VIIRS will add 0.85 of MODIS FRP
      if (iYear >= 2012) {fireMat$satFRP <- fireMat$aquaFRP + fireMat$viirsFRPboost}
      if (iYear < 2012) {fireMat$satFRP <- fireMat$aquaFRP + viirsAquaRatio*fireMat$aquaFRP}
      
      # Use detection rate for initial boost from cloud/haze gaps
      initCloudBoost <- fireMat$satFRP/obsFrac - fireMat$satFRP
      initCloudBoost[is.na(initCloudBoost)] <- 0
      satBoost <- fireMat$satFRP + initCloudBoost
    }
    
    if (inSat == "Terra") {
      terraFRPBoost <- fireMat$terraFRP + viirsAquaRatio*fireMat$terraFRP
      fireMat$satFRP <- fireMat$satFRP + terraFRPBoost
      
      initCloudBoost <- terraFRPBoost/obsFrac - terraFRPBoost
      initCloudBoost[is.na(initCloudBoost)] <- 0
      satBoost <- satBoost + terraFRPBoost + initCloudBoost
    }
  }
  
  # Step 2: Gap fill for clouds/haze
  
  # weighted Aqua/Terra detection rate
  wObsFrac <- (fireMat$terraObsFrac*fireMat$terraFRP + fireMat$aquaObsFrac*fireMat$aquaFRP)/
    (fireMat$terraFRP+fireMat$aquaFRP)
  wObsFrac[is.na(wObsFrac)] <- apply(fireMat[,c("terraObsFrac","aquaObsFrac")],1,
                                     function(x) {weighted.mean(x,c(sum(fireMat$terraFRP),
                                                                    sum(fireMat$aquaFRP)))})[is.na(wObsFrac)]
  fireMat$satObsFrac <- wObsFrac
  
  satFRPmax <- sum(satBoost)
  
  # cumulative FRP, normalized by maximum
  cFRP <- cumsum(satBoost/satFRPmax)
  x1 <- max(which(round(cFRP,2)==0))
  x2 <- min(which(round(cFRP,2)==1))
  
  cFRP_iIter <- rep(1,length(cFRP))
  
  # clouds/gap fill by iteration until convergence
  while (round(sum(cFRP_iIter),1) > 0) {
    cFRPprior <- cFRP
    dFRP <- c(0,cFRP[-1]-cFRP[-nObs])
    frpPeak <- weighted.mean(seq_along(cFRP),dFRP,na.rm=T) # midpoint of fire season
    frpNPS <- get_d_nls(dFRP)
    
    for (iFRP in seq_along(cFRP)[-c(1:nWindow,(nObs-nWindow):nObs)]) {
      # only adjust FRP for days within the fire season
      if (iFRP > x1 & iFRP < x2) {
        neighbor_cFRP <- cFRP[(iFRP-nWindow):(iFRP+nWindow)]
        neighbor_dFRP <- dFRP[(iFRP-nWindow):(iFRP+nWindow)]
        neighbor_nps <- frpNPS[(iFRP-nWindow):(iFRP+nWindow)]
        
        # slopes of cumulative distribution derived from linear regression
        slope1 <- coef(lm(neighbor_cFRP[c(windowBlock1,midWindow)]~
                            seq_along(neighbor_cFRP)[c(windowBlock1,midWindow)]))[2]
        slope2 <- coef(lm(neighbor_cFRP[c(midWindow,windowBlock2)]~
                            seq_along(neighbor_cFRP)[c(midWindow,windowBlock2)]))[2]
        slope <- coef(lm(neighbor_cFRP~seq_along(neighbor_cFRP)))[2]
        
        neighbor_off1 <- max(neighbor_dFRP[windowBlock1])
        neighbor_off2 <- max(neighbor_dFRP[windowBlock2])
        neighbor_off <- max(neighbor_dFRP)
        
        neighborSens <- neighbor_nps[midWindow]/max(neighbor_nps)
        
        # for days at the edges of the distribution, use slopes that were
        # calculated excluding days outside the fire season
        inSlope <- slope
        if (iFRP < x1+nWindow) {inSlope <- slope2}
        if (iFRP > x2-nWindow) {inSlope <- slope1}
        
        # stop adjustments if criteria are met
        if (neighbor_dFRP[midWindow] == 0 | 
            iFRP < frpPeak & neighbor_off1 > 0 & neighbor_dFRP[midWindow]/neighbor_off1 > neighborSens*neighborSensAdj |
            iFRP > frpPeak & neighbor_off2 > 0 & neighbor_dFRP[midWindow]/neighbor_off > neighborSens*neighborSensAdj) {
          inSlope <- 0
        }
        
        # adjust cFRP at a rate proportional to (1 - detection rate) so that
        # cloudy/hazy days have priority
        cFRP[iFRP:nObs] <- cFRP[iFRP:nObs] + (inSlope * (1-wObsFrac[iFRP]))
      }
    }
    cFRP_iIter <- abs(cFRP - cFRPprior)
  }
  
  satBoostCloud <- c(0,cFRP[-1]-cFRP[-length(cFRP)])*satFRPmax
  adjFactor <- (satBoost/satBoostCloud)[which(satBoost==max(satBoost))]
  fireMat$adjFRP <- satBoostCloud*adjFactor
  
  # gap fill zero values outside the main fire season by taking the mean
  # of FRP before and after the given day
  fireAbsSeas <- which(fireMat$adjFRP>0)
  zeroHoles <- which(fireMat$adjFRP==0)
  zeroHoles <- zeroHoles[which(zeroHoles > min(fireAbsSeas) & zeroHoles < max(fireAbsSeas))]
  
  while(length(zeroHoles) > 0) {
    fillGap <- cbind(fireMat$adjFRP[zeroHoles-1],fireMat$adjFRP[zeroHoles+1])
    fireMat$adjFRP[zeroHoles] <- rowMeans(fillGap,na.rm=T)
    
    zeroHoles <- which(fireMat$adjFRP==0)
    zeroHoles <- zeroHoles[which(zeroHoles > min(fireAbsSeas) & zeroHoles < max(fireAbsSeas))]
  }
  
  # cloud/haze gap-fill boost
  fireMat$cloudBoost <- fireMat$adjFRP-fireMat$satFRP
  fireMat$cloudBoost[fireMat$cloudBoost<0] <- 0
  
  # time series fraction, plus more adj. ts nudged to original MODIS+VIIRS FRP ts
  tsFrac <- fireMat$adjFRP/sum(fireMat$adjFRP)
  
  cloud_consAdj <- tapply(1:nObs,1:nObs,function(x) {
    weighted.mean(c(fireMat$satFRP[x],fireMat$adjFRP[x]),
                  c(wObsFrac[x],1-wObsFrac[x]))})
  
  tsFrac_adj <- cloud_consAdj/sum(cloud_consAdj)
  
  fireMat$cloudBoost <- tsFrac_adj/tsFrac * fireMat$cloudBoost
  fireMat$adjFRP <- fireMat$satFRP + fireMat$cloudBoost
  
  # Step 3: partial burning
  fireMat$pBurnBoost <- fireMat$adjFRP*(0.75*partialBurnRatio+wholeBurnRatio)/wholeBurnRatio - fireMat$adjFRP
  fireMat$adjFRP <- fireMat$adjFRP + fireMat$pBurnBoost
  
  # Step 4: diurnal cycle
  terraFrac <- sum(fireMat$terraFRP)/sum(fireMat$terraFRP+fireMat$aquaFRP)
  aquaFrac <- 1-terraFrac
  fireMat$diurnalBoost <- (fireMat$adjFRP*2.67)/diurnalSat - fireMat$adjFRP
  fireMat$adjFRP <- fireMat$adjFRP + fireMat$diurnalBoost
  
  fireMat[is.na(fireMat)] <- 0
  
  # remove anomalous spikes in FRP
  rollWinHalf <- 2
  rollIdx <- 1:nObs
  
  adjFRP_noSpike <- fireMat$adjFRP
  rollStats <- data.frame(do.call(rbind,tapply(1:nObs,1:nObs,function(x) {
    idx <- (x-rollWinHalf):(x+rollWinHalf)
    idx <- idx[which(idx>0 & idx <= nObs)]
    idx <- idx[-which(idx==x)]
    return(c(mean(adjFRP_noSpike[idx]),sd(adjFRP_noSpike[idx]),
             max(adjFRP_noSpike[idx])))
  })))
  
  colnames(rollStats) <- c("mean","sd","max")
  
  spikeValid <- adjFRP_noSpike/rollStats$max
  spikeRmv <- quantile(adjFRP_noSpike,0.25)
  
  adjFRP_noSpike[which(spikeValid>3 & adjFRP_noSpike>spikeRmv)] <- 
    (rollStats$mean + rollStats$sd)[which(spikeValid>3 & adjFRP_noSpike>spikeRmv)]
  
  fireMat$adjFactor <- adjFRP_noSpike/fireMat$adjFRP
  fireMat$adjFRP <- adjFRP_noSpike
  
  fireMat[is.na(fireMat)] <- 0
  
  write.table(fireMat,paste0("tables/MODISadjFRP/MODISadjFRP_",inState,"_",iYear,".csv"),sep=",",row.names=F)
  timestamp(prefix=paste("Year",iYear,": ##------ "))
}
