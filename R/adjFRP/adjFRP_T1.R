# ---------------------------------------------------
# Adjusted FRP for Tier 1 states - Punjab, Haryana
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: October 11, 2020
# ---------------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')

xMonths <- 9:12
inState <- "Punjab"
xYears <- 2003:2018

nWindow <- 7 # size of window to gap fill for clouds/haze
neighborSensAdj <- 0.75 # sensitivity of cloud-gap fill procedure
dupRatio <- 0.97 # scaling for duplicate observations

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
cutoff_end <- 0.375
b_coefs <- c(-1.142,17.713,-69.968)

# viirs-aqua ratio
viirsAquaRatio <- c(0.88,0.85,1.64,1.43)[iState]

# a1-ml ratio (MODIS/Aqua)
a1_mlRatio <- c(1.25,1.34,0.93,0.80)[iState]*(375^2/1000^2)

# clouds/haze rate of detection approximation
getObsFrac <- function(SRred) {
  
  if (SRred >= cutoff_end) {obsFrac <- 0}
  if (SRred <= cutoff_st) {obsFrac <- 1}
  if (SRred < cutoff_end & SRred > cutoff_st) {
    obsFrac <- exp(b_coefs[1]+(SRred*b_coefs[2])+(SRred^2*b_coefs[3]))
    obsFrac[obsFrac>1] <- 1
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
    
    # account for possible duplicate observations
    fireMat$terraFRP <- fireMat$terraFRP*dupRatio
    fireMat$aquaFRP <- fireMat$aquaFRP*dupRatio
    fireMat$viirsFRPboost <- fireMat$viirsFRPboost*dupRatio*(1+a1_mlRatio) # account for the ml/a1 difference
    
    # add boost from overall lack of detection from clouds/haze
    obsFrac <- fireMat[,paste0(tolower(inSat),"ObsFrac")]
    obsFrac[obsFrac<0.01] <- 0.01 # set minimum to 0.01
    
    if (inSat == "Aqua") {
      # for years prior to 2012, approximate that VIIRS will add an avg boost of MODIS FRP
      if (iYear >= 2012) {fireMat$satFRP <- fireMat$aquaFRP + fireMat$viirsFRPboost}
      if (iYear < 2012) {
        viirsBoostAqua <- fireMat$aquaFRP*dupRatio*(viirsAquaRatio+a1_mlRatio)
        fireMat$viirsFRPboost <- viirsBoostAqua
        fireMat$satFRP <- fireMat$aquaFRP + viirsBoostAqua
      }
      
      # Use detection rate for initial boost from cloud/haze gaps
      initCloudBoost <- fireMat$satFRP/obsFrac - fireMat$satFRP
      initCloudBoost[is.na(initCloudBoost)] <- 0
      satBoost <- fireMat$satFRP + initCloudBoost
    }
    
    if (inSat == "Terra") {
      viirsBoostTerra <- fireMat$terraFRP*dupRatio*(viirsAquaRatio+a1_mlRatio)
      fireMat$viirsFRPboost <- fireMat$viirsFRPboost+viirsBoostTerra
      
      terraFRPBoost <- fireMat$terraFRP + viirsBoostTerra
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
  
  satFRPmax <- sum(satBoost)
  
  # cumulative FRP, normalized by maximum
  cFRP <- cumsum(satBoost/satFRPmax)
  dFRP <- c(0,cFRP[-1]-cFRP[-nObs])
  
  x1 <- getBrkPtdate(satBoost,0.025)
  x2 <- getBrkPtdate(satBoost,0.975)
  
  betaPts <- c(getBrkPtdate(satBoost,0.1),getBrkPtdate(satBoost,0.5),getBrkPtdate(satBoost,0.9))
  gauss_satBoost <- getGaussShape(seq_along(satBoost),betaPts)
  target_gFRP <- gauss_satBoost*max(satBoost)
  
  ratioToGauss <- sum(satBoost)/sum(target_gFRP)
  ratioToGaussNew <- 0
  
  iIter <- 0
  
  # clouds/gap fill by iteration until convergence
  while (iIter <= 50 & iIter >= 0 & abs(ratioToGauss-ratioToGaussNew) != 0 & ratioToGaussNew < 0.8) {
    
    iIter <- iIter + 1    
    
    if (iIter == 1) {
      dFRP <- c(0,cFRP[-1]-cFRP[-nObs])
      tFRP <- dFRP*satFRPmax
    }
    
    frpPeak <- weighted.mean(seq_along(cFRP),dFRP,na.rm=T) # midpoint of fire season
    frpNPS <- get_d_nls(dFRP)

    betaPts <- c(getBrkPtdate(tFRP,0.1),getBrkPtdate(tFRP,0.5),getBrkPtdate(tFRP,0.9))
    gauss_satBoost <- getGaussShape(seq_along(satBoost),betaPts)
    target_gFRP <- gauss_satBoost*max(satBoost)
    ratioToGauss <- sum(tFRP)/sum(target_gFRP)
    
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
        
        neighbor_off1 <- neighbor_dFRP[midWindow]/max(neighbor_dFRP[windowBlock1])
        neighbor_off <- neighbor_dFRP[midWindow]/max(neighbor_dFRP)
        
        neighborSens_off <- neighbor_nps[midWindow]/max(neighbor_nps)*neighborSensAdj
        
        # for days at the edges of the distribution, use slopes that were
        # calculated excluding days outside the fire season
        inSlope <- slope
        if (iFRP < x1+nWindow) {inSlope <- slope1}
        if (iFRP > x2-nWindow) {inSlope <- slope2}
        
        # stop adjustments if criteria are met
        if (neighbor_dFRP[midWindow] == 0 | 
            iFRP < frpPeak & neighbor_off1 > 0 & neighbor_off1 > neighborSens_off |
            iFRP > frpPeak & neighbor_off > 0 & neighbor_off > neighborSens_off) {
          inSlope <- 0
        }
        
        # adjust cFRP at a rate proportional to (1 - detection rate) so that
        # cloudy/hazy days have priority
        cFRP[iFRP:nObs] <- cFRP[iFRP:nObs] + (inSlope * (1-wObsFrac[iFRP]))
        dFRP <- c(0,cFRP[-1]-cFRP[-nObs])
        tFRP <- dFRP*sum(satFRPmax)
        
        maxLimFRP <- apply(cbind(target_gFRP,satBoost),1,max)
        tFRP[tFRP>maxLimFRP] <- maxLimFRP[tFRP>maxLimFRP]
        ratioToGaussNew <- sum(tFRP)/sum(target_gFRP)
      }
    }
  }
  
  fireMat$adjFRP <- tFRP
  
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
  fireMat$diurnalBoost <- (fireMat$adjFRP*terraFrac*2.5/0.75 + fireMat$adjFRP*aquaFrac*1.5/0.75)/
    diurnalSat - fireMat$adjFRP
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
