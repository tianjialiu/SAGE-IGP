# ----------------------------------------------
# Calculate VIIRS scaling factor for each state
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# ----------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')

xYears <- 2012:2018
xMonths <- 9:12
inState <- "Punjab"

modisFRPgrid <- read.csv(paste0("gee/MxD14A1/MxD14A1_FRP_",inState,".csv"))

fireMat <- blankDates(xMonths[1],xMonths[length(xMonths)],xYears,
                      NAcols=c("modisFRP","viirsFRPboost"))
fireMat <- fireMat[fireMat$Month %in% xMonths,]

ratioVM <- matrix(NA,length(xYears),2)
for (iYear in xYears) {
  for (iMonth in xMonths) {
    viirs_frp <- read.csv(paste0("gee/VNP14IMGML_FRPboost/VNP14IMGML_FRP_",inState,"_",iYear,"_",sprintf("%02d",iMonth),".csv"))
    fireMat$viirsFRPboost[which(fireMat$Year==iYear & fireMat$Month==iMonth)] <- viirs_frp$FRPboost
    
    monthday <- which(fireMat$Year==iYear & fireMat$Month==iMonth)
    
    for (iDay in seq_along(monthday)) {
      frp <- modisFRPgrid$MYD14A1[which(modisFRPgrid$YYYYMMDD==(iYear*1e4+iMonth*1e2+iDay))]
      fireMat$modisFRP[which(fireMat$Year==iYear & fireMat$Month==iMonth & fireMat$Day==iDay)] <- frp
    }
  }
  ratioVM[iYear-xYears[1]+1,] <- c(sum(fireMat$viirsFRPboost[fireMat$Year==iYear]),
                                   sum(fireMat$modisFRP[fireMat$Year==iYear]))
}

ratio_mean <- wt.mean(ratioVM[,1]/ratioVM[,2],ratioVM[,1]+ratioVM[,2])
ratio_sd <- wt.sd(ratioVM[,1]/ratioVM[,2],ratioVM[,1]+ratioVM[,2])

print(round(ratio_mean,2))
print(round(ratio_sd,2))
print(round(ratio_sd/ratio_mean*100))
