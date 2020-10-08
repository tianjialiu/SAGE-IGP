# --------------------------------------------
# Cloud/haze gap cutoff values
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# --------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')
setwd(file.path(geeHome,"MCD14ML_SR"))

xYears <- 2003:2018
xMonths <- 9:12
xSat <- c("Aqua","Terra")
xStates <- c("Punjab","Haryana","Uttar_Pradesh","Bihar")
band_name <- "sur_refl_b01"

blankMat <- blankDates(1,12,xYears,NAcols=xSat)
blankMat <- blankMat[blankMat$Month%in%xMonths,]
hist_prob <- array(NA,dim=c(100,dim(blankMat)[1],length(xSat)))
fire_sum <- matrix(NA,dim(blankMat)[1],length(xSat))

for (iSat in seq_along(xSat)) {
  inSat <- xSat[iSat]
  
  for (iYear in xYears) {
    blankYr <- blankDates(xMonths[1],xMonths[length(xMonths)],iYear)
    
    for (iMonth in xMonths) {
      
      fire_srSt <- list()
      for (iState in seq_along(xStates)) {
        fire_sr_name <- paste0("MCD14ML_SR_",inSat,"_",xStates[iState],"_",
                               iYear,"_",sprintf("%02d",iMonth),".csv")
        fire_srSt[[iState]] <- read.csv(fire_sr_name)
      }
      fire_sr <- do.call(rbind,fire_srSt)
      
      monthday <- which(blankYr$Month==iMonth)
      
      for (iDay in seq_along(monthday)) {
        dayIdx <- which(blankMat$Year==iYear & blankMat$Month==iMonth & blankMat$Day==iDay)
        
        fire_sr_day <- fire_sr[which(as.numeric(substr(fire_sr$YYYYMMDD,7,8)) == iDay),]
        fire_sr_day <- fire_sr_day[which(fire_sr_day[,band_name]>0),]
        
        if (dim(fire_sr_day)[1] > 0) {
          fire_sum[dayIdx,iSat] <- sum(fire_sr_day$FRP)
          fire_sr_day_band <- fire_sr_day[,band_name]
          fire_sr_day_band[fire_sr_day_band > 1] <- 1
          fire_sr_day_band[fire_sr_day_band < 0] <- 0
          hist_fire_sr_day <- hist(fire_sr_day_band,breaks=seq(0,1,0.01),plot=F)
          hist_prob[,dayIdx,iSat] <- hist_fire_sr_day$counts/max(hist_fire_sr_day$counts)
        }
      }
    }
  }
}

hist_prob[is.na(hist_prob)] <- 0
fire_sum[is.na(fire_sum)] <- 0

hist_prob_sat <- cbind(hist_prob[,,1],hist_prob[,,2])
fire_sum_sat <- c(fire_sum[,1],fire_sum[,2])

whist_prob <- apply(hist_prob_sat,1,function(x) {weighted.mean(x,fire_sum_sat)})
whist_prob <- whist_prob/max(whist_prob)

histSteps <- seq(0.005,0.995,0.01)
plot(histSteps,whist_prob,t="l")

cutoff_st <- histSteps[which(whist_prob==max(whist_prob))]
cutoff_end <- histSteps[which(round(whist_prob,2)==0 & histSteps > cutoff_st)[1]]

x <- histSteps[which(round(whist_prob,2)>0 & histSteps > cutoff_st)]
y <- whist_prob[which(round(whist_prob,2)>0 & histSteps > cutoff_st)]
x2 <- x^2

model_lm <- lm(log(y)~x+x2)
b0 <- coef(model_lm)[1]
b1 <- coef(model_lm)[2]
b2 <- coef(model_lm)[3]

print(model_lm)
print(c(cutoff_st,cutoff_end))
abline(v=c(cutoff_st,cutoff_end))
