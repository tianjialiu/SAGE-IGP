# --------------------------------------------
# Convert adjusted FRP to dry matter
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# --------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')

xYears <- 2003:2018
xStates <- c("Haryana","Punjab","Uttar_Pradesh","Bihar","Rajasthan")

alpha <- 0.41 # FRE to DM, kg/MJ
t_fire <- 60*60*0.5 # duration of fire, sec

adjDM_st <- matrix(NA,length(xYears),length(xStates))

for (iState in seq_along(xStates)) {
  
  setwd(file.path(tabHome,"MODISadjFRP"))
  adjDMdaily <- matrix(NA,122,length(xYears))
  for (iYear in seq_along(xYears)) {
    fireMat <- read.csv(paste0("MODISadjFRP_",xStates[iState],"_",xYears[iYear],".csv"))
    adjDMdaily[,iYear] <- fireMat$adjFRP*t_fire*alpha/1e9
  }
  colnames(adjDMdaily) <- paste0("Y",xYears)
  
  setwd(file.path(tabHome,"MODISadjFRP_DM"))
  write.table(adjDMdaily,paste0("MODISadjFRP_DM_",xStates[iState],".csv"),sep=",",row.names=F)
  
  adjDM_st[,iState] <- colSums(adjDMdaily)
}

print(round(mean(rowSums(adjDM_st)),1))
print(round(sd(rowSums(adjDM_st)),1))
