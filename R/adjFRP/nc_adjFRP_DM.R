# -----------------------------------------------------
# Make annual netCDF files for the SAGE-IGP inventory
# @author Tianjia Liu (tianjialiu@g.harvard)
# Last updated: September 17, 2020
# -----------------------------------------------------

source('~/Google Drive/India/methods/IndiaSurvey/globalParams.R')

doVIIRS <- F

if (doVIIRS == T) {
  xYears <- 2012:2018
} else {
  xYears <- 2003:2018
}

file_suffix <- "_"
if (doVIIRS == T) {file_suffix <- "VIIRS_"}

for (inYear in xYears) {
  
  # define time dimension
  time_out <- seq(as.Date(paste0(inYear,"-09-01")),as.Date(paste0(inYear,"-12-31")),"day")
  time_out <- as.POSIXct(paste(time_out,"00:00:00 UTC"),tz="UTC")
  time_out <- floor(unclass(time_out)) 
  
  # read raster file
  setwd(file.path(rasHome,"adjFRP_Inv/"))

  inFileName <- paste0("adjFRP_DM_daily",file_suffix,inYear)
  outFileName <- paste0("SAGE-IGP_daily",file_suffix,inYear)
  
  DMyr <- stack(paste0("DM/",inFileName,".tif"))
  DMaer_yr <- stack(paste0("DMaer/adjFRP_DMaer_daily",file_suffix,inYear,".tif"))
  
  coordsRas <- coordinates(DMyr)
  xLon <- unique(coordsRas[,1]); xLat <- unique(coordsRas[,2])
  
  areaRas <- raster::area(DMyr[[1]])*1e6 # m2
  areaMat <- t(as.matrix(areaRas))
  
  DMmat <- array(NA,dim=c(length(xLon), length(xLat), length(time_out)))
  DMaer_mat <- DMmat
  for (iDay in seq_along(time_out)) {
    DMmat[,,iDay] <- t(as.matrix(DMyr[[iDay]]))*1e9 # kg
    DMaer_mat[,,iDay] <- t(as.matrix(DMaer_yr[[iDay]]))*1e9 # kg
  }
  
  # netCDF output
  # xy dimensions in lat/lon
  xdim <- ncdim_def('lon', 'degrees_east', xLon)
  ydim <- ncdim_def('lat', 'degrees_north', xLat)
  
  tdim <- ncdim_def('time', 'seconds since 1970-01-01',
                    as.integer(time_out))
  
  dm_var <- ncvar_def('DM', 'kg', list(xdim, ydim, tdim))
  dmaer_var <- ncvar_def('DMaer', 'kg', list(xdim, ydim, tdim))
  area_var <- ncvar_def('area', 'm2', list(xdim, ydim))
  
  # Projection specific xy definitions
  setwd(file.path(rasHome,"adjFRP_Inv/nc/"))
  nc <- nc_create(paste0(outFileName,".nc"), list(dm_var, dmaer_var, area_var), force_v4 = T)
  
  ncatt_put(nc, 'lon', 'standard_name', 'longitude')
  ncatt_put(nc, 'lon', 'long_name', 'longitude')
  ncatt_put(nc, 'lat', 'standard_name', 'latitude')
  ncatt_put(nc, 'lat', 'long_name', 'latitude')
  
  # Insert DM data
  ncvar_put(nc, dm_var, DMmat)
  ncvar_put(nc, dmaer_var, DMaer_mat)
  ncvar_put(nc, area_var, areaMat)
  
  ncatt_put(nc, 'time', 'standard_name', 'time')
  ncatt_put(nc, 'time', 'calendar', 'standard')
  
  ncatt_put(nc, 'DM', 'standard_name', 'DM')
  ncatt_put(nc, 'DM', 'long_name', 'Dry Matter burned')
  
  ncatt_put(nc, 'DMaer', 'standard_name', 'DMaer')
  ncatt_put(nc, 'DMaer', 'long_name', 'Dry Matter burned, aerosols')
  
  ncatt_put(nc, 'area', 'long_name', 'grid cell area')
  
  ncatt_put(nc, 0, 'crs', '+proj=longlat')
  ncatt_put(nc, 0, 'crs_format', 'PROJ.4')
  ncatt_put(nc, 0, 'citation', 'Liu et al. Crop residue burning practices across north India inferred from household survey data: bridging gaps in satellite observations.')
  ncatt_put(nc, 0, 'title', 'SAGE-IGP agricultural fire emissions in north India')
  ncatt_put(nc, 0, 'time_created', format(Sys.time(), tz = 'UTC'))
  
  nc_close(nc)
  
  timestamp(prefix=paste("Year",inYear,": ##------ "))
}
