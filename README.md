# SAGE-IGP

SAGE-IGP (Survey Constraints on FRP-based AGricultural Fire Emissions in the Indo-Gangetic Plain): agricultural fire emissions in north India 

![banner image](https://github.com/tianjialiu/SAGE-IGP/blob/main/docs/imgs/adjFRP_DMyr.png)

## Input Datasets
We use the following datasets:

* MODIS/MxD14A1 Fire Radiative Power, 1km
* MODIS/MCD14ML Active Fires
* VIIRS/VNP14IMGML Active Fires
* MODIS/MOD09GA Daily Surface Reflectance, 500m
* MODIS/MCD12Q1 Land Cover, 500 m

## Workflow

### Google Earth Engine
#### Preprocess input datasets at state-level
GEE/Preprocess/ > 

1. `MxD14A1_FRP.js`: export daily state-level MODIS FRP
2. `MxD14A1_FireMask.js`: export a MODIS fire mask to EE assets based on the aggregate of all fire detections
3. `VNP14IMGML_FRPboost.js`: export daily state-level VIIRS FRP boost relative to MODIS/Aqua FRP
4. `MCD14ML_SR.js`: export surface reflectance values associated with each MODIS fire detection for degree of cloudiness/haziness estimation
5. `MCD14ML_SR_fill.js`: export surface reflectance values associated with each MODIS fire detection for degree of cloudiness/haziness estimation (back-fill for months with no fires)

#### Preprocess input datasets on a 0.25°x0.25° grid
GEE/GriddedFRP/ > 

6. `MxD14A1_FRP_Grid.js`: export daily state-level MODIS FRP on a 0.25°x0.25° grid
7. `VNP14IMGML_FRPboost_Grid.js`: export daily state-level VIIRS FRP boost relative to MODIS/Aqua FRP on a 0.25° x 0.25° grid

### R
#### Preprocess output tables from GEE
R/processGEE/ >

1. `MCD14ML_SR_process.R`: preprocess MCD14ML SR files from EE
2. `MODIS_FRP_Grid_process.R`: preprocess MxD14A1 FRP files from EE
3. `VIIRS_FRPboost_Grid_process.R`: preprocess VIIRS FRP boost files from EE

#### Adjusted FRP algorithm and construction of SAGE-IGP Inventory
R/adjFRP/ > 

4. `cloud_gap_cutoff.R`: cloud/haze gap cutoff values
5. `viirs_boost_coef.R`: calculate VIIRS scaling factor for each state
6. `adjFRP_T1.R`: adjusted FRP for Tier 1 states - Punjab, Haryana
7. `adjFRP_T2.R`: adjusted FRP for Tier 2 states - UP, Bihar
8. `adjFRP_T3.R`: adjusted FRP for Tier 3 states - Rajasthan
9. `adjFRP_DM.R`: convert adjusted FRP to dry matter
10. `adjFRP_DMgridST.R`: grid state-level dry matter at 0.25°x0.25° resolution
11. `adjFRP_DMgrid.R`: combine gridded state-level dry matter
12. `adjFRP_DMaer.R`: combine gridded state-level dry matter, aerosols
13. `nc_adjFRP_DM.R`: make annual netCDF files for the SAGE-IGP inventory

## SAGE-IGP Dataset
Liu T., L.J. Mickley, S. Singh, M. Jain, R.S. DeFries, and M.E. Marlier (2020). SAGE-IGP agricultural fire emissions in north India, https://doi.org/10.7910/DVN/JUMXOL, Harvard Dataverse, V1

### Usage Notes
* Use agricultural emissions factors to convert dry matter (DM) to other chemical species, e.g. from [Andreae (2019, ACP)](https://doi.org/10.5194/acp-2019-303)
* Use `DMaer` for aerosol species and `DM` for all other chemical species

```
# ======================================
# Example R Script for Reading SAGE-IGP
# ======================================

library(ncdf4); library(raster)

# --------------------------
# Convert netCDF to raster
# --------------------------
# read SAGE-IGP netCDF for 2017
nc <- nc_open("adjFRP_Inv/nc/SAGE-IGP_daily_2017.nc")

# retrieve variables from netCDF
DMdaily <- ncvar_get(nc,"DM")

# total DM in 2017 from Sep-Dec
DMtotal <- apply(DMdaily,c(1,2),sum)

# read lat, lon
lat <- ncvar_get(nc,"lat")
lon <- ncvar_get(nc,"lon")

# extent of SAGE-IGP bounds
regionExtent <- extent(c(72,89,23,33))

# convert 'DMtotal' (a matrix) to a raster
DMtotal_ras <- raster(t(DMtotal)[rev(order(lat)),]) / 1e9 # convert from kg to Tg per grid cell
extent(DMtotal_ras) <- regionExtent # set extent of raster to regionExtent
crs(DMtotal_ras) <- crs(raster()) # EPSG:4326 - default lat/lon projection

# a simple plot of total DM in 2017
plot(DMtotal_ras)

# raster where each grid cell represents its area in sq. meters
area_m2 <- raster::area(DMtotal_ras) * 1e6

# convert Nov 1, 2017 emissions from kg to kg/m2/s
DMdaily_ras <- raster(t(DMdaily[,,62])[rev(order(lat)),])
DMdaily_ras <- DMtotal_ras / area_m2 / (24*24*60)
extent(DMdaily_ras) <- regionExtent
crs(DMdaily_ras) <- crs(raster())

# --------------------------------------
# Convert DM to other chemical species
# --------------------------------------
# for aerosol species, e.g. OC, BC, PM25, use 'DMaer'
# for all other species, use 'DM'

# retrieve variables from netCDF
DMdailyAer <- ncvar_get(nc,"DMaer")

# define emissions factors as g / kg DM, from Andreae (2019, ACP)
efs_andreae <- data.frame(OC=4.9,BC=0.42)

OCdaily <- DMdailyAer * efs_andreae$OC / 1e9 # daily OC, in Gg
BCdaily <- DMdailyAer * efs_andreae$BC / 1e9 # daily BC, in Gg
```

### GEOS-Chem - HEMCO Compatibility
SAGE-IGP can be combined with the GFEDv4s inventory by replacing GFEDv4s with SAGE-IGP agricultural fire emissions over north India. The resulting GFEDv4s+SAGE-IGP fire emissions follows the netCDF format of the original GFEDv4s emissions in HEMCO. Please contact me via email if you're interested in running GEOS-Chem with these files.

### Issues
* The original script `adjFRP_T1.R` used to produce Punjab/Haryana emissions in SAGE-IGP was accidentally overwritten and cannot be recovered. The updated `adjFRP_T1.R` yields the same mean budget of dry matter burned overall, but there may be slight differences in annual and daily variability.

## Publication
Liu T., L.J. Mickley, S. Singh, M. Jain, R.S. DeFries, and M.E. Marlier (2020). Crop residue burning practices across north India inferred from household survey data: bridging gaps in satellite observations. *Atmos. Environ. X* 8, 100091. https://doi.org/10.1016/j.aeaoa.2020.100091
