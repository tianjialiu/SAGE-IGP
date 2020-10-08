/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var indiaShp = ee.FeatureCollection("projects/GlobalFires/IndiaAgFires/IND_adm1"),
    mcd12q1 = ee.ImageCollection("MODIS/006/MCD12Q1"),
    fireMask = ee.Image("projects/GlobalFires/IndiaAgFires/FireMask_MxD14A1");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// ==============================================
// Export surface reflectance values
// associated with each MODIS fire detection
// for degree of cloudiness/haziness estimation
// ==============================================
// @author: Tianjia Liu
// Last updated: August 12, 2020
// ----------------------------------------------

// Input Parameters
var projFolder = 'projects/GlobalFires/';
var ST_NM = 'Punjab';

var sYear = 2003; var eYear = 2018;
var sMonth = 9; var eMonth = 12;
var satMODIS = 'T'; // Aqua: 'A' or Terra 'T'

// Global Parameters
var params = require('users/tl2581/SAGE-IGP:InputParams.js');
var proj = params.modis500m.projection();
var nDayMonthList = params.nDayMonthList;

// Region Boundaries
var Shp = indiaShp.filter(ee.Filter.eq('STATE',ST_NM.toUpperCase()));

if (satMODIS == 'A') {var satMODISsr = 'MYD09GA'; var satName = 'Aqua'}
if (satMODIS == 'T') {var satMODISsr = 'MOD09GA'; var satName = 'Terra'}

for (var inYear = sYear; inYear <= eYear; inYear++) {
  
  var days_of_month = nDayMonthList[inYear];
  var filterYr = ee.Filter.calendarRange(inYear,inYear,'year');
  
  // Agricultural Mask
  if (inYear < 2018) {
    var agMask = ee.Image(mcd12q1.filter(filterYr).first())
      .select('LC_Type2').eq(12)
      .reproject({crs: proj, scale: proj.nominalScale()});
  } else {
    var agMask = ee.Image(mcd12q1
      .filter(ee.Filter.calendarRange(2018,2018,'year')).first())
      .select('LC_Type2').eq(12)
      .reproject({crs: proj, scale: proj.nominalScale()});
  }
  
  for (var inMonth = sMonth; inMonth <= eMonth; inMonth++) {
    var inMonthStr = ee.Number(inMonth).format('%02d').getInfo();
    
    var fire_sr_month = []; var fire_sr_area_month = [];
    for (var inDay = 1; inDay <= days_of_month[inMonth-1]; inDay++) {
      
      var yyyymmdd = inYear*1e4+inMonth*1e2+inDay;
      
      var mcd14ml_day = ee.FeatureCollection(projFolder + 'MCD14ML/MCD14ML_' + inYear + '_' + inMonthStr)
        .filter(ee.Filter.eq('sat',satMODIS)).filter(ee.Filter.eq('type',0))
        .filter(ee.Filter.eq('dn','D'))
        .filter(ee.Filter.eq('YYYYMMDD',yyyymmdd))
        .filterBounds(Shp);

      var mxd09ga_day = ee.Image(ee.ImageCollection('MODIS/006/' + satMODISsr)
        .filter(ee.Filter.calendarRange(inYear,inYear,'year'))
        .filter(ee.Filter.calendarRange(inMonth,inMonth,'month'))
        .filter(ee.Filter.calendarRange(inDay,inDay,'day_of_month')).first())
        .select('sur_refl_b0[1-7]*').multiply(0.0001)
        .updateMask(agMask).updateMask(fireMask);
      
      var fire_sr_incm = mxd09ga_day.select('sur_refl_b01').clamp(0,1)
        .reduceRegions({
          collection: Shp,
          reducer: ee.Reducer.fixedHistogram(0,1,100),
          crs: proj,
          scale: proj.nominalScale()
        });
      
      fire_sr_incm = ee.Feature(fire_sr_incm.first()).get('histogram');
      fire_sr_area_month[inDay-1] = ee.Feature(null,{YYYYMMDD: yyyymmdd, histogram: fire_sr_incm});
      
      var fire_sr = mxd09ga_day.reduceRegions({
        collection: mcd14ml_day,
        reducer: ee.Reducer.max().unweighted(),
        crs: proj,
        scale: proj.nominalScale()
      });
      
      fire_sr_month[inDay-1] = fire_sr;
    }
    
    fire_sr_month = ee.FeatureCollection(fire_sr_month).flatten();
    fire_sr_area_month = ee.FeatureCollection(fire_sr_area_month);
    
    Export.table.toDrive({
      collection: fire_sr_month,
      folder: 'MCD14ML_SR',
      description: 'MCD14ML_SR_' + satName + '_' + ST_NM.replace(' ','_') + '_' + inYear + '_' + inMonthStr,
      selectors: ['YYYYMMDD','HHMM','conf','sat','FRP',
        'sur_refl_b01','sur_refl_b02','sur_refl_b03',
        'sur_refl_b04','sur_refl_b05','sur_refl_b06',
        'sur_refl_b07','.geo']
    });
    
    Export.table.toDrive({
      collection: fire_sr_area_month,
      folder: 'MCD14ML_SR',
      description: 'MCD14ML_SR_Area_' + satName + '_' + ST_NM.split(' ').join('_') + '_' + inYear + '_' + inMonthStr,
      selectors: ['YYYYMMDD','histogram']
    });
  }
}
