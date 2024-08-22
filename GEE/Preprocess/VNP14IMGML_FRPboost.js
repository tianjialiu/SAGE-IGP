/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var indiaShp = ee.FeatureCollection("projects/GlobalFires/IndiaAgFires/IND_adm1"),
    mcd12q1 = ee.ImageCollection("MODIS/006/MCD12Q1"),
    myd14a1 = ee.ImageCollection("MODIS/006/MYD14A1");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// =========================================
// Export daily state-level VIIRS FRP boost
// relative to MODIS/Aqua FRP
// =========================================
// @author: Tianjia Liu
// Last updated: August 17, 2020
// -----------------------------------------

// Input Parameters
var projFolder = 'projects/GlobalFires/';
var ST_NM = 'Punjab';

var sYear = 2012; var eYear = 2018;
var sMonth = 9; var eMonth = 12;
var satMODIS = 'A'; // Aqua: 'A' or Terra 'T'

// Global Parameters
var params = require('users/embrslab/SAGE-IGP:InputParams.js');
var proj = params.modis1km.projection();
var bufferSize = proj.nominalScale();
var nDayMonthList = params.nDayMonthList;
var monthStrList = params.monthStrList;

// Region Boundaries
var Shp = indiaShp.filter(ee.Filter.eq('STATE',ST_NM.toUpperCase()));
if (ST_NM == 'Rajasthan') {
  Shp = ee.FeatureCollection(ee.List(['Ganganagar','Hanumangarh'])
    .map(params.filterDistricts)).union();
}

if (satMODIS == 'A') {var satMODISsr = 'MYD09GA'; var satName = 'Aqua'}
if (satMODIS == 'T') {var satMODISsr = 'MOD09GA'; var satName = 'Terra'}

var addBuffer = function(feat) {return feat.buffer(bufferSize)};
var doesIntersect = function(feat,rightFeat) {return feat.intersects(rightFeat)};

for (var inYear = sYear; inYear <= eYear; inYear++) {
  
  var days_of_month = nDayMonthList[inYear];
  
  for (var inMonth = sMonth; inMonth <= eMonth; inMonth++) {
    var inMonthStr = monthStrList[inMonth-1];
    
    var fire_month = [];
    for (var inDay = 1; inDay <= days_of_month[inMonth-1]; inDay++) {
      
      var yyyymmdd = inYear*1e4+inMonth*1e2+inDay;
      
      var mcd14ml_day = ee.FeatureCollection(projFolder + 'MCD14ML/MCD14ML_' + inYear + '_' + inMonthStr)
        .filter(ee.Filter.eq('sat',satMODIS)).filter(ee.Filter.eq('type',0))
        .filter(ee.Filter.eq('dn','D'))
        .filter(ee.Filter.eq('YYYYMMDD',yyyymmdd))
        .filterBounds(Shp);
        
      var vnp14ml_day = ee.FeatureCollection(projFolder + 'VNP14IMGML/VNP14IMGML_' + inYear + '_' + inMonthStr)
        .filter(ee.Filter.eq('type',0))
        .filter(ee.Filter.lt('HHMM',1000))
        .filter(ee.Filter.eq('YYYYMMDD',yyyymmdd))
        .filterBounds(Shp);
      
      var modisImg = mcd14ml_day.reduceToImage({
        properties: ['HHMM'],
        reducer: ee.Reducer.max()
      }).reproject({crs: proj, scale: proj.nominalScale()});
      
      var viirsBuf = vnp14ml_day.map(addBuffer);

      var viirsIntMODIS = modisImg.unmask(0).gt(0).reduceRegions({
        collection: viirsBuf,
        reducer: ee.Reducer.max(),
        crs: proj,
        scale: proj.nominalScale()
      }).filter(ee.Filter.eq('max',0));

      fire_month[inDay-1] = ee.Feature(null,{
        YYYYMMDD: yyyymmdd,
        HHMM: vnp14ml_day.reduceColumns(ee.Reducer.mode(),['HHMM']).get('mode'),
        FRP: vnp14ml_day.aggregate_sum('FRP'),
        FRPboost: viirsIntMODIS.aggregate_sum('FRP')
      });
      
    }

    Export.table.toDrive({
      collection: ee.FeatureCollection(fire_month),
      folder: 'VNP14IMGML_FRPboost',
      description: 'VNP14IMGML_FRP_' + ST_NM.split(' ').join('_') + '_' + inYear + '_' + inMonthStr,
      selectors: ['YYYYMMDD','HHMM','FRP','FRPboost']
    });
  }
}
