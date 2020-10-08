/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var indiaShp = ee.FeatureCollection("projects/GlobalFires/IndiaAgFires/IND_adm1"),
    gfedGrid = ee.FeatureCollection("projects/GlobalFires/GFEDv4poly"),
    myd14a1 = ee.ImageCollection("MODIS/006/MYD14A1");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// =========================================
// Export daily state-level VIIRS FRP boost
// relative to MODIS/Aqua FRP
// on a 0.25deg x 0.25deg grid
// =========================================
// @author: Tianjia Liu
// Last updated: May 9, 2020
// -----------------------------------------

// Input Parameters
var projFolder = 'projects/GlobalFires/';
var ST_NM = 'Punjab';

var sYear = 2012; var eYear = 2018;
var sMonth = 9; var eMonth = 12;
var satMODIS = 'A'; // Aqua: 'A' or Terra 'T'

// Global Parameters
var params = require('users/tl2581/SAGE-IGP:InputParams.js');
var proj = params.modis1km.projection();
var nDayMonthList = params.nDayMonthList;

// Region Boundaries
var Shp = indiaShp.filter(ee.Filter.eq('STATE',ST_NM.toUpperCase()));
var ShpGrid = gfedGrid.filterBounds(Shp).sort('id');
if (ST_NM == 'Rajasthan') {
  var ShpDist = ee.FeatureCollection(ee.List(['Ganganagar','Hanumangarh'])
    .map(params.filterDistricts));
  ShpGrid = gfedGrid.filterBounds(ShpDist).sort('id');
}
var ShpGridList = ShpGrid.toList(500,0);

if (satMODIS == 'A') {var satMODISsr = 'MYD09GA'; var satName = 'Aqua'}
if (satMODIS == 'T') {var satMODISsr = 'MOD09GA'; var satName = 'Terra'}

var addBuffer = function(feat) {return feat.buffer(proj.nominalScale())};
var doesIntersect = function(feat,rightFeat) {return feat.intersects(rightFeat)};

var nGridList = {
  'Punjab': 106,
  'Haryana': 103,
  'Uttar Pradesh': 428,
  'Bihar': 181,
  'Rajasthan': 52
};

// Note: may need to chunk exports to prevent timeouts
for (var iGrid = 0; iGrid < nGridList[ST_NM]; iGrid++) {
  
  var inShpGrid = ee.Feature(ShpGridList.get(iGrid));
  var gridID = inShpGrid.get('id').getInfo();
  var inShpGrid = inShpGrid.geometry();
  
  var fire_year = [];
  for (var inYear = sYear; inYear <= eYear; inYear++) {
    
    var days_of_month = nDayMonthList[inYear];
  
    for (var inMonth = sMonth; inMonth <= eMonth; inMonth++) {
      var inMonthStr = ee.Number(inMonth).format('%02d').getInfo();
      
      var fire_month = [];
      for (var inDay = 1; inDay <= days_of_month[inMonth-1]; inDay++) {
        
        var yyyymmdd = inYear*1e4+inMonth*1e2+inDay;
        
        var mcd14ml_day = ee.FeatureCollection(projFolder + 'MCD14ML/MCD14ML_' + inYear + '_' + inMonthStr)
          .filter(ee.Filter.eq('sat',satMODIS)).filter(ee.Filter.eq('type',0))
          .filter(ee.Filter.eq('dn','D'))
          .filter(ee.Filter.eq('YYYYMMDD',yyyymmdd))
          .filterBounds(Shp).filterBounds(inShpGrid);
          
        var vnp14ml_day = ee.FeatureCollection(projFolder + 'VNP14IMGML/VNP14IMGML_' + inYear + '_' + inMonthStr)
          .filter(ee.Filter.eq('type',0))
          .filter(ee.Filter.lt('HHMM',1000))
          .filter(ee.Filter.eq('YYYYMMDD',yyyymmdd))
          .filterBounds(Shp).filterBounds(inShpGrid);
        
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
          FRP: viirsBuf.aggregate_sum('FRP'),
          FRPboost: viirsIntMODIS.aggregate_sum('FRP')
        });
      }
      
      var fire_year = ee.FeatureCollection(fire_year).merge(ee.FeatureCollection(fire_month));
    }
  }

  Export.table.toDrive({
    collection: ee.FeatureCollection(fire_year),
    folder: 'VNP14IMGML_FRPboost_Grid',
    description: 'VNP14IMGML_FRP_' + ST_NM.split(' ').join('_') + '_Grid_' + gridID,
    selectors: ['YYYYMMDD','FRP','FRPboost']
  });
}