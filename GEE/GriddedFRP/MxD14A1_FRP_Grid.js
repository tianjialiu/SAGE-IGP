/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var indiaShp = ee.FeatureCollection("projects/GlobalFires/IndiaAgFires/IND_adm1"),
    gfedGrid = ee.FeatureCollection("projects/GlobalFires/GFEDv4poly"),
    mcd12q1 = ee.ImageCollection("MODIS/006/MCD12Q1"),
    mod14a1 = ee.ImageCollection("MODIS/006/MOD14A1"),
    myd14a1 = ee.ImageCollection("MODIS/006/MYD14A1");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// ====================================
// Export daily state-level MODIS FRP
// on a 0.25deg x 0.25deg grid
// ====================================
// @author: Tianjia Liu
// Last updated: August 17, 2020
// ------------------------------------

// Input Parameters
var ST_NM = 'Punjab';

// Global Parameters
var params = require('users/embrslab/SAGE-IGP:InputParams.js');
var proj = params.modis1km.projection();
var sYear = params.sYear;
var eYear = params.eYear;
var nDayList = params.nDayList;
var gridIDsList = params.gridIDsList;

// Region Boundaries
var Shp = indiaShp.filter(ee.Filter.eq('STATE',ST_NM.toUpperCase()));
var ShpGrid = gfedGrid.filterBounds(Shp).sort('id');
if (ST_NM == 'Rajasthan') {
  var ShpDist = ee.FeatureCollection(ee.List(['Ganganagar','Hanumangarh'])
    .map(params.filterDistricts));
  ShpGrid = gfedGrid.filterBounds(ShpDist).sort('id');
}
var ShpGridList = ShpGrid.toList(500);

// Calculate FRP
var getFRP = function(image) {
  var FireMask = image.select('FireMask').gte(7);
  var FRP = image.select('MaxFRP').unmask(0).multiply(0.1);
  var frpMasked = FRP.updateMask(agMask).updateMask(FireMask);
  return frpMasked;
};

var getFRPregion = function(inDay) {
  
  var filterJDay = ee.Filter.calendarRange(inDay,inDay,'day_of_year');
  var mod14a1Day = mod14a1Yr.filter(filterJDay);
  var myd14a1Day = myd14a1Yr.filter(filterJDay);
  
  var n_mod14a1 = mod14a1Day.size();
  mod14a1Day = getFRP(mod14a1Day.first());
  
  var n_myd14a1 = myd14a1Day.size();
  myd14a1Day = getFRP(myd14a1Day.first());
  
  var emptyFRP = ee.Image(0).toFloat().rename('MaxFRP');
  mod14a1Day = ee.Algorithms.If(n_mod14a1.gt(0),mod14a1Day,emptyFRP);
  myd14a1Day = ee.Algorithms.If(n_myd14a1.gt(0),myd14a1Day,emptyFRP);

  var combinedDay = ee.ImageCollection([
    ee.Image(mod14a1Day).toFloat(),
    ee.Image(myd14a1Day).toFloat()
  ]).sum();
  
  var frpDay = ee.Image(mod14a1Day).rename('MOD14A1')
    .addBands(ee.Image(myd14a1Day).rename('MYD14A1'))
    .addBands(combinedDay.rename('MxD14A1'))
    .reproject({crs: proj, scale: proj.nominalScale()})
    .clip(Shp);
    
  var frpRegionDay = frpDay
    .reduceRegions({
      reducer: ee.Reducer.sum(),
      collection: inShpGrid,
      crs: proj,
      scale: proj.nominalScale(),
    }).first();
  
  var inDate = ee.Date.fromYMD(inYear,1,1)
    .advance(ee.Number(inDay).subtract(1),'day');
    
  var year = inDate.get('Year');
  var month = inDate.get('Month');
  var day = inDate.get('Day');
  var yyyymmdd = year.multiply(1e4).add(month.multiply(1e2)).add(day);
    
  return ee.Feature(null,{
    YYYYMMDD: year.multiply(1e4).add(month.multiply(1e2)).add(day),
    Year: year, Month: month, Day: day
  }).copyProperties(frpRegionDay,['MOD14A1','MYD14A1','MxD14A1']);
};

var nGridList = {
  'Punjab': 106,
  'Haryana': 103,
  'Uttar Pradesh': 428,
  'Bihar': 181,
  'Rajasthan': 52
};

var gridIDs = gridIDsList[ST_NM];

// Note: may need to chunk exports to prevent timeouts
for (var iGrid = 0; iGrid < nGridList[ST_NM]; iGrid++) {
  var gridID = gridIDs[iGrid];
  var inShpGrid = ee.Feature(ShpGridList.get(iGrid));

  var FRPtable = [];
  for (var inYear = sYear; inYear <= eYear; inYear++) {
    var filterYr = ee.Filter.calendarRange(inYear,inYear,'year');
    var mod14a1Yr = mod14a1.filter(filterYr);
    var myd14a1Yr = myd14a1.filter(filterYr);
  
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
    
    var nDay = nDayList[inYear];
    
    var FRPyr = ee.FeatureCollection(ee.List.sequence(1,nDay,1).map(getFRPregion));
    var FRPtable = ee.FeatureCollection(FRPtable).merge(FRPyr);
  }
  
  Export.table.toDrive({
    collection: ee.FeatureCollection(FRPtable),
    description: 'MxD14A1_FRP_' + ST_NM.split(' ').join('_') + '_Grid_' + gridID,
    selectors: ['YYYYMMDD','Year','Month','Day',
      'MOD14A1','MYD14A1','MxD14A1']
  });
}
