/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var indiaShp = ee.FeatureCollection("projects/GlobalFires/IndiaAgFires/IND_adm1");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// ==============================================
// Export a MODIS fire mask to EE assets
// based on the aggregate of all fire detections
// ==============================================
// @author: Tianjia Liu
// Last updated: August 11, 2020
// ----------------------------------------------

// Global Parameters
var params = require('users/embrslab/SAGE-IGP:InputParams.js');
var proj = params.modis1km.projection();
var sYear = params.sYear;
var eYear = params.eYear;

var outputRegion = indiaShp.geometry().bounds();

// Calculate FRP
var getFRP = function(image) {
  var FireMask = image.select('FireMask').gte(7);
  var FRP = image.select('MaxFRP').unmask(0).multiply(0.1);

  return FRP.updateMask(FireMask);
};

var terraMask = ee.ImageCollection('MODIS/006/MOD14A1')
  .filter(ee.Filter.calendarRange(sYear,eYear,'year'))
  .filter(ee.Filter.calendarRange(9,12,'month'))
  .map(getFRP).max().gt(0).unmask(0);

var aquaMask = ee.ImageCollection('MODIS/006/MYD14A1')
  .filter(ee.Filter.calendarRange(sYear,eYear,'year'))
  .filter(ee.Filter.calendarRange(9,12,'month'))
  .map(getFRP).max().gt(0).unmask(0);
  
var fireMask = terraMask.add(aquaMask).gt(0).selfMask()
  .clip(indiaShp).reproject({crs: proj, scale: proj.nominalScale()});

Map.setCenter(77,30,7);
Map.addLayer(fireMask, {palette: '#FF0000'});

Export.image.toAsset({
  image: fireMask,
  region: outputRegion,
  description: 'fireMask',
  assetId: 'projects/GlobalFires/IndiaAgFires/FireMask_MxD14A1',
  crs: 'SR-ORG:6974',
  crsTransform: [926.625433055833,0,-20015109.354,0,926.6254330558334,10007554.677003],
  maxPixels: 1e12
});
