/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var indiaShp = ee.FeatureCollection("projects/GlobalFires/IndiaAgFires/IND_adm1"),
    indiaShpDist = ee.FeatureCollection("projects/GlobalFires/IndiaAgFires/IND_adm2"),
    mod14a1 = ee.ImageCollection("MODIS/006/MOD14A1"),
    mod09ga = ee.ImageCollection("MODIS/006/MOD09GA"),
    mod44b = ee.ImageCollection("MODIS/006/MOD44B");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// ====================================
// Input parameters for SAGE-IGP
// ====================================
// @author: Tianjia Liu
// Last updated: May 9, 2020
// ------------------------------------

// Time Parameters
var sYear = 2003;
var eYear = 2018;

var days_of_month_nonleap = [31,28,31,30,31,30,31,31,30,31,30,31];
var days_of_month_leap = [31,29,31,30,31,30,31,31,30,31,30,31];

var nDayList = {
  2000: 366, 2001: 365, 2002: 365, 2003: 365,
  2004: 366, 2005: 365, 2006: 365, 2007: 365,
  2008: 366, 2009: 365, 2010: 365, 2011: 365,
  2012: 366, 2013: 365, 2014: 365, 2015: 365,
  2016: 366, 2017: 365, 2018: 365, 2019: 365,
};

var nDayMonthList = {
  2000: days_of_month_leap, 2001: days_of_month_nonleap, 2002: days_of_month_nonleap, 2003: days_of_month_nonleap,
  2004: days_of_month_leap, 2005: days_of_month_nonleap, 2006: days_of_month_nonleap, 2007: days_of_month_nonleap,
  2008: days_of_month_leap, 2009: days_of_month_nonleap, 2010: days_of_month_nonleap, 2011: days_of_month_nonleap,
  2012: days_of_month_leap, 2013: days_of_month_nonleap, 2014: days_of_month_nonleap, 2015: days_of_month_nonleap,
  2016: days_of_month_leap, 2017: days_of_month_nonleap, 2018: days_of_month_nonleap, 2019: days_of_month_nonleap,
};

// Projections
var modis250m = ee.Image(mod44b.select('Percent_NonVegetated').first());
var modis500m = ee.Image(mod09ga.select('sur_refl_b07').first());
var modis1km = ee.Image(mod14a1.select('MaxFRP').first());

// Global Functions
var getQABits = function(image, start, end, newName) {
  // Compute the bits we need to extract
  var pattern = 0;
  for (var i = start; i <= end; i++) {
     pattern += Math.pow(2, i);
  }
  return image.select([0], [newName])
                .bitwiseAnd(pattern)
                .rightShift(start);
};

var filterDistricts = function(district) {
   return indiaShpDist.filterMetadata('DISTRICT','equals',district).first();
};

var filterStates = function(state) {
   return indiaShp.filterMetadata('STATE','equals',state).first();
};

// Export Variables
exports.getQABits = getQABits;
exports.filterDistricts = filterDistricts;
exports.filterStates = filterStates;
exports.modis250m = modis250m;
exports.modis500m = modis500m;
exports.modis1km = modis1km;
exports.sYear = sYear;
exports.eYear = eYear;
exports.nDayMonthList = nDayMonthList;
exports.nDayList = nDayList;