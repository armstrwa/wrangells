'''
New script for swath profiler that does not require GRASS
Simplified to bare bones, just input filenames of raster and polyline to sample
begun 16 nov 2016
William Armstrong

Inputs: line geometry, raster
Outputs: raster swath profiled along line geometry

'''

## -------------- ##
## IMPORT MODULES ##
## -------------- ##

import os
import sys
import numpy as np
import wrangellsGeospatialTools as wgt 
import scipy as sci 
import csv
import glob
import matplotlib.pyplot as plt
import json
from shapely import geometry
import fiona
from osgeo import ogr, osr, gdal
from pyproj import Proj, transform
from fiona.crs import from_epsg
from gdalconst import *
import datetime as dt
import time

## --------- ##
## FUNCTIONS ##
## --------- ##

# functions loaded from wrangellsGeospatialTools.pyc

## ----------- ##
## USER INPUTS ##
## ----------- ##

shapefileIn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/qgis/keinholz_centerlines_paredDownToStudyArea.shp' 
shapefileEPSG = 3338 # spatial reference code for profile to sample
rasterIn = '/Users/wiar9509/Google Drive/wrangells/corr/LC80640172014057LGN0_LC80640172014089LGN0_2014-073_32_24_40_10_vv_filt_corr0.15_delcorr0.05_epsg3338.tif'
dDist = 500 # interpolation distance along transect line
boxWidth = 1000 # swath profiler width across transect line
# featNum = 66 for rootCL
featNum = 66 # if shp is multiline, use id number to choose 1 line 
outDirectory = '/Users/wiar9509/' # output directory

## -------- ##
## ANALYSIS ##
## -------- ##

# Bring in feature from shapefile
shp = fiona.open(shapefileIn)
originalSrsProj = shp.crs	
feat = shp[featNum]
shapelyLineGeom = geometry.shape(feat['geometry'])

# Make shapefile for swath profiling
xInt,yInt = wgt.interpolateShapelyLineGeom(shapelyLineGeom,dDist) # interpolate vertices using shapely

wgt.makePolylineShpFromCoords(xInt,yInt,shapefileEPSG,outDirectory+'tmpSwathSamplePoints.shp')
slopes=wgt.calculateSlopes(xInt,yInt) # calculate "slope" (azimuth) of line segments
perpSlopes=wgt.calculatePerpendicularSlopes(slopes) # calculate perpendicular "slope" (azimuth)
allBoxVerts = wgt.makeAllBoxVertices(xInt,yInt,boxWidth) # get box coordinates for each interpolated point, with a specified box width (cross-transect coord)

# Create the spatial reference
srs = osr.SpatialReference()
srs.ImportFromEPSG(shapefileEPSG) # This hard-coded to UTM 7N w/ WGS84 datum. Would be nice to have this defined based of input shapefile

wgt.makeMultipolygonFromBoxVertices(allBoxVerts,srs,outDirectory + 'tmpSamplePolygons.shp') # Make shapefile for swath profiling

profileCoords = wgt.readLineCoordsFromMultiLineShp(shapefileIn,featNum) # get coordinates for PRE-INTERPOLATED centerline profile
geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']} # make shapely-style geometry dict
profGeom = geometry.shape(geomToPass) # make into a shapely geometry for testing intersection in wgt.getPathRowList()

dist = wgt.distanceAlongLine(xInt,yInt) # get distance along line
dist = np.delete(dist,0) # delete the first entry, which we don't have a value for b/c don't have a slope for the first entry (b/c need to points), so interpolator doesn't make a 0,0 point .. this quick and dirty

interpProfileCoords = wgt.readVertexCoordsFromPolyline(outDirectory + 'tmpSwathSamplePoints.shp') # read vertex coordinates from interpolated line

# Get information about correlation
corrDict = wgt.getCorrelationInfoFromFn(rasterIn)			
raster = gdal.Open(rasterIn) # open one correlation
xyTfm = wgt.createTransformation(rasterIn) # transformation parameters for converting raster pixels to ground coordinates

# Rename spatial variables for writing to profile_data dict
sampled_on_date_iso = dt.datetime.now().isoformat()
sample_pts_lon_lat = interpProfileCoords['samplePtsLatLon']
sample_pts_frontdist = dist
sample_pts_PS = interpProfileCoords['samplePtsXy']
profile_shapefile_dir, trash , profile_shapefile = shapefileIn.rpartition('/')
profile_feature = featNum 
profile_shapefile_projection_Proj = originalSrsProj

# Swath profile
try:
	stats=wgt.zonal_stats(outDirectory + 'tmpSamplePolygons.shp', rasterIn, nodata_value=-1) # swath profile along transect
	statsDict = wgt.readZonalStatsOutput(stats) # read swath profile output
	
	profile_dict = wgt.getCorrelationInfoFromFn(rasterIn) # strip filename for correlation parameters
	profile_dict['mean'] = statsDict['means'] # add mean speed (duplicate, but just to be obvious in future
	profile_dict['max'] = statsDict['maxs'] # add maximum speed within swath profile polygons
	profile_dict['min'] = statsDict['mins'] # add minimum speed within swath profile polygons
	profile_dict['count'] = statsDict['counts'] # add count (num datapoints) within swath profile polygons
	profile_dict['std'] = statsDict['stds'] # add standard deviation of speed within swath profile polygons

	# Collect all transect data to write out
	timeNow = dt.datetime.now()
	profile_data={ 'sample_pts_lon_lat':sample_pts_lon_lat, 'sample_pts_PS':sample_pts_PS, 'sample_pts_frontdist':sample_pts_frontdist.tolist(),\
	'sampled_on_date_iso':timeNow.isoformat(), 'profile_shapefile':profile_shapefile, 'profile_shapefile_dir':profile_shapefile_dir,\
	'profile_shapefile_projection_Proj':originalSrsProj, 'profiles':profile_dict, 'featNum':featNum }

	# Write transect velocity data to json file
	jsonFnOut = 'tmpSwathVelocitySampling.json'
	wgt.writeJsonFile(profile_data,jsonFnOut)
	
	print "Output swath profile data to: " + jsonFnOut
	
except: # not sure why this error occuring
	print "Something failed while profiling transect " + str(featNum)
	print "For image " + rasterIn 
	




