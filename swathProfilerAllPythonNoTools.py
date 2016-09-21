'''
New script for swath profiler that does not require GRASS
begun 06 Sep 2016
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
# Switch directories just to get the wgt ``module''
# Don't get why need '%', but got from https://ipython.org/ipython-doc/2/interactive/reference.html
%cd /Users/wiar9509/git/wrangells/
import wrangellsGeospatialTools as wgt
%cd ~
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



dDist = 500 # interpolation distance along transect line
boxWidth = 1000 # swath profiler width across transect line
epsgNum = 32607 # spatial reference EPSG code
		
#shapefileIn='/Users/wiar9509/git/pySwath/stEliasTransectsUtm7n.shp'		
shapefileIn='/Users/wiar9509/git/pySwath/wrangellsTransectsUtm7n.shp'
rasterFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/AK_BA_62_17/akv2ba_S8_062_017_016_2014_043_2014_059_hp_filt_3.0_vv.tif'
rasterDirectory = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/test/'
jsonDataFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/kaskWB_evenlySpace_100m_profiles_sampled_2015-12-04.json'
wrs2Fn = '/Users/wiar9509/Documents/generalScripts/fromMarkFahnestock/wrs2_descending/wrs2_descending.shp'
				
## -------- ##
## ANALYSIS ##
## -------- ##


# Bringing in a test centerline geometry
#shapefileIn='/Users/wiar9509/Documents/generalScripts/swath/shapefileInputs/newWStETransects_19nov15_fromQGIS_EGPSG102006.shp'

shp = ogr.Open(shapefileIn)
lyr = shp.GetLayer(0)

if 0: # to get names and associated index
	for i in range(0,lyr.GetFeatureCount()):
		feat = lyr.GetFeature(i)
		print i, feat.GetField(1)
	
#feat = lyr.GetFeature(59) # kaskWB if using the epsg102006 transects
#feat = lyr.GetFeature(19) # kaskWB if using the st elias utm7n transects
feat = lyr.GetFeature(8) # nabeEB if using the wrangells utm7n transects
lineGeom = feat.geometry()
transName = feat.GetField(1) # transect name

# Make shapefile for swath profiling
xInt,yInt = wgt.interpolateVerticesOfPolyline(lineGeom,dDist) # interpolate vertices
slopes=wgt.calculateSlopes(xInt,yInt) # calculate "slope" (azimuth) of line segments
perpSlopes=wgt.calculatePerpendicularSlopes(slopes) # calculate perpendicular "slope" (azimuth)
allBoxVerts = wgt.makeAllBoxVertices(xInt,yInt,boxWidth) # get box coordinates for each interpolated point, with a specified box width (cross-transect coord)

# Create the spatial reference
srs = osr.SpatialReference()
srs.ImportFromEPSG(epsgNum) # This hard-coded to UTM 7N w/ WGS84 datum. Would be nice to have this defined based of input shapefile

wgt.makeMultipolygonFromBoxVertices(allBoxVerts,srs,transName+'_samplePolygons.shp') # Make shapefile for swath profiling


profileCoords = wgt.readLineCoordsFromMultiLineShp(shapefileIn,14)
geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']}
profGeom = geometry.shape(geomToPass)


# now run through the WRS2 outlines and find path,rows that intersect with the MultiPoint object
# get list of path/rows overlapping with profile geometry
pr_list = wgt.getPathRowList(profGeom,wrs2Fn)
# get list of images in folder
all_frames = wgt.listImagesToSample(pr_list,rasterDirectory)

profiles=[]
for frame in all_frames:
	path=frame['path']
	print path
	row=frame['row']
	for inf in frame['infiles']:
		

		
		raster = gdal.Open(inf)
		xyTfm = wgt.createTransformation(inf) # transformation parameters for converting raster pixels to ground coordinates

		stats=wgt.zonal_stats(transName+'_samplePolygons.shp', inf, nodata_value=-1)
		statsDict = wgt.readZonalStatsOutput(stats)
		dist = wgt.distanceAlongLine(xInt,yInt)
		dist = np.delete(dist,0) # delete the first entry, which we don't have a value for b/c don't have a slope for the first entry (b/c need to points), so interpolator doesn't make a 0,0 point .. this quick and dirty
		
				

		profile_dict = wgt.getCorrelationInfoFromFn(inf)
		profile_dict['dist'] = dist
		profiles.append(profile_dict)






# Brute force, just to every file in directory
fileList = glob.glob(rasterDirectory+'*064017*')
for rasterFn in fileList:
	print "Profiling file: " + fileList[0].split('/')[-1]
	raster = gdal.Open(rasterFn)
	xyTfm = createTransformation(rasterFn) # transformation parameters for converting raster pixels to ground coordinates

	stats=zonal_stats('test.shp', rasterFn, nodata_value=-1)
	statsDict = readZonalStatsOutput(stats)
	dist = distanceAlongLine(xInt,yInt)
	dist = np.delete(dist,0) # delete the first entry, which we don't have a value for b/c don't have a slope for the first entry (b/c need to points), so interpolator doesn't make a 0,0 point .. this quick and dirty
	statsDict['dist'] = dist

	plt.plot(statsDict['dist'],statsDict['means'],lw=2)
	plt.plot(statsDict['dist'],statsDict['counts'],c='gray',ls='--')

plt.show()
plt.close()






		
		
	


