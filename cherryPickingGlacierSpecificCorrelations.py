'''
Glacier specific correlations
09 oct 2016
Cherry picking as intensely as possible
'''

## -------------- ##
## IMPORT MODULES ##
## -------------- ##

import os
import sys
import numpy as np
# Switch directories just to get the wgt ``module''
# Don't get why need '%', but got from https://ipython.org/ipython-doc/2/interactive/reference.html
#%cd /Users/wiar9509/git/wrangells/ # billy cpu
#%cd /Users/anderson/Desktop/ARMSTRONG/git/wrangells # lab cpu, added to python path on lab cpu so not needed
import wrangellsGeospatialTools as wgt 
#%cd ~
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

#shapefileIn='/Users/wiar9509/git/pySwath/stEliasTransectsUtm7n.shp'		
#shapefileIn='/Users/wiar9509/git/pySwath/wrangellsTransectsUtm7n.shp'
#shapefileIn='/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/newWStETransectsUtm7n_24sep15.shp'
wgs84centerlinesFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/new_centerlines_withChugach_06oct2016.shp'
albersCenterlinesFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/new_centerlines_withChugach_epsg102006.shp'
#shapefileIn='/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/centerline_transects_onlyUtm6n_08oct2016.shp'

#	shapefileIn='/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/centerline_transects_onlyUtm7n_08oct2016.shp'
#	shapefileEPSG = 32607 # proj doesn't know EPSG 102006?	
rasterDirectory = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/'
outDirectory = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/'
#jsonDataFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/kaskWB_evenlySpace_100m_profiles_sampled_2015-12-04.json'
wrs2Fn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/wrs2_descending/wrs2_descending.shp'
pathRowImageFolder = '/Users/anderson/Desktop/ARMSTRONG/wrangells/imagery/L8/'
wgsDemFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/dtm/alaska_yukon_merge_epsg4326.tif'
albersDemFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/dtm/alaska_yukon_merge_epsg102006.tif'
utmZoneShpFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/UTM_Zone_Boundaries/UTM_Zone_Boundaries.shp'
glacierSpecificCorrListFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/selectSummerCorr_glacierSpecific_064_017.csv'
utm7onlyCenterlinesFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/centerline_transects_onlyUtm7n_08oct2016.shp'
shapefileEPSG = 32607 # proj doesn't know EPSG 102006?

## -------- ##
## ANALYSIS ##
## -------- ##

shapefileIn = utm7onlyCenterlinesFn

# List of select correlations, cherry picking for specific glaciers
corrList = np.genfromtxt(glacierSpecificCorrListFn,delimiter=',',skip_header=1,dtype=None)
uniqueNames = np.unique(corrList[:,2]) # just get transect names

shp = fiona.open(utm7onlyCenterlinesFn)

for transName in uniqueNames:

	# iterate to find current centerline 
	for i in shp:
		if i['properties']['name'] == transName:
			fionaCl = i # get fiona record for current centerline profile
			
	# Make shapefile for swath profiling
	shapelyLineGeom = geometry.shape(fionaCl['geometry'])
	xInt,yInt = wgt.interpolateShapelyLineGeom(shapelyLineGeom,dDist) # interpolate vertices using shapely

	wgt.makePolylineShpFromCoords(xInt,yInt,shapefileEPSG,transName+'_swathSamplePoints.shp')
	slopes=wgt.calculateSlopes(xInt,yInt) # calculate "slope" (azimuth) of line segments
	perpSlopes=wgt.calculatePerpendicularSlopes(slopes) # calculate perpendicular "slope" (azimuth)
	allBoxVerts = wgt.makeAllBoxVertices(xInt,yInt,boxWidth) # get box coordinates for each interpolated point, with a specified box width (cross-transect coord)

	# Create the spatial reference
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(shapefileEPSG) # This hard-coded to UTM 7N w/ WGS84 datum. Would be nice to have this defined based of input shapefile

	wgt.makeMultipolygonFromBoxVertices(allBoxVerts,srs,outDirectory + transName+'_samplePolygons.shp') # Make shapefile for swath profiling

#	profileCoords = wgt.readLineCoordsFromMultiLineShp(shapefileIn,featNum) # get coordinates for PRE-INTERPOLATED centerline profile
#	geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']} # make shapely-style geometry dict
#	profGeom = geometry.shape(geomToPass) # make into a shapely geometry for testing intersection in wgt.getPathRowList()

	dist = wgt.distanceAlongLine(xInt,yInt) # get distance along line
	dist = np.delete(dist,0) # delete the first entry, which we don't have a value for b/c don't have a slope for the first entry (b/c need to points), so interpolator doesn't make a 0,0 point .. this quick and dirty

	interpProfileCoords = wgt.readVertexCoordsFromPolyline(transName+'_swathSamplePoints.shp') # read vertex coordinates from interpolated line

	# Rename spatial variables for writing to profile_data dict
	sampled_on_date_iso = dt.datetime.now().isoformat()
	sample_pts_lon_lat = interpProfileCoords['samplePtsLatLon']
	sample_pts_frontdist = dist
	sample_pts_PS = interpProfileCoords['samplePtsXy']
	profile_shapefile_dir, trash , profile_shapefile = shapefileIn.rpartition('/')
	profile_feature = transName 
	profile_shapefile_projection_Proj = srs.ExportToProj4()

	# Get correlation list for this centerline
	clInd = np.where(corrList[:,2]==clNow)
	corrSubList = corrList[clInd,0]			

	profiles=[] # initialize
	for frame in corrSubList[0]: # iterate over all specified images
		
		# Get information about correlation
		corrDict = wgt.getCorrelationInfoFromFn(frame)			
		raster = gdal.Open(frame) # open one correlation
		xyTfm = wgt.createTransformation(frame) # transformation parameters for converting raster pixels to ground coordinates

		path=corrDict['path']
		row=corrDict['row']
		
		stats=wgt.zonal_stats(outDirectory + transName+'_samplePolygons.shp', frame, nodata_value=-1) # swath profile along transect
		statsDict = wgt.readZonalStatsOutput(stats)
		
		profile_dict = corrDict # strip filename for correlation parameters
		profile_dict['speed'] = statsDict['means'] # add mean speed in swath profile to dict. keep 'speed' for consistency with Mark's old routine
		profile_dict['speedMean'] = statsDict['means'] # add mean speed (duplicate, but just to be obvious in future
		profile_dict['speedMax'] = statsDict['maxs'] # add maximum speed within swath profile polygons
		profile_dict['speedMin'] = statsDict['mins'] # add minimum speed within swath profile polygons
		profile_dict['speedCount'] = statsDict['counts'] # add count (num datapoints) within swath profile polygons
		profile_dict['speedStd'] = statsDict['stds'] # add standard deviation of speed within swath profile polygons
		#profile_dict['image1_snowIndex'] = ndsi1mean # add mean snow index profile for image 1
		#profile_dict['image2_snowIndex'] = ndsi2mean # add mean snow index profile for image 1
			
		profiles.append(profile_dict) # add to list

	# Translate timeslices across space (profiles) to timeseries across time at a point (pts)
	pts = wgt.convertProfilesListToPointTimeseries(profiles)

	# Collect all transect data to write out
	timeNow = dt.datetime.now()
	profile_data={ 'sample_pts_lon_lat':sample_pts_lon_lat, 'sample_pts_PS':sample_pts_PS, 'sample_pts_frontdist':sample_pts_frontdist.tolist(),\
					'sampled_on_date_iso':timeNow.isoformat(), 'profile_shapefile':profile_shapefile, 'profile_shapefile_dir':profile_shapefile_dir,\
					'profile_shapefile_projection_Proj':profile_shapefile_projection_Proj, 'profiles':profiles, 'pts':pts, 'transect_name':transName }

	# Write transect velocity data to json file
	jsonFnOut = transName + '_cherryPicking_swathVelocitySampling_' + timeNow.strftime("%Y-%m-%d") + '.json'
	wgt.writeJsonFile(profile_data,jsonFnOut)


	# Plot raw data

	#transName = 'tanaWB'
	ending = '_cherryPicking_swathVelocitySampling_2016-10-09.json'
	profile_data = wgt.readJsonFile(transName + ending)
	nProfiles = len(profile_data['profiles'])
	dist = profile_data['sample_pts_frontdist']
	stdMax = 1
	countMin = 10

	for i in range(0,nProfiles):
		u = profile_data['profiles'][i]['speed']
		counts = profile_data['profiles'][i]['speedCount']
		stds = profile_data['profiles'][i]['speedStd']	
	
		badDataInd = np.logical_or(np.array(counts)<countMin, np.array(stds)>stdMax)
		uMa = np.copy(u)
		uMa[badDataInd] = np.nan
		#uMa = np.ma.masked_where(np.logical_or(np.array(counts)<countMin, np.array(stds)>stdMax),u)
	
		plt.subplot(4,1,1)
		plt.plot(dist,u,lw=2)
		plt.ylabel('vel')
		plt.title(transName)
		plt.xlim((np.min(dist),np.max(dist)))	
		plt.subplot(4,1,2)	
		plt.plot(dist,counts,lw=2)
		plt.plot((dist[0],dist[-1]),(countMin,countMin),c='k',ls='--',lw=1)			
		plt.ylabel('counts')
		plt.xlim((np.min(dist),np.max(dist)))		
		plt.subplot(4,1,3)	
		plt.plot(dist,stds,lw=2)
		plt.plot((dist[0],dist[-1]),(stdMax,stdMax),c='k',ls='--',lw=1)		
		plt.ylabel('std')		
		plt.xlim((np.min(dist),np.max(dist)))
		plt.subplot(4,1,4)	
		plt.plot(dist,uMa,lw=2)
		plt.xlim((np.min(dist),np.max(dist)))
		plt.ylabel('uMask')
		plt.xlabel('Downglacier distance')
	
	plt.savefig(transName + 'filteringVelocity_stdCounts.pdf')
	#plt.show()
	plt.close()

	
	







