'''
Making elevation profiles, glacier area, etc figures.
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

labCpuToggle = 1 # for ease in switching between laptop and lab cpu
billyCpuToggle = 0 # billy's laptop
dDist = 500 # interpolation distance along transect line
boxWidth = 1000 # swath profiler width across transect line

if billyCpuToggle == 1 and labCpuToggle == 1:
	print "Billy cpu and lab cpu toggle both on; have to choose one"
	exit()
	
if billyCpuToggle == 1:
	#shapefileIn='/Users/wiar9509/git/pySwath/stEliasTransectsUtm7n.shp'		
	#shapefileIn='/Users/wiar9509/git/pySwath/wrangellsTransectsUtm7n.shp'
	shapefileIn='/Users/wiar9509/Documents/generalScripts/swath/shapefileInputs/newWStETransects_19nov15_fromQGIS_EGPSG102006.shp'
	shapefileEPSG = 3338 # proj doesn't know EPSG 102006?
	#epsgNum = 32607 # spatial reference EPSG code
	#rasterFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/AK_BA_62_17/akv2ba_S8_062_017_016_2014_043_2014_059_hp_filt_3.0_vv.tif'
	#rasterDirectory = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/test/'
	rasterDirectory = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/'
	outDirectory = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/swath/'
	jsonDataFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/kaskWB_evenlySpace_100m_profiles_sampled_2015-12-04.json'
	wrs2Fn = '/Users/wiar9509/Documents/generalScripts/fromMarkFahnestock/wrs2_descending/wrs2_descending.shp'
elif labCpuToggle == 1:
	#shapefileIn='/Users/wiar9509/git/pySwath/stEliasTransectsUtm7n.shp'		
	#shapefileIn='/Users/wiar9509/git/pySwath/wrangellsTransectsUtm7n.shp'
	shapefileIn='/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/newWStETransectsUtm7n_24sep15.shp'
	wgs84centerlinesFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/new_centerlines_withChugach_06oct2016.shp'
	albersCenterlinesFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/new_centerlines_withChugach_epsg102006.shp'
	shapefileEPSG = 32607 # proj doesn't know EPSG 102006?
	rasterDirectory = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/'
	outDirectory = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/'
	#jsonDataFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/kaskWB_evenlySpace_100m_profiles_sampled_2015-12-04.json'
	wrs2Fn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/wrs2_descending/wrs2_descending.shp'
	pathRowImageFolder = '/Users/anderson/Desktop/ARMSTRONG/wrangells/imagery/L8/'
	wgsDemFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/dtm/alaska_yukon_merge_epsg4326.tif'
	albersDemFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/dtm/alaska_yukon_merge_epsg102006.tif'
	utmZoneShpFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/UTM_Zone_Boundaries/UTM_Zone_Boundaries.shp'
	albersRGIfn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/rgi_50_alaska_epsg102006.shp'
	rgiStudySubsetFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/rgi_subsetToStudyGlaciers_withChugach_epsg102006.shp'
## ELEVATION PROFILES
# elev profiles not yet working as of 07oct2016 wha

clines = fiona.open(albersCenterlinesFn)
shapefileEPSG = 102006
elevDict = {}
for featNum in range(0,len(clines)): # iterate over all transects
	# just hacking these for testing - will change to loop	
	#feat = lyr.GetFeature(59) # kaskWB if using the epsg102006 transects
	#feat = lyr.GetFeature(19) # kaskWB if using the st elias utm7n transects
	#feat = lyr.GetFeature(8) # nabeEB if using the wrangells utm7n transects
	#featNum = 17 # nabeEB if using new WStE transects
	feat = clines[featNum] # get feature
	#lineGeom = feat.geometry() # get geometry of profile line
	shapelyLineGeom = geometry.shape(feat['geometry'])
	lonMin,latMin,lonMax,latMax = shapelyLineGeom.bounds
	transName = feat['properties']['name'] # transect name
		
	# Make shapefile for swath profiling
	xInt,yInt = wgt.interpolateShapelyLineGeom(shapelyLineGeom,dDist) # interpolate vertices using shapely

	wgt.makePolylineShpFromCoords(xInt,yInt,shapefileEPSG,transName+'_swathSamplePoints.shp')
	slopes=wgt.calculateSlopes(xInt,yInt) # calculate "slope" (azimuth) of line segments
	perpSlopes=wgt.calculatePerpendicularSlopes(slopes) # calculate perpendicular "slope" (azimuth)
	allBoxVerts = wgt.makeAllBoxVertices(xInt,yInt,boxWidth) # get box coordinates for each interpolated point, with a specified box width (cross-transect coord)

	# Create the spatial reference
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(shapefileEPSG) # This hard-coded to UTM 7N w/ WGS84 datum. Would be nice to have this defined based of input shapefile

	wgt.makeMultipolygonFromBoxVertices(allBoxVerts,srs,outDirectory + transName+'_samplePolygons.shp') # Make shapefile for swath profiling

	#profileCoords = wgt.readLineCoordsFromMultiLineShp(shapefileIn,featNum) # get coordinates for PRE-INTERPOLATED centerline profile
	#geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']} # make shapely-style geometry dict
	#profGeom = geometry.shape(geomToPass) # make into a shapely geometry for testing intersection in wgt.getPathRowList()

	dist = wgt.distanceAlongLine(xInt,yInt) # get distance along line
	dist = np.delete(dist,0) # delete the first entry, which we don't have a value for b/c don't have a slope for the first entry (b/c need to points), so interpolator doesn't make a 0,0 point .. this quick and dirty

	interpProfileCoords = wgt.readVertexCoordsFromPolyline(transName+'_swathSamplePoints.shp') # read vertex coordinates from interpolated line

	stats = wgt.zonal_stats(outDirectory + transName+'_samplePolygons.shp', albersDemFn) # swath profile along transect
	zStats = wgt.readZonalStatsOutput(stats) # read swath profile output
	
	dz = np.diff(zStats['means'])
	dx = np.diff(dist)
	
	alpha = np.nan
	alpha = np.append(alpha,dz/dx)
		
	elevDict[transName] = {'dist':dist.tolist(),'slope':alpha.tolist(),'elevMean':zStats['means'],'elevMin':zStats['mins'],'elevMax':zStats['maxs'],'elevStd':zStats['stds'],'elevCount':zStats['counts']}
	
# Write transect velocity data to json file
jsonFnOut = os.path.split(albersCenterlinesFn)[1][:-4] + '_elevationProfiles.json'
wgt.writeJsonFile(elevDict,jsonFnOut)	


## PLOTTING ELEVATION PROFILES
color_idx = np.linspace(0, 1, len(elevDict))
cols = plt.cm.GnBu(color_idx)
transIter = elevDict.iterkeys()
meanMins = []
meanMeans = []
meanMaxs = []
meanDists = []


for i in range(0,len(elevDict)):
	transName = transIter.next()
	elevNow = elevDict[transName]
	randInd = np.random.randint(0,len(elevDict))
	plt.plot(np.array(elevNow['dist'])/1e3,elevNow['elevMean'],lw=2,alpha=0.6,c=cols[randInd])
	
	meanMins.append(elevNow['elevMean'][-1])
	meanMaxs.append(elevNow['elevMean'][0])
	meanMeans.append(np.nanmean(elevNow['elevMean']))
	meanDists.append(elevNow['dist'][-1])	

		
	
meanDist = np.mean(meanDists)/1e3		
medDist = np.median(meanDists)/1e3
meanMax = np.mean(meanMaxs)
meanMin = np.mean(meanMins)
meanMean = np.mean(meanMeans)
medMean = np.median(meanMeans)
medMin = np.median(meanMins)
medMax = np.median(meanMaxs)

xLims = (0,87)
yLims = (-50,4600)
plt.xlim(xLims)
plt.ylim(yLims)
plt.plot(np.array((medDist,medDist)),yLims,c='gray',ls='--')
plt.plot(xLims,(medMax,medMax),c='gray',ls='--')
plt.plot(xLims,(medMin,medMin),c='gray',ls='--')
plt.plot(xLims,(medMean,medMean),c='k',ls='--',lw=2)
plt.text(medDist+1,yLims[1]-300,'Median x$_{term}$ = ' + "%.0f" % medDist + ' km'  ,fontsize=14)
plt.text(xLims[1]-32,medMax+100,'Median z$_{max}$ = ' + "%.0f" % medMax + ' m'  ,fontsize=14)
plt.text(xLims[1]-33,medMean+100,'Median z$_{mean}$ = ' + "%.0f" % medMean + ' m'  ,fontsize=14)
plt.text(xLims[1]-30,medMin-300,'Median z$_{min}$ = ' + "%.0f" % medMin + ' m'  ,fontsize=14)
plt.ylabel('Elevation [m]',fontsize=18)
plt.xlabel('Distance from headwall [km]',fontsize=18)
plt.savefig('allTransects_elevationProfiles_medians.pdf')
plt.show()
plt.close()


## RGI comparison

with fiona.open(albersRGIfn) as rgi:
	areaList = []
	zminList = []
	zmaxList = []
	zmedList = []
	lmaxList = []

	for g in rgi:
		areaList.append(g['properties']['Area'])
		zminList.append(g['properties']['Zmin'])
		zmaxList.append(g['properties']['Zmax'])	
		zmedList.append(g['properties']['Zmed'])
		lmaxList.append(g['properties']['Lmax'])
		
with fiona.open(rgiStudySubsetFn) as rgiSub:	
	subAreaList = []
	subZminList = []
	subZmaxList = []
	subZmedList = []
	subLmaxList = []
	
	for g in rgiSub:
		subAreaList.append(g['properties']['Area'])
		subZminList.append(g['properties']['Zmin'])
		subZmaxList.append(g['properties']['Zmax'])	
		subZmedList.append(g['properties']['Zmed'])
		subLmaxList.append(g['properties']['Lmax'])	

plt.scatter(lmaxList,areaList,c='gray',alpha=0.6,lw=0,label='Alaska RGI 5.0')
plt.scatter(subLmaxList,subAreaList,s=50,c='c',alpha=0.8,label='Study glaciers')
plt.xlim((1e2,2.5e5))
plt.ylim((1e-1,1e4))
ax = plt.gca()
ax.set_yscale('log')	
ax.set_xscale('log')
ax.grid(True)
plt.ylabel('Area [m$^2$]',fontsize=16)
plt.xlabel('Length [m]',fontsize=16)
plt.legend(loc=2,frameon=False,numpoints=1,fontsize=14)	
plt.savefig('areaLength_vsRGI.pdf',dpi=300)
plt.show()
plt.close()

## Elevation histogram
rgiHist,edges = np.histogram(zmedList,range=(0,5e3),bins=1e2,density=1)
subHistr,subEdges = np.histogram(subZmedList,range=(0,5e3),bins=0.5e2,density=1)

fig, ax = plt.subplots()

rgiBars = ax.bar(edges[1:], rgiHist, 45, facecolor='gray',alpha=0.7,edgecolor=None,lw=0,label='Alaska RGI 5.0')
rgiBars = ax.bar(subEdges[1:], subHistr, 60, facecolor='c',alpha=0.8,edgecolor=None,lw=0,label='Study glaciers')
plt.xlim((0,4e3))
plt.ylim((0,1.1e-3))
plt.ylabel('Probability(z$_{med}$) [-]',fontsize=16)
plt.xlabel('Median elevation, z$_{med}$ [m]',fontsize=16)
plt.legend(loc=2,frameon=False,numpoints=1,fontsize=14)
plt.savefig('elevationHistogram_vsRGI.png',dpi=300)
plt.show()
plt.close()

#plt.scatter(lmaxList,zmedList,c='gray',alpha=0.6,lw=0,label='Alaska RGI 5.0')
#plt.scatter(subLmaxList,subZmedList,s=50,c='c',alpha=0.8,label='Study glaciers')

#plt.ylim((1e-1,1e4))
ax = plt.gca()
#ax.set_yscale('log')	
#ax.set_xscale('log')

#plt.savefig('areaLength_vsRGI.png',dpi=300)
plt.show()
plt.close()

	