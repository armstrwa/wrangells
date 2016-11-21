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
import matplotlib.gridspec as gridspec

## --------- ##
## FUNCTIONS ##
## --------- ##

# functions loaded from wrangellsGeospatialTools.pyc


def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-0):-(window_len/2)]


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
glacierSpecificCorrListFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/selectCorrForSpecificGlaciersSummerAndWinter_064_017.csv'
utm7onlyCenterlinesFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/centerline_transects_onlyUtm7n_08oct2016.shp'
shapefileEPSG = 32607 # proj doesn't know EPSG 102006?
elevationProfileFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/new_centerlines_withChugach_epsg102006_elevationProfiles.json'

ndsiFn1 = '/Users/anderson/Desktop/ARMSTRONG/wrangells/imagery/L8/64_17/LC80640172013214LGN00/LC80640172013214LGN00_ndsi.TIF'
ndsiFn2 = '/Users/anderson/Desktop/ARMSTRONG/wrangells/imagery/L8/64_17/LC80640172014233LGN00/LC80640172014233LGN00_ndsi.TIF'

## -------- ##
## ANALYSIS ##
## -------- ##

velocityChangeDict = {}

shapefileIn = utm7onlyCenterlinesFn

# Load elevation data
elevDict=wgt.readJsonFile(elevationProfileFn)

# List of select correlations, cherry picking for specific glaciers
corrList = np.genfromtxt(glacierSpecificCorrListFn,delimiter=',',skip_header=1,dtype=None)
uniqueNames = np.unique(corrList[:,2]) # just get transect names

shp = fiona.open(utm7onlyCenterlinesFn)

for transName in uniqueNames:
	if transName == '': # empty transname for some reason
		continue
	elif transName in ('chitSB','wForCL'): # no winter correlations for chitSB or wForCl
		print "No winter correlations: ", transName
		continue
	
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
	clInd = np.where(corrList[:,2]==transName)
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


	
	# get ndsi
	stats=wgt.zonal_stats(outDirectory + transName+'_samplePolygons.shp',ndsiFn1)
	ndsi1=wgt.readZonalStatsOutput(stats)
	stats=wgt.zonal_stats(outDirectory + transName+'_samplePolygons.shp',ndsiFn2)
	ndsi2=wgt.readZonalStatsOutput(stats)
	ndsiToget = np.vstack((np.array(ndsi1['means']),np.array(ndsi2['means'])))
	ndsiMean = np.nanmean(ndsiToget,axis=0)
	
 
	
	# read elevation dict
	elevData = elevDict[transName]

	# Plot raw data
	ending = '_cherryPicking_swathVelocitySampling_2016-10-10.json'
	profile_data = wgt.readJsonFile(outDirectory + transName + ending)
	nProfiles = len(profile_data['profiles'])
	dist = profile_data['sample_pts_frontdist']
	stdMax = 1
	countMin = 10

	# Classify as summer/winter correlations and clean bad data
	winterCorr = [] # intializing
	summerCorr = []
	winterSpeeds = np.nan * np.ones((len(profile_data['profiles'][0]['speed']))) # intializing arrays for stats
	summerSpeeds = np.nan * np.ones((len(profile_data['profiles'][0]['speed'])))
	
	for i in range(0,nProfiles):
		profileNow = profile_data['profiles'][i]
		u = profileNow['speed']
		counts = profileNow['speedCount']
		stds = profileNow['speedStd']
		midYearFrac = profileNow['mid_dec_year']
		
		midDoy =  wgt.yearFracToDoy(midYearFrac)
		
		badDataInd = np.logical_or(np.array(counts)<countMin, np.array(stds)>stdMax)
		uMa = np.copy(u)
		uMa[badDataInd] = np.nan
		
		profileNow['stdMax']=stdMax
		profileNow['countMin']=countMin
		profileNow['speedMeanClean']=uMa
				
		if midDoy > 250 or midDoy < 120:
			winterCorr.append(profileNow)
			winterSpeeds = np.vstack((winterSpeeds,uMa))
		else:
			summerCorr.append(profileNow)
			summerSpeeds = np.vstack((summerSpeeds,uMa))	
	
	# Seasonal statistics
	w25,w33,w50,w66,w75 = np.nanpercentile(winterSpeeds,(25,33,50,66,75),axis=0)
	s25,s50,s75 = np.nanpercentile(summerSpeeds,(25,50,75),axis=0)
	winterDataNum = np.sum(np.isfinite(winterSpeeds),axis=0)	
	summerDataNum = np.sum(np.isfinite(summerSpeeds),axis=0)	

	# Flipping distance axis	
	termDist = flipDistance(dist)
	normTermDist = normalizeDistance(termDist)
	# Velocity change
	diff50 = s50-w50
	diff7525 = s75 - w25
	diff2575 = s25 - w75	
	
	velocityChangeDict[transName] = {'distTerm':termDist.tolist(),'distTermNorm':normTermDist.tolist(),'diffMedians':diff50.tolist(),'minDiff':diff2575.tolist(),'maxDiff':diff7525.tolist()}
	
	# Plotting seasonal velocity
	winterCind = np.linspace(0.25, 1, len(winterCorr))
	winterCols = plt.cm.GnBu(winterCind)

	summerCind = np.linspace(0, 1, len(summerCorr))
	summerCols = plt.cm.OrRd(summerCind)
	
	alphaVal = 0.6

	if 0: # plot counts, velocoty, elevation, ndsi	
		gs = gridspec.GridSpec(4,1,height_ratios=[1,3,1,1])
		ax0=plt.subplot(gs[3])
		ax1=plt.subplot(gs[0])
		ax2=plt.subplot(gs[2])
		ax2b = ax2.twinx()	
		ax3=plt.subplot(gs[1])
	
		for i in range(0,len(winterCorr)):
			profileNow = winterCorr[i]
			u = profileNow['speed']
			counts = profileNow['speedCount']
			stds = profileNow['speedStd']	
			uClean = profileNow['speedMeanClean']
			countMin = profileNow['countMin']
			stdMax = profileNow['stdMax']
			randInd = np.random.randint(0,len(winterCorr))
		
	#		ax0.plot(np.array(dist)/1e3,u,lw=2,c=winterCols[randInd],alpha=alphaVal)
			#ax1.plot(np.array(dist)/1e3,counts,lw=2,c=winterCols[randInd],alpha=alphaVal)
			#ax1.plot((np.array(dist)[0]/1e3,np.array(dist)[-1]/1e3),(countMin,countMin),c='k',ls='--',lw=1)			
	#		ax2.plot(np.array(dist)/1e3,stds,lw=2,c=winterCols[randInd],alpha=alphaVal)
	#		ax2.plot((np.array(dist)[0]/1e3,np.array(dist)[-1]/1e3),(stdMax,stdMax),c='k',ls='--',lw=1)		
			ax3.scatter(np.array(dist)/1e3,uClean,lw=0,c=winterCols[randInd],alpha=alphaVal*0.6)
		
		for i in range(0,len(summerCorr)):
			profileNow = summerCorr[i]
			u = profileNow['speed']
			counts = profileNow['speedCount']
			stds = profileNow['speedStd']	
			uClean = profileNow['speedMeanClean']
			countMin = profileNow['countMin']
			stdMax = profileNow['stdMax']
			randInd = np.random.randint(0,len(summerCorr))
		
	#		ax0.plot(np.array(dist)/1e3,u,lw=2,c=summerCols[randInd],alpha=alphaVal)	
			#ax1.plot((np.array(dist)[0]/1e3,np.array(dist)[-1]/1e3),(countMin,countMin),c='k',ls='--',lw=1)			
	#		ax2.plot(np.array(dist)/1e3,stds,lw=2,c=summerCols[randInd],alpha=alphaVal)
	#		ax2.plot((np.array(dist)[0]/1e3,np.array(dist)[-1]/1e3),(stdMax,stdMax),c='k',ls='--',lw=1)		
			ax3.scatter(np.array(dist)/1e3,uClean,lw=0,c=summerCols[randInd],alpha=alphaVal*0.6)
	

		# plot seasonal statistics
		ax0.plot(np.array(dist)/1e3,ndsi1['means'],lw=1,c='palegreen')
		ax0.plot(np.array(dist)/1e3,ndsi2['means'],lw=1,c='mediumspringgreen')	
		ax0.plot(np.array(dist)/1e3,ndsiMean,lw=2,c='dodgerblue')		
		ax0.plot((np.array(dist)[0]/1e3,np.array(dist)[-1]/1e3),(0.6,0.6),lw=1,ls='--',c='k')		
		ax1.plot(np.array(dist)/1e3,winterDataNum,lw=2,c='cornflowerblue')
		ax1.plot(np.array(dist)/1e3,summerDataNum,lw=2,c='darkred')
		ax2b.scatter(np.array(elevData['dist'])/1e3,elevData['slope'],lw=2,color='lightgray')
		ax2b.plot((np.array(dist)[0]/1e3,np.array(dist)[-1]/1e3),(0,0),ls='--',lw=1,color='gray')
		ax2.plot(np.array(elevData['dist'])/1e3,elevData['elevMean'],lw=2,color='k')
	#	ax3.plot(np.array(dist)/1e3,w25,c='cornflowerblue',lw=1)	
		ax3.plot(np.array(dist)/1e3,w50,c='cornflowerblue',lw=2)
	#	ax3.plot(np.array(dist)/1e3,w75,c='cornflowerblue',lw=1)	
		ax3.fill_between(np.array(dist)/1e3,w25,w75,facecolor='powderblue',color='cornflowerblue',alpha=0.7)		
		ax3.plot(np.array(dist)/1e3,s50,c='darkred',lw=2)
		ax3.fill_between(np.array(dist)/1e3,s25,s75,facecolor='salmon',color='darkred',alpha=0.7)		
	#	ax3.plot(np.array(dist)/1e3,s75,c='darkred',lw=1)
		
		ax1.set_title(transName,fontsize=20)
		ax0.set_ylabel('NDSI [-]')
		ax0.set_xlim((np.array(dist)[0]/1e3,np.array(dist)[-1]/1e3))
		ax1.set_ylabel('Counts [-]')
		ax1.set_xlim((np.array(dist)[0]/1e3,np.array(dist)[-1]/1e3))	
		ax2.set_ylabel('Elevation [m]')
		ax2b.set_ylabel('Slope [m/m]')	
		ax2.set_xlim((np.array(dist)[0]/1e3,np.array(dist)[-1]/1e3))	
		ax3.set_ylabel('Vel. [m/d]')	
		maxSpeed = np.nanmax((np.nanmax(summerSpeeds),np.nanmax(winterSpeeds)))
		ax0.set_xlabel('Distance from headwall [km]')
		ax3.set_ylim((-0.05,np.nanmin((2.5,maxSpeed))))
		ax3.set_xlim((np.array(dist)[0]/1e3,np.array(dist)[-1]/1e3))	
		plt.savefig(transName + '_velocityElevationNdsiProfiles.pdf')
		#plt.show()
		plt.close()
	
	if 1: # plot counts and velocity
		gs = gridspec.GridSpec(2,1,height_ratios=[1,4])
				
		ax0=plt.subplot(gs[0])
		ax1=plt.subplot(gs[1])
		#ax2=plt.subplot(gs[2])		

		for i in range(0,len(winterCorr)):
			profileNow = winterCorr[i]
			u = profileNow['speed']
			counts = profileNow['speedCount']
			stds = profileNow['speedStd']	
			uClean = profileNow['speedMeanClean']
			countMin = profileNow['countMin']
			stdMax = profileNow['stdMax']
			randInd = np.random.randint(0,len(winterCorr))
		
			ax1.scatter(np.array(termDist)/1e3,uClean,lw=0,c=winterCols[randInd],alpha=alphaVal)
		
		
		for i in range(0,len(summerCorr)):
			profileNow = summerCorr[i]
			u = profileNow['speed']
			counts = profileNow['speedCount']
			stds = profileNow['speedStd']	
			uClean = profileNow['speedMeanClean']
			countMin = profileNow['countMin']
			stdMax = profileNow['stdMax']
			randInd = np.random.randint(0,len(summerCorr))
			
			ax1.scatter(np.array(termDist)/1e3,uClean,lw=0,c=summerCols[randInd],alpha=alphaVal)

		ax0.plot(np.array(termDist)/1e3,winterDataNum,lw=2,c='cornflowerblue')
		ax0.plot(np.array(termDist)/1e3,summerDataNum,lw=2,c='darkred')			
		ax0.set_ylabel('Counts [-]')
		ax0.set_title(transName,fontsize=20)
		ax0.set_xlim((np.array(termDist)[0]/1e3,np.array(termDist)[-1]/1e3))
		ax0.set_ylim((-1,15))
		ax1.plot(np.array(termDist)/1e3,w50,c='cornflowerblue',lw=2)	
		ax1.fill_between(np.array(termDist)/1e3,w25,w75,facecolor='powderblue',color='cornflowerblue',alpha=0.7)		
		ax1.plot(np.array(termDist)/1e3,s50,c='darkred',lw=2)
		ax1.fill_between(np.array(termDist)/1e3,s25,s75,facecolor='salmon',color='darkred',alpha=0.7)					
		ax1.set_xlim((np.array(termDist)[0]/1e3,np.array(termDist)[-1]/1e3))
		maxSpeed = np.nanmax((np.nanmax(summerSpeeds),np.nanmax(winterSpeeds)))
		ax1.set_xlabel('Distance from terminus [km]')
		ax1.set_ylabel('Velocity [m d$^{-1}$]')		
		ax1.set_ylim((-0.05,np.nanmin((1.75,maxSpeed))))
		plt.savefig(transName + '_justCountsAndVelocity.pdf')
		plt.close()
						
	if 0: # plot velocity difference, ndsi, elev/slope
		gs = gridspec.GridSpec(3,1,height_ratios=[1,3,1])
				
		ax0=plt.subplot(gs[0])
		ax0b=ax0.twinx()
		ax1=plt.subplot(gs[1])
		ax2=plt.subplot(gs[2])		
		

		
		ax0b.scatter(np.array(flipDistance(elevData['dist']))/1e3,elevData['slope'],lw=2,color='lightgray')
		ax0b.plot((np.array(termDist)[0]/1e3,np.array(termDist)[-1]/1e3),(0,0),ls='--',lw=1,color='gray')
		ax0.plot(np.array(flipDistance(elevData['dist']))/1e3,elevData['elevMean'],lw=2,color='k')
		ax1.plot(np.array(termDist)/1e3,diff50,c='darkred',lw=2)
		ax1.plot((np.array(termDist[0])/1e3,np.array(termDist[-1])/1e3),(0,0),c='k',lw=1,ls='--')	
		ax1.fill_between(np.array(termDist)/1e3,diff2575,diff7525,facecolor='powderblue',color='cornflowerblue',alpha=0.7)
		ax2.plot(np.array(termDist)/1e3,ndsi1['means'],lw=1,c='palegreen')
		ax2.plot(np.array(termDist)/1e3,ndsi2['means'],lw=1,c='mediumspringgreen')	
		ax2.plot(np.array(termDist)/1e3,ndsiMean,lw=2,c='dodgerblue')
		ax1.set_xlim((np.min(termDist)/1e3,np.max(termDist)/1e3))
		ax0.set_xlim((np.min(termDist)/1e3,np.max(termDist)/1e3))		
		ax2.set_xlim((np.min(termDist)/1e3,np.max(termDist)/1e3))	
		ax2.plot((np.array(termDist)[0]/1e3,np.array(termDist)[-1]/1e3),(0.6,0.6),lw=1,ls='--',c='k')			
		ax1.set_ylim((-0.4,0.6))	
		ax2.set_ylim((0,0.8))
		ax2.set_xlabel('Distance from terminus [km]',fontsize=20)
		ax2.set_ylabel('NDSI [-]',fontsize=16)						
		ax1.set_ylabel('Speedup [m d$^{-1}$]',fontsize=16)								
		ax0.set_ylabel('Elevation [m]',fontsize=16)										
		ax0b.set_ylabel('Slope [m/m]',fontsize=16)
		ax0.set_title(transName,fontsize=20)												
		ax1.invert_xaxis()
		ax0.invert_xaxis()
		ax2.invert_xaxis()				
		#plt.show()
		plt.savefig(transName + '_velocityDifferenceNDSIandElevation.pdf')		
		plt.close()
		
		
		
### PLOTTING DIFFERENCE FOR ALL GLACIERS


speedupJsonFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/speedups.json'
speedDict = wgt.readJsonFile(speedupJsonFn)

goodTransects = ('barnEB','tanaWB','rootCL','kennCL','nabeWB','nabeEB','gateCL','hawkCL','chitNB','chisEB')
transNums = np.linspace(0.25, 1, len(goodTransects))
transCols = plt.cm.Accent(transNums)
colIter = 0
for transName in speedDict:
	if transName in goodTransects:
		transNow = speedDict[transName]
		termDist = transNow['distTerm']
		diff50 = transNow['diffMedians']
		minDiff = transNow['minDiff']
		maxDiff = transNow['maxDiff']		
		#intFun = interp1d(termDist,diff50,kind='quadratic')
		#xnew = np.linspace(np.min(termDist),np.max(termDist),300)
		#intFun = si.UnivariateSpline(termDist,diff50,k=5)
		diff50int = smooth(np.array(diff50))
		
		plt.scatter(np.array(termDist)/1e3,diff50,alpha=0.2,lw=0,color=transCols[colIter])
		plt.plot(np.array(termDist)/1e3,diff50int,lw=2.5,alpha=0.9,label=transName,c=transCols[colIter])		
		plt.fill_between(np.array(termDist)/1e3,minDiff,maxDiff,color=transCols[colIter],alpha=0.1,lw=0)
		colIter+=1

plt.plot((0,80),(0,0),lw=2,c='k',ls='--')
plt.ylabel('Speedup [m d$^{-1}$]',fontsize=18)
plt.xlabel('Distance from terminus [km]',fontsize=18)
plt.legend(loc=3,frameon=False,ncol=2)
plt.ylim((-.5,.5))
plt.xlim((0,80))
plt.gca().invert_xaxis()
plt.savefig('speedupSquiggles_withFill.pdf')
plt.show()
plt.close()

# normalized distance
colIter = 0
for transName in speedDict:
	if transName in goodTransects:
		transNow = speedDict[transName]
		termDist = transNow['distTermNorm']
		diff50 = transNow['diffMedians']
		minDiff = transNow['minDiff']
		maxDiff = transNow['maxDiff']
				
		#intFun = interp1d(termDist,diff50,kind='quadratic')
		#xnew = np.linspace(np.min(termDist),np.max(termDist),300)
		#intFun = si.UnivariateSpline(termDist,diff50,k=5)
		diff50int = smooth(np.array(diff50))
		
		plt.scatter(np.array(termDist),diff50,alpha=0.2,lw=0,color=transCols[colIter])
		plt.plot(np.array(termDist),diff50int,lw=2,alpha=0.9,label=transName,c=transCols[colIter])		
		plt.fill_between(np.array(termDist),minDiff,maxDiff,color=transCols[colIter],alpha=0.1,lw=0)
		colIter+=1

plt.plot((0,1),(0,0),lw=2,c='k',ls='--')
plt.ylabel('Speedup [m d$^{-1}$]',fontsize=18)
plt.xlabel('Normalized distance from terminus [-]',fontsize=18)
plt.legend(loc=3,frameon=False,ncol=2)
plt.ylim((-.5,.5))
plt.xlim((0,1))
plt.gca().invert_xaxis()
plt.savefig('speedupSquiggles_normalizedDist_withFill.pdf')
plt.show()
plt.close()

		
