'''
Filtering figure
07 dec 2016
'''

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
import matplotlib.gridspec as gridspec



import json
import matplotlib.pyplot as plt
import numpy as np
import wrangellsGeospatialTools as wgt 
import matplotlib.gridspec as gridspec
import glob


## USER DEFINED
#jsonFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/kennCL_swathVelocitySampling_2016-10-03.json'
#jsonFolder =  '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/'
jsonFolder = '/Users/wiar9509/Google Drive/wrangells/json/'
countMin = 10 # minimum pixels within swath profiler box to be considered good data
stdMax = 1 # maximum standard devuation of velocity within swath profiler box to be considered good data
summerStartEnd=((150,250)) # start and end DOYs to define the summer season for seasonal velocity split
elevationProfileFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/new_centerlines_withChugach_epsg102006_elevationProfiles.json'

## PROCESSING
## Initializing
velocityChangeDict = {}
velocityChangeDict2 = {}
jsonList = glob.glob(jsonFolder + '*2016-10-05.json')

jsonFn = jsonList[0]

profile_data = wgt.readJsonFile(jsonFn) # read in json velocity data
transName = profile_data['transect_name'] # get transect name	
print "Processing " + transName
nProfiles = len(profile_data['profiles']) # number of profiles (time slices) within json data
dist = profile_data['sample_pts_frontdist'] # distance from headwall [m]
termDist = wgt.flipDistance(dist) # Flipping distance axis, so zero at terminus
normTermDist = wgt.normalizeDistance(termDist) # Normalizing distance by glacier length
profile_data  = wgt.swathFilter(profile_data,countMin,stdMax) # filter velocity results based on counts and standard deviation within swath profiler boxes


# Temporal filtering
uMatxRaw =wgt.makeVelocityMatrix(profile_data,'speed') # make velocity matrix (nCols = len(space), nRows = len(time) )
uMatxSwath = wgt.makeVelocityMatrix(profile_data,'speedMeanClean') # make velocity matrix (nCols = len(space), nRows = len(time) )
uMatxRawOut = wgt.removeOutliers(uMatxRaw) # remove outliers
uMatxSwathOut = wgt.removeOutliers(uMatxSwath) # remove outliers
uMatxSwathOut2 = wgt.removeOutliers(uMatxSwathOut) # second iteration removing outliers
# All 2 suffixes below indicate 2 iterations of outerlier removal

# Put temporally filtered data back into profile_data dict form
profileDataSwathOutliers = wgt.putTempFiltDataBack(profile_data,uMatxSwathOut)
profileDataSwathOutliers2 = wgt.putTempFiltDataBack(profile_data,uMatxSwathOut2)

# Calculate temporal mean and standard deviation of velocity at each spatial location
meanRaw = np.nanmean(uMatxRaw,axis=0)
stdRaw = np.nanstd(uMatxRaw,axis=0)
meanSwath = np.nanmean(uMatxSwath,axis=0)
stdSwath = np.nanstd(uMatxSwath,axis=0)
meanTemp = np.nanmean(uMatxSwathOut,axis=0)
stdTemp = np.nanstd(uMatxSwathOut,axis=0)

# Break velocities into summer and winter categories
seasonalVelocityDict = wgt.splitSummerWinter(profileDataSwathOutliers,summerStartEnd)
seasonalVelocityDict2 = wgt.splitSummerWinter(profileDataSwathOutliers2,summerStartEnd)
speedupDictNow = wgt.calculateSeasonalSpeedup(seasonalVelocityDict,termDist)
speedupDictNow2 = wgt.calculateSeasonalSpeedup(seasonalVelocityDict2,termDist)
velocityChangeDict[transName] = speedupDictNow
velocityChangeDict2[transName] = speedupDictNow2



nProfiles = len(profileDataSwathOutliers['profiles'])
dist = profileDataSwathOutliers['sample_pts_frontdist']
termDist = wgt.flipDistance(dist)

if 0:
	fig = plt.figure()	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,4])
	
	ax1=plt.subplot(gs[0])
	ax0=plt.subplot(gs[1])

	for i in range(0,nProfiles):
		profNow = profileDataSwathOutliers['profiles'][i]
		ax0.scatter(np.array(termDist)/1e3,profNow['speed'],s=15,alpha=0.5,lw=0,c='darkorange')

	for i in range(0,nProfiles):
		profNow = profileDataSwathOutliers['profiles'][i]
		ax0.scatter(np.array(termDist)/1e3,profNow['speedMeanClean'],s=15,alpha=0.5,lw=0,c='darkred')

	for i in range(0,nProfiles):
		profNow = profileDataSwathOutliers['profiles'][i]
		ax0.scatter(np.array(termDist)/1e3,profNow['speedTemporalFilter'],s=15,alpha=0.5,lw=0,c='dodgerblue')

	ax0.plot(np.array(termDist)/1e3,meanRaw,lw=3,c='orange',alpha=0.8)
	ax0.plot(np.array(termDist)/1e3,meanSwath,lw=3,c='darkred',alpha=0.8)
	ax0.plot(np.array(termDist)/1e3,meanTemp,lw=3,c='mediumblue',alpha=0.8)
	ax0.fill_between(np.array(termDist)/1e3,meanRaw-stdRaw,meanRaw+stdRaw,lw=0,facecolor='orange',alpha=0.3)
	ax0.fill_between(np.array(termDist)/1e3,meanSwath-stdSwath,meanSwath+stdSwath,lw=0,facecolor='darkred',alpha=0.3)
	ax0.fill_between(np.array(termDist)/1e3,meanTemp-stdTemp,meanTemp+stdTemp,lw=0,facecolor='dodgerblue',alpha=0.5)


	ax1.plot(np.array(termDist)/1e3,np.sum(np.isfinite(uMatxRaw),axis=0),lw=3,c='orange')
	ax1.plot(np.array(termDist)/1e3,np.sum(np.isfinite(uMatxSwath),axis=0),lw=3,c='darkred')
	ax1.plot(np.array(termDist)/1e3,np.sum(np.isfinite(uMatxSwathOut),axis=0),lw=3,c='mediumblue')

	ax0.scatter((-10,-10),(-10,-10),alpha=0.5,lw=0,c='orange',label='Corr')
	ax0.scatter((-10,-10),(-10,-10),alpha=0.5,lw=0,c='darkred',label='Corr+Swath')
	ax0.scatter((-10,-10),(-10,-10),alpha=0.5,lw=0,c='dodgerblue',label='Corr+Swath+Temporal')

	ax0.set_ylim((-.1,3))
	ax0.set_xlim((0,38))
	#ax0.legend(loc=1,frameon=False)
	#plt.xlabel('Distance from terminus [km]',fontsize=18)
	ax0.set_ylabel('Velocity [m d$^{-1}$]',fontsize=18)
	ax0.invert_xaxis()
	ax0.set_xlabel('Distance from terminus [km]',fontsize=18)
	ax1.set_ylim((18,33))
	ax1.set_xlim((0,38))
	ax1.set_yticks(np.arange(20,35,4))
	ax1.set_xticks([])
	ax1.set_ylabel('# Obs. [-]',fontsize=18)
	#ax1.legend(loc=1,frameon=False)
	#plt.xlabel('Distance from terminus [km]',fontsize=18)
	#ax0.set_ylabel('Velocity [m d$^{-1}$]',fontsize=18)
	ax1.invert_xaxis()

	plt.savefig('filteringAndDataQuality.pdf')
	plt.show()
	plt.close()


	fig = plt.figure()	
	gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
	
	ax0=plt.subplot(gs[0])
	ax1=plt.subplot(gs[1])

	im0 = ax0.imshow(uMatxRaw[1:,:],aspect='auto',vmin=0,vmax=3,interpolation='none')
	#fig.colorbar(im0, orientation='vertical')

	ax1.fill_between(np.array(termDist)/1e3,0,stdRaw,lw=0.5,color='k',facecolor='khaki',alpha=0.9)
	ax1.plot(np.array(termDist)/1e3,meanRaw,c='orange',lw=3)
	ax1.invert_xaxis()

	
	ax0.set_xticks([])
	ax0.set_yticks([])
	ax0.set_ylabel('Time',fontsize=18)
	ax1.set_ylabel('$\mu$, $\sigma$ [m d$^{-1}$]',fontsize=18)
	ax1.set_xlabel('Distance from terminus [km]',fontsize=18)
	plt.savefig('test.pdf',dpi=500)
	plt.show()
	plt.close()

# swath filtering
profNow = profileDataSwathOutliers['profiles'][i]
gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	
ax0=plt.subplot(gs[0])
ax1=plt.subplot(gs[1])

ax0.scatter(np.array(termDist)/1e3,profNow['speedCount'],s=15,facecolor='cornflowerblue')
ax0b = ax0.twinx()
ax0b.plot(np.array(termDist)/1e3,profNow['speedStd'],lw=2,c='orange')
#ax0b.plot(np.array(termDist)/1e3,np.ones((len(termDist),1)),lw=1,ls='--',c='orange')
ax0.plot(np.array(termDist)/1e3,10*np.ones((len(termDist),1)),lw=1,ls='--',c='k')


ax1.plot(np.array(termDist)/1e3,profNow['speedMean'],ls='--',lw=0,c='darkred',marker='o',label='Raw')
ax1.plot(np.array(termDist)/1e3,profNow['speedMeanClean'],ls='-',lw=2,c='mediumblue',marker='o',label='Swath filter')

ax1.legend(loc=2,frameon=False)
ax1.set_xlim((0,38))
ax0.set_xlim((0,38))
ax0.set_ylim((-15,30))
ax0.set_xticklabels([])
ax1.set_ylabel('Velocity [m d$^{-1}$]',fontsize=18,color='k')
ax1.set_xlabel('Distance from terminus [km]',fontsize=18,color='k')
ax0.set_ylabel('Observation count [-]',fontsize=18,color='cornflowerblue')
ax0.tick_params(axis='y', colors='cornflowerblue')
ax0.spines['left'].set_color('cornflowerblue')
ax0b.set_ylabel('Standard deviation [m d$^{-1}$]',fontsize=18,color='orange')
ax0b.tick_params(axis='y', colors='orange')
ax0b.spines['right'].set_color('orange')


ax1.invert_xaxis()
ax0.invert_xaxis()
#ax0b.invert_xaxis()
ax0.set_yticks(np.arange(0,30,5))
plt.savefig('test.pdf')
plt.show()
plt.close()


