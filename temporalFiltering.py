'''
Developing temporal filtering for correlation quality
22 november 2016
wha
'''


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

jsonFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/rootCL_evenlySpace_100m_profiles_sampled_2015-12-04.json'

data = wgt.readJsonFile(jsonFn)


# Make velocity matrix with len(time) rows and len(space) cols
uMatx = np.ones((1,len(data['sample_pts_frontdist'])))*np.nan
numProfiles = len(data['profiles'])
midDoy = []
for i in range(0,numProfiles):
	profNow = data['profiles'][i]
	uNow = np.array(profNow['speed'])
	dDays = profNow['delts_days']
	midDoy.append(profNow['mid_dec_year'])
	uNow[uNow==-1] = np.nan
	uMatx = np.vstack((uMatx,uNow))
	

du_time = np.diff(uMatx,axis=0)
du_space = np.diff(uMatx,axis=1)

if 0:
	plt.subplot(3,1,1)
	plt.title('u')
	plt.imshow(uMatx,clim=(0,2))
	#plt.yticks(midDoy)
	plt.colorbar()
	plt.subplot(3,1,2)
	plt.title('du_time')
	plt.imshow(du_time,clim=(-.5,.5))
	plt.colorbar()
	plt.subplot(3,1,3)
	plt.title('du_space')
	plt.imshow(du_space,clim=(-.25,.25))
	plt.colorbar()
	plt.savefig('test.png')
	plt.show()
	plt.close()

# iterative removal of outliers by mean/stdev
def iterateMeanRemove(uMatx):
	mean1 = np.nanmean(uMatx,axis=0)
	std1 = np.nanstd(uMatx,axis=0)
	uMatx2 = np.copy(uMatx)
	tol = 2*std1
	for rowNum in range(0,np.size(uMatx,0)):
		ind1 = np.where(np.logical_or(uMatx[rowNum,:]<=mean1-tol,uMatx[rowNum,:]>=mean1+tol))
		uMatx2[rowNum,ind1] = np.nan

	return uMatx2

mean1 = np.nanmean(uMatx,axis=0)
std1 = np.nanstd(uMatx,axis=0)
uMatx2 = iterateMeanRemove(uMatx)
mean2 = np.nanmean(uMatx2,axis=0)
std2 = np.nanstd(uMatx2,axis=0)
uMatx3 = iterateMeanRemove(uMatx2)
mean3 = np.nanmean(uMatx3,axis=0)
std3 = np.nanstd(uMatx3,axis=0)


plt.plot(mean1,lw=2,c='r',label='$\mu$ raw')
plt.plot(mean2,lw=2,c='g',label='$\mu$ 1 iter')
plt.plot(mean3,lw=2,c='b',label='$\mu$ 2 iter')

plt.plot(std1,lw=2,c='r',ls='--',label='$\sigma$ raw')
plt.plot(std2,lw=2,c='g',ls='--',label='$\sigma$ 1 iter')
plt.plot(std3,lw=2,c='b',ls='--',label='$\sigma$ 2 iter')

plt.title('temporal filtering: iterative mean and stdev')
plt.legend(frameon=False)
plt.savefig('temporalFiltering_iterativeMeanAndStdev.pdf')
plt.show()
plt.close()

plt.subplot(3,1,1)
plt.imshow(uMatx,clim=(0,2))
plt.title('Raw data')
plt.subplot(3,1,2)
plt.imshow(uMatx2,clim=(0,2))
plt.title('One iteration, remove 2$\sigma$')
plt.subplot(3,1,3)
plt.imshow(uMatx3,clim=(0,2))
plt.title('Two iterations, remove 2$\sigma$')
plt.savefig('temporalFiltering_iterativeMeanAndStdev_heatmap.png')
plt.show()
plt.close()
