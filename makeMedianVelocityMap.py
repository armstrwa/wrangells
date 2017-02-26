import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal

rasterPath = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/select/'

pathRow = '064017'

fileList = glob.glob(rasterPath + '*' + pathRow + '*vv' + '*.tif')

# create raster stack
for i in range(0,len(fileList)):

	fileNow = fileList[i]
	
	rast = gdal.Open(fileNow)
	arr = rast.ReadAsArray()
	
	# initialize
	if i == 0:
		rasterStack = arr
	else:
		rasterStack = np.dstack((rasterStack,arr))
	
medianPixelValues = np.median(rasterStack,axis=2)



#computeMedianVelocitiesPixelwise(rasterList)