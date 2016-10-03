## SCRIPT TO ITERATE FILTERING PROCEDURE
#  Written 09 Nov 15 by William Armstrong
#  Meat of filtering procedure written by Karen Alley.
#  Mostly just modified to source WHA file directory


##### USER INPUTS #######
# Filter values
corrMin = 0.15 # Set a minimum correlation peak height
delcorrMin = 0.05 # Set a minimum difference between strongest and next strongest correlation peak height
doySpanMin = 0 # minimum days elapsed between images (set to 0 to not consider)
doySpanMax = 96 # Maximum days elapsed between images
noDataValue = -1.0 # No data value

# Filter toggles
# Turn to 1 to filter based on correlation height, difference in correlation height between 1st and 2nd strongest peaks
corrFlag = 1
delcorrFlag = 1
doySpanFlag = 1 # to only process if days elapsed between


# Where to read vv viles from
vvIndir = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/'
#accessoryIndirPrefix = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/' # root directory for correlations
 
# Where to write filtered files to 
outdir = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/'

outTag = 'filt_corr'+str(corrMin)+'_delcorr'+str(delcorrMin)+'.tif' # This will be the end of the filename of each output 

###### PROCESSING ######

# LOADING MODULES

import numpy as np
import gdal
import os
import matplotlib.pyplot as plt
from matplotlib import colors
import glob

def writeGeotiff(geoImageFilepath,writeArray,outFn):
	noDataValue = -1

	# Opening files			
	geoImage = gdal.Open(geoImageFilepath) # just used for georeferencing information
	#imageArray = image.ReadAsArray().astype(np.float32)

	# Generating output
	[cols,rows] = writeArray.shape
	proj        = geoImage.GetProjection()
	trans       = geoImage.GetGeoTransform()

	# Create the file, using the information from the original file
	outDriver = gdal.GetDriverByName("GTiff")
	outData   = outDriver.Create(str(outFn), rows, cols, 1, gdal.GDT_Float32)
	# Write the array to the file
	outData.GetRasterBand(1).WriteArray(writeArray)
	# Set a no data value if required
	outData.GetRasterBand(1).SetNoDataValue(noDataValue)
	# Georeference the image
	outData.SetGeoTransform(trans)
	# Write projection information
	outData.SetProjection(proj)
	# clear memory
	outData.FlushCache() 

def writeGeotiffRasterio(geoimageFn,arrayVals,outFn):
	'''
	Write a geotiff using rasterio. GDAL version wasn't working on lab cpu..
	Actually may have just been because the folder I was trying to write to didn't exist..
	But this much cleaner than old version using GDAL bindings
	'''
	with rasterio.open(geoimageFn) as src:
		prof = src.profile

	prof.update(nodata=-1,)

	with rasterio.open(outFn, 'w', **prof) as dst:
		dst.write(arrayVals, 1)


# FILTERING FUNCTION
def velocityFilter(corrMin,delcorrMin,noDataValue,vvPath,outdir):

	vvFn = vvPath.split('/')[-1]
	
	print "Filtering file: " + vvFn
	print "Filter values: corrMin = " + str(corrMin) + "; delcorrMin = " + str(delcorrMin) 
	
	fn1 = vvFn.split('_')[0] # first image filename for later manipulation
	fn2 = vvFn.split('_')[1] # first image filename
	yearDoy = vvFn.split('_')[2] # correlation midpoint year-doy
	doySpan = vvFn.split('_')[3] # days elapsed between images
	wrsPath = fn1[3:6] # path number from filename
	wrsRow = fn1[6:9] # row number from filename
	sourceChip = vvFn.split('_')[4] # size of half source chip
	targetChip = vvFn.split('_')[5] # size of half source chip
	increment = vvFn.split('_')[6] # increment step size
	gfiltSigma = 3 # HARDCODED b/c this not in the filename but does exist in the directory structure
	
	# Constructing correlation filepath
	imageImagePath = fn1+'_'+fn2+'_'+str(yearDoy)+'_'+str(doySpan)+'_hp' 
	correlationSuffix = '_' + str(gfiltSigma) + '_' + str(sourceChip) + '_' + str(targetChip) + '_' + str(increment)

	# Get associated corr and delcorr files
	corrFileFullpath = glob.glob(vvIndir + fn1 + "_" + fn2 + '*corr.tif')[0]
	delcorrFileFullpath = glob.glob(vvIndir + fn1 + "_" + fn2 + '*delcorr.tif')[0]
	vxFileFullpath = glob.glob(vvIndir + fn1 + "_" + fn2 + '*vx.tif')[0]
	vyFileFullpath = glob.glob(vvIndir + fn1 + "_" + fn2 + '*vy.tif')[0]
			
	# Opening files			
	gdCorr = gdal.Open(corrFileFullpath)
	gdDelcorr = gdal.Open(delcorrFileFullpath)
	gdVx = gdal.Open(vxFileFullpath)
	gdVy = gdal.Open(vyFileFullpath)
	gdVv = gdal.Open(vvPath)
	
	# Reading files
	arrCorr = gdCorr.ReadAsArray().astype(np.float32)
	arrDelcorr = gdDelcorr.ReadAsArray().astype(np.float32)
	arrVx = gdVx.ReadAsArray().astype(np.float32)
	arrVy = gdVy.ReadAsArray().astype(np.float32)
	arrVv = gdVv.ReadAsArray().astype(np.float32)
	arrVvUnfiltered = gdVv.ReadAsArray().astype(np.float32)    
	arrVx[np.where(arrVx == 0)] = noDataValue  # Sets no data values from image to script-specified no data value
	arrVy[np.where(arrVy == 0)] = noDataValue
	arrVv[np.where(arrVv == 0)] = noDataValue
	
	# Filtering on correlation peak height
	if corrFlag == 1:
		arrVx[np.where(arrCorr<corrMin)] = noDataValue
		arrVy[np.where(arrCorr<corrMin)] = noDataValue
		arrVv[np.where(arrCorr<corrMin)] = noDataValue

	# Filtering on difference in correlation peak height
	if delcorrFlag == 1:
		arrVx[np.where(arrDelcorr<delcorrMin)] = noDataValue
		arrVy[np.where(arrDelcorr<delcorrMin)] = noDataValue
		arrVv[np.where(arrDelcorr<delcorrMin)] = noDataValue

	outVxFn = os.path.join(outdir,vxFileFullpath.split('/')[-1].split('hp_')[0] + 'vx_' + outTag)
	outVyFn = os.path.join(outdir,vyFileFullpath.split('/')[-1].split('hp_')[0] + 'vy_' + outTag)
	outVvFn = os.path.join(outdir,vvFn.split('hp_filt_3.0_vv.tif')[0] + 'vv_' + outTag)
	
	writeGeotiffRasterio(corrFileFullpath,arrVx,outVxFn)
	writeGeotiffRasterio(corrFileFullpath,arrVy,outVyFn)
	writeGeotiffRasterio(corrFileFullpath,arrVv,outVvFn)
			

# IDENTIFYING FILES FOR ITERATION
allFiles = os.listdir(vvIndir)

corrFiles = []
delcorrFiles = []
vxFiles = []
vyFiles = []
vvFiles = []
				
# ITERATING
for i in range(len(allFiles)):
#for i in range(0,10):
	# Identify vv files
	if allFiles[i].endswith('vv.tif'):
		vvFn = allFiles[i]
		vvPath = os.path.join(vvIndir,vvFn)
		
		# Run the filter
		velocityFilter(corrMin,delcorrMin,noDataValue,vvPath,outdir)
	