## SCRIPT TO CALCULATE NORMALIZED DIFFERENCE SNOW INDEX
#  29 Nov 2015
#  William Armstrong

### PROCESSING ###

from osgeo import gdal, ogr
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from matplotlib import colors

def calculateNDVI(imageDirectory,saveImageToggle=0,overwrite=0):
	'''
	Calculate normalized difference snow index for landsat multiband data
	Input: imageDirectory = directory containing landsat bands for one day at one path/row;
	Input: saveImageToggle = 0 to not output a geotiff
	Input: overwrite = 0 to not overwrite an existing ndsi image;
	Output: ndsi = np.array of snow index values
	'''
	# Test for existence
	outFn = imageDirectory + imageDirectory.split('/')[-2] + '_ndsi.TIF'
	
	if os.path.isfile(outFn) and overwrite == 0 and saveImageToggle == 1:
		print "File already exists: " + outFn	
	else:
		# Locate green and swir files in directory
		greenFn = glob.glob(imageDirectory+'*_B3.TIF')[0]
		swirFn = glob.glob(imageDirectory+'*_B6.TIF')[0]

		# Also more bands needed for true color
		redFn = glob.glob(imageDirectory+'*_B4.TIF')[0]
		blueFn = glob.glob(imageDirectory+'*_B2.TIF')[0]

		# Let's just to panchromatic for simplicity
		#panFn = glob.glob(imageDirectory+'*_B8.TIF')[0]

		# Open files
		greenRast = gdal.Open(greenFn)
		swirRast = gdal.Open(swirFn)
		#panRast = gdal.Open(panFn)

		green = greenRast.ReadAsArray().astype(np.float)
		swir = swirRast.ReadAsArray().astype(np.float)
		#pan = panRast.ReadAsArray().astype(np.float)

		# Calculate NDSI
		ndsi = (green-swir)/(green+swir)

		# NDSI histogram
		#hist,edges = np.histogram(ndsi,range=(-1,1),bins=100)
		#plt.plot(edges[1:],hist)
		#plt.show()

		# Delineate snow
		snowInd = np.where(ndsi>=snowThreshold)
		snow = np.copy(ndsi)
		snow[:,:] = 0
		snow[snowInd] = 1

		# Delineate ice
		iceInd = np.logical_and(ndsi>=iceThreshold,ndsi<snowThreshold)
		ice = np.copy(ndsi)
		ice[:,:] = 0
		ice[iceInd] = 0.5

		# mask
		#G = np.copy(ndsi)
		G = np.zeros((ndsi.shape[0],ndsi.shape[1],3))
		G[:,:] = [0,0,0]
		G[ice>0] = [1,0,0]
		G[snow>0] = [0,0,1]

		### SAVING AND PLOTTING ###

		if plotImageToggle == 1:
			# Plot
			plt.subplot(121)
			plt.imshow(ndsi,vmin=0,vmax=1)
			plt.title('NDSI')
			#plt.colorbar()

			plt.subplot(122)
			plt.imshow(G,interpolation='nearest')
			plt.title('Snow = blue, ice = red')
			plt.show()
			plt.close()
	
		if saveImageToggle == 1:
					[cols,rows] = ndsi.shape
					trans       = greenRast.GetGeoTransform()
					proj        = greenRast.GetProjection()
			
					outFn = imageDirectory + imageDirectory.split('/')[-2] + '_ndsi.TIF'
					print "Writing out file: ", outFn
					
					# Create the file, using the information from the original file
					outdriver = gdal.GetDriverByName("GTiff")
					outData   = outdriver.Create(str(outFn), rows, cols, 1, gdal.GDT_Float32, options = [ 'COMPRESS=LZW' ])


					# Write the array to the file, which is the original array in this example
					outData.GetRasterBand(1).WriteArray(ndsi)

					# Set a no data value if required
					outData.GetRasterBand(1).SetNoDataValue(-1)

					# Georeference the image
					outData.SetGeoTransform(trans)

					# Write projection information
					outData.SetProjection(proj)
			
					outData.FlushCache() 
			
		return ndsi 



## USER INPUTS ###
#imageDirectory = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/imagery/landsat8/imagery/64_17/LC80640172013166LGN00/'
#imageDirectory = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/imagery/landsat8/imagery/64_17/LC80640172014121LGN00/'
#imageDirectory = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/imagery/landsat8/imagery/64_17/LC80640172014233LGN00/'
#imageDirectory = '/Volumes/oldEthan/l8_imagery/imagery/63_17/LC80630172013239LGN00/'
snowThreshold = 0.604
iceThreshold = 0.25
saveImageToggle = 1
plotImageToggle = 0


folderPath = '/Users/anderson/Desktop/ARMSTRONG/wrangells/imagery/L8/62_17/'

folderList = os.listdir(folderPath)

for folder in folderList:
	if os.path.isdir(folderPath + folder):
		ndsi = calculateNDVI(folderPath + folder + '/')



