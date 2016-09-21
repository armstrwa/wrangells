'''
Script to filter velocity based on downhill direction for glacier

01 Mar 2016
William Armstrong

'''

from osgeo import ogr, gdal
import numpy as np
import sys, os
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt

def gaussian_blur(in_array, size):
# Function from http://gis.stackexchange.com/questions/9431/what-raster-smoothing-generalization-tools-are-available
    # expand in_array to fit edge of kernel
    padded_array = np.pad(in_array, size, 'symmetric')
    # build kernel
    x, y = np.mgrid[-size:size + 1, -size:size + 1]
    g = np.exp(-(x**2 / float(size) + y**2 / float(size)))
    #g = (g / g.sum()).astype(in_array.dtype)
    g = (g / g.sum()).astype('float32') # modified b/c floor divide
    # do the Gaussian blur
    return fftconvolve(padded_array, g, mode='valid')

# Function to write a geotiff
def writeGeotiff(geoImageFilepath,writeArray,outFn):
	noDataValue = -999

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

# Function to thin vector field
# Inputs = x coord vector, y coord vector, vx matx, vy matx, vv matx, plot frequency in x, plot frequency in y
def thinVectors(x,y,vx,vy,vv,plotEvery):
	xThin = x[::plotEvery]
	yThin = y[::plotEvery]
	YT, XT = np.meshgrid(xThin,yThin)
	vxT = vx[::plotEvery,::plotEvery]
	vyT = vy[::plotEvery,::plotEvery]	
	vvT = vv[::plotEvery,::plotEvery]	
	
	return XT, YT, vxT, vyT, vvT

# Creates transformation parameters for future work
# Outputs xyTfm vector that is [minX,maxY,pxSizeX,pxSizeY]
def createTransformation(rasterIn):
	rast = gdal.Open(rasterIn) # open raster
	gt = rast.GetGeoTransform() # get image to ground geotransformation params
	numPxX = rast.RasterXSize # number of pixels in x direction
	numPxY = rast.RasterYSize # number of pixels in y direction
	pxSizeX = gt[1] # pixel size in x
	pxSizeY = gt[5] # pixel size in y

	# upper left coordinates
	minX = gt[0] # minimum x value
	maxY = gt[3] # maximum y value

	maxX = minX + numPxX*pxSizeX # maximum x value
	minY = maxY + numPxY*pxSizeY # minimum y value

	xyTfm = np.array([minX,maxY,pxSizeX,pxSizeY]).astype('float')

	return xyTfm
	
# Function to convert pixel location to ground location
def imageToXy(i,j):
	
	outX = []
	outY = []
	
	it = np.nditer([i,j])
	
	for ii,jj in it:
		outX.append( xyTfm[0] + ( (ii+0.5) * xyTfm[2] ) )
		outY.append( xyTfm[1] + ( (jj+0.5) * xyTfm[3] ) )
	
	outArray = np.array((outX,outY))
	outArrayTranspose = outArray.transpose()
	return outArrayTranspose

# Function to convert ground coordinates to image pixel location
def xyToImageIj(x,y):
	
	outI = []
	outJ = []
	
	it = np.nditer([x,y])
	
	for xx,yy in it:
		outI.append( (xx - xyTfm[0])/xyTfm[2] - 0.5 )
		outJ.append( (yy - xyTfm[1])/xyTfm[3] - 0.5 )
	
	outArray = np.array((outI,outJ))
	outArrayTranspose = outArray.transpose()
	return outArrayTranspose

# Function to get coordinates from shapefile from study glaciers shapefile
# Outputs array of x, y coordinates of shapefile vertices
def getShapefileCoordinatesFromMultiline(shapefileIn,transName=None):
	driver = ogr.GetDriverByName("ESRI Shapefile")
	ds = driver.Open(shapefileIn,0)
	lyr = ds.GetLayer(0)
	
	numFeats = lyr.GetFeatureCount()
	
	for i in range(0,numFeats):
		featNow = lyr.GetFeature(i)
		transNow = featNow.GetField(1)
		
		if transName is None or transName == transNow:
			
			geom = featNow.geometry()
			numPoints = geom.GetPointCount()
			
			x = []
			y = []
			
			for j in range(0,numPoints):
				x.append( geom.GetX(j) )
				y.append( geom.GetY(j) )
		
	coordArray = np.array((x,y))
	coordArrayTranspose = coordArray.transpose()
		
	return coordArrayTranspose
	
# Function to get coordinates from shapefile from evenly space shapefile
# Outputs array of x, y coordinates of shapefile vertices
def getShapefileCoordinatesFromMultipoint(shapefileIn):
	driver = ogr.GetDriverByName("ESRI Shapefile")
	ds = driver.Open(shapefileIn,0)
	lyr = ds.GetLayer(0)
	
	numPoints = lyr.GetFeatureCount()
			
	x = []
	y = []
	
	for i in range(0,numPoints):
		pointNow = lyr.GetFeature(i)
		geom = pointNow.geometry()
		x.append( geom.GetX() )
		y.append( geom.GetY() )
		
	coordArray = np.array((x,y))
	coordArrayTranspose = coordArray.transpose()
		
	return coordArrayTranspose

# Function to sample raster at specified coordinates
# Returns values of raster at x,y
def sampleRasterAtXY(rasterIn,x,y):
	
	z = []
	
	imageIj = xyToImageIj(x,y)
	
	rast = gdal.Open(rasterIn)
	band = rast.GetRasterBand(1) # assumes single band raster
	bandArr = band.ReadAsArray() # convert band into array
	
	samplePtNum = imageIj.shape[0] # number of sample points
	
	for i in range(0,samplePtNum):
		z.append( bandArr[ np.round(imageIj[i,1]), np.round(imageIj[i,0] ) ] )
	
	return z
	
# Calculate the distance along a line, given vertex coordinates
def distanceAlongLine(x,y):
	
	dist = [0] # starts at 0 distance
	
	numPoints = len(x)
	
	for i in range(1,numPoints):
		oldEasting = x[i-1]
		oldNorthing = y[i-1]
		nowEasting = x[i]
		nowNorthing = y[i]
		dE = nowEasting - oldEasting
		dN = nowNorthing - oldNorthing
		distChange = np.sqrt(dE**2 + dN**2)
		dist.append( distChange + dist[i-1] )
		
	return dist


# shapefile containing glacier profiles
#shapefileIn='/Users/wiar9509/Documents/generalScripts/swath/shapefileInputs/newWStETransects_19nov15_fromQGIS_EGPSG102006.shp'
#demIn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/dtm/wrangellStElias_gdem_EPSG102006.tif'

path = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/'
vxFile = path + 'LC80630172014194LGN0_LC80630172014290LGN0_2014-242_96_24_40_10_vx_filt_corr0.15_delcorr0.05_EPSG102006.tif'
vyFile = path + 'LC80630172014194LGN0_LC80630172014290LGN0_2014-242_96_24_40_10_vy_filt_corr0.15_delcorr0.05_EPSG102006.tif'

# Flow direction
vxImage = gdal.Open(vxFile)
vyImage = gdal.Open(vyFile)
vx = vxImage.ReadAsArray()
vy = vyImage.ReadAsArray()
vv = np.sqrt(vx**2 + vy**2)
theta = np.arctan2(vx,vy)

dx = 150 # pixel size in x dirn
dy = 150 # pixel size in y dirn
plotEvery=10 # plot every 10th velocity vector
xObs = vxImage.RasterXSize
yObs = vxImage.RasterYSize
#x = np.round(np.linspace(0,xObs*dx,np.ceil(float(xObs)/plotEvery)))
x = np.round(np.linspace(0,xObs*dx,xObs))
y = np.round(np.linspace(0,yObs*dy,yObs))
Y, X = np.meshgrid(x,y)

XThin, YThin, vxThin, vyThin, vvThin = thinVectors(x,y,vx,vy,vv,5)

uMax = 1
vxThinClean = np.copy(vxThin)
vyThinClean = np.copy(vyThin)
vvThinClean = np.copy(vvThin)
vxThinClean[np.where(vvThinClean>uMax)] = np.nan
vyThinClean[np.where(vvThinClean>uMax)] = np.nan
vvThinClean[np.where(vvThinClean>uMax)] = np.nan

scaleFactor = 10
#lws = np.linspace(0,2,YThin.size)
plt.imshow(theta,extent=[x.min(), x.max(), y.min(), y.max()])
plt.colorbar()
#Q = plt.quiver(XThin,YThin,vxThinClean,vyThinClean,vvThinClean,width=1e2,units='xy',angles='xy',scale=2.5e-4)
Q = plt.quiver(XThin,YThin,vxThinClean,vyThinClean,color='white',width=1e2,units='xy',angles='xy',scale=2.5e-4)
# qk = plt.quiverkey(Q, 0.5, 0.03, 1, r'$1 \frac{m}{s}$',
#                    fontproperties={'weight': 'bold'})
#plt.colorbar()
plt.show()
plt.close()


# DEM things
#demIn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/dtm/stelias_gdem.tif'
demIn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/dtm/tanaDemMerge.tif'
dem = gdal.Open(demIn) # Open DEM
demArr = dem.ReadAsArray() # convert tif to array
demSm = gaussian_blur(demArr,20)

meterInDeg = 1/111412.24 # length of a meter in degrees at 60 deg N
spatialRes = 500

# Resample DEM
os.system('gdalwarp -tr ' + str(meterInDeg*spatialRes) + ' ' + str(meterInDeg*spatialRes) + ' ' + demIn + ' ' +demIn[:-4] + '_resample' + str(spatialRes) + 'm.tif')

if 1:
	fig = plt.figure(0)
	plt.subplot(1,2,1)
	plt.imshow(demArr)
	plt.clim(0,2000)
	plt.colorbar()
	
	plt.subplot(1,2,2)
	plt.imshow(demSm)
	plt.clim(0,2000)
	plt.colorbar()	
	plt.show()
	plt.close()

outPrefix = demIn[:-4]
writeGeotiff(demIn,demSm,outPrefix+'_gblur20.tif')

# Processing
aspectOutFn = demIn[0:-4] + '_aspect_rs' + str(spatialRes) + '.tif'
slopeOutFn = demIn[0:-4] + '_slope_rs'  + str(spatialRes) + '.tif'
#hillshadeOutFn = demIn[0:-4] + '_hsd_rs750.tif'
demFn = outPrefix + '_resample'  + str(spatialRes) + 'm.tif'
	
# Create slope and aspect maps
if not os.path.isfile(aspectOutFn): # Don't run if file already exists
	print 'Running gdaldem aspect'
	#os.system('gdaldem aspect ' + demIn + ' ' + aspectOutFn)
	os.system('gdaldem aspect ' + demFn + ' ' + aspectOutFn)
if not os.path.isfile(slopeOutFn): # Don't run if file already exists	
	print 'Running gdaldem slope'
	#os.system('gdaldem slope -s 111120 ' + demIn + ' ' + slopeOutFn) # -s is scale factor if horizontal and vertical units differ
	os.system('gdaldem slope -s 111120 ' + demFn + ' ' + slopeOutFn) # -s is scale factor if horizontal and vertical units differ
if not os.path.isfile(hillshadeOutFn): # Don't run if file already exists	
	os.system('gdaldem hillshade -z 0.00001792 ' + demFn + ' ' + hillshadeOutFn) # output hillshade
	
# Read in aspect map
asp = gdal.Open(aspectOutFn) # open file
aspArr = asp.ReadAsArray() # turn to array
northness = np.sin(np.deg2rad(aspArr)) # calculate northness (sine of aspect)
eastness = np.cos(np.deg2rad(aspArr)) # calculate eastness (cosine of aspect)

	