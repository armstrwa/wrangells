'''
Testing script to filter velocities based of downhill-direction.
Requires digitized centerline profiles.

Input = x-direction velocity, y-direction velocity, centerline profiles
Output = "cleaned x, y, and total velocities

05 Sep 2016
William Armstrong

'''

from osgeo import ogr, gdal
import numpy as np
import sys, os
import matplotlib.pyplot as plt

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
def thinVectors(x,y,vx,vy,vv,plotEveryX,plotEveryY=plotEveryX):
	xThin = x[::plotEveryX]
	yThin = y[::plotEveryY]
	YT, XT = np.meshgrid(xThin,yThin)
	vxT = vx[::plotEveryY,::plotEveryX]
	vyT = vy[::plotEveryY,::plotEveryX]	
	vvT = vv[::plotEveryY,::plotEveryX]	
	
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
shapefileIn='/Users/wiar9509/Documents/generalScripts/swath/shapefileInputs/newWStETransects_19nov15_fromQGIS_EGPSG102006.shp'

# Open shapefile
shp = ogr.Open(shapefileIn)
lyr = shp.GetLayer(0)
nFeat = lyr.GetFeatureCount()

# Iterate over centerline profiles
for i in range(0,nFeat):
	feat = lyr.GetFeature(i)
	transName = feat.GetField(1)
	print i, transName



