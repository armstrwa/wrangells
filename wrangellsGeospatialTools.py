'''
New script for swath profiler that does not require GRASS
begun 06 Sep 2016
William Armstrong

Inputs: line geometry, raster
Outputs: raster swath profiled along line geometry

'''

## -------------- ##
## IMPORT MODULES ##
## -------------- ##

import os
# import sys
import numpy as np
# import scipy as sci 
# import csv
import glob
from functools import partial
import matplotlib.pyplot as plt
import json
from shapely import geometry
from shapely.ops import transform
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

# Read json file
# Input = full file path for json data
# Output = json dict
def readJsonFile(jsonFn):
	with open(jsonFn) as data_file:
		data = json.load(data_file)
	return data

# Write a dict to a json file
def writeJsonFile(dataToWrite,jsonFn):
	with open(jsonFn, 'w') as f: 
		f.write(json.dumps(dataToWrite))


# Get coordinates for a line given its centroid
# Stolen from krb_vectools
def makeLine(centerPt, slope, halfLength, slopeUnits="radians"):
    if slopeUnits=="degrees":
        # convert to radians
        slope=slope/360*2*np.pi
    
    dx=np.cos(slope)*halfLength
    dy=np.sin(slope)*halfLength    
    firstPt=[centerPt[0]-dx, centerPt[1]-dy]
    lastPt=[centerPt[0]+dx, centerPt[1]+dy]
    line=np.vstack((firstPt, centerPt, lastPt))
    return line


# Get coordinates of line geometry vertices
def getLineGeometryVerticesXY(lineGeom):
	nPoints = lineGeom.GetPointCount()
	clX = []
	clY = []
	
	for j in range(0,nPoints):
		clX.append(lineGeom.GetX(j))
		clY.append(lineGeom.GetY(j))	
		
	return clX, clY

# Get coordinates of multipoint geometry 
def getMultipointGeometryVerticesXY(multipointGeom):
	nGeom = multipointGeom.GetGeometryCount()
	
	clX = []
	clY = []
	
	for i in range(0,nGeom):
		pointGeom = multipointGeom.GetGeometryRef(i)
		clX.append(pointGeom.GetX(0))
		clY.append(pointGeom.GetY(0))	
		
	return clX, clY


# Fn to interpolate evenly spaced vertices on a polyline
# As of 13 Jul 2016, these not evenly spaced, but I'm not super concerned about that right now.
# As of 06 Sep 2016, xInt and yInt aren't the same length? yInt one longer.. Appears to be dependent on dDist. Works correctly for some..
# 02 oct 2016, this causing problems... wha
def interpolateVerticesOfPolyline(lineGeom,dDist):
	# initialize/clear if exist
	clX = []
	clY = []
	
	clX,clY = getLineGeometryVerticesXY(lineGeom) # line vertices in utm
	dist = distanceAlongLine(clX,clY) # distance along centerline
	nP = len(clX) # number of vertices
	
	for i in range(1,nP): # iterate over vertices
#		print i
			
		nPointSeg = np.round(np.array(dist[i]-dist[i-1])/dDist) # number of points that will be on line segment
		dE = clX[i] - clX[i-1] # change in easting across line seg
		dN = clY[i] - clY[i-1] # change in northing across line seg
		dx = np.array(dE)/nPointSeg # change in easting between interpolated points
		dy = np.array(dN)/nPointSeg # change in easting between interpolated points
		xSeg = clX[i-1] + np.arange(0,dE,dx) # x-coordinates for interpolated points on line
		ySeg = clY[i-1] + np.arange(0,dN,dy) # y-coordinates for interpolated points on line
		
		# This very hackish way to delete with np.arange not behaving like expected. Sometimes len(xSeg) != len(ySeg) because of last point.
		# Seemed like this happened on the first point (i.e., i = 1). wha 26sep2016
		if len(ySeg) > len(xSeg):
			ySeg = np.delete(ySeg,-1)
		elif len(xSeg) > len(ySeg):
			xSeg = np.delete(xSeg,-1)
			
		#print i, len(xSeg), len(ySeg)
		
		if i == 1: # special treatment for first point
			xInt = xSeg
			yInt = ySeg
		else: # append new points to old
			xInt = np.append(xInt,xSeg)
			yInt = np.append(yInt,ySeg)
		
	return xInt, yInt


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
def imageToXy(i,j,xyTfm):
	
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
# Returns i,j pairs. I.e., outArrayTranspose[:,0] = i (col) coord,
# outArrayTransport[:,1] = j (row) coord
def xyToImageIj(x,y,xyTfm):
	
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


# Function to get spatial reference from a raster
def getSrsFromRaster(rasterFn):
	raster = gdal.Open(rasterFn)
	rastSrs = osr.SpatialReference()
	rastProj = raster.GetProjectionRef()
	rastSrs.ImportFromWkt(rastProj)
	return rastSrs

# Function to get spatial reference from a shapefile
def getSrsFromShapefile(shapefileIn):
	shp = ogr.Open(shapefileIn)
	lyr = shp.GetLayer(0)
	lyrSrs = lyr.GetSpatialRef()
	return lyrSrs

# Function to test for matching spatial reference
# Returns 1 if true (matching), 0 if false (not matching)
def testForMatchingSpatialReference(fn1,fn2):

	# Get spatial reference of first file
	if fn1[-4:] == '.shp':
		srs1 = getSrsFromShapefile(fn1) # spatial reference if vector
	elif fn1[-4:] == '.tif' or fn1[-4:] == 'tiff':
		srs1 = getSrsFromRaster(fn1) # spatial reference if raster
	
	# Get spatial reference of second file
	if fn2[-4:] == '.shp':
		srs2 = getSrsFromShapefile(fn2) # spatial reference if vector
	elif fn2[-4:] == '.tif' or fn2[-4:] == 'tiff':
		srs2 = getSrsFromRaster(fn2) # spatial reference if raster
	
	return srs1.IsSame(srs2) 

# Function to subset a raster by specified ground coordinates
def subsetRasterByXy(boundX,boundY,rasterFn):
	xyTfm = createTransformation(rasterFn) # transformation parameters for converting raster pixels to ground coordinates
	boundIj = xyToImageIj(boundX,boundY,xyTfm) # image coordinates of bounding box
	boundI = boundIj[:,0] # i (col) coords
	boundJ = boundIj[:,1] # j (row coords)
	dI = int(np.round(boundI[1] - boundI[0])) # change in i
	dJ = np.abs(int(np.round(boundJ[1] - boundJ[0]))) # change in j
	rastArr = raster.ReadAsArray(int(boundI[0]),int(boundJ[1]),dI,dJ) # subsetted raster

	return rastArr

# stolen from krb_vectools
def calculateSlope(pt1,pt2, units="radians"):
    
    ''' calculates slope between horizontal and line defined by p1=[x1,y1] and 
    p2=[x2,y2] at the point defined by [x1, y1] using np.arctan2'''
    
    if (np.asarray(pt1).shape[0]==np.asarray(pt2).shape[0]) and (np.asarray(pt1).size%2==0 ) and(np.asarray(pt2).size%2==0 ):   
    
        if np.asarray(pt1).size>2: 
            # reqires unzipping...
            x1, y1=zip(*pt1)    
            x2, y2=zip(*pt2)
            x=[]
            y=[]
            # f-you tuples...
            for i in range(len(x1)):
                x.append(x2[i]-x1[i])
                y.append(y2[i]-y1[i])                       
        else:
            x=pt2[0]-pt1[0]
            y=pt2[1]-pt1[1]  
            
        slope=np.arctan2(y,x)    
        if units=="degrees":
            # convert to degrees
            slope=slope/(2*np.pi)*360        
        return slope


# Calculate slope of line segments
def calculateSlopes(xArr,yArr):
	slopes = []
	for i in range(1,len(xArr)):
		pt2 = (xArr[i],yArr[i])
		pt1 = (xArr[i-1],yArr[i-1])
		slopes.append(calculateSlope(pt1,pt2))

	return slopes
	
# Calculate perpendicular slope of line segments
def calculatePerpendicularSlopes(slopes):
	if type(slopes) == list: # turn to array if list
		slopes = np.array(slopes)
	thetaPerp = ((slopes+2*np.pi)-(np.pi/2))%(2*np.pi) # rotate slope by pi/2
	return thetaPerp

# Make vertices for one box for sampling
def makeOneBoxVertices(centersX,centersY,perpSlope,boxWidth):
	topCoords = makeLine((centersX[0],centersY[0]), perpSlope,boxWidth/2.)
	bottomCoords = makeLine((centersX[1],centersY[1]), perpSlope,boxWidth/2.)
	boxCoords = np.array((topCoords[0],topCoords[2],bottomCoords[2],bottomCoords[0]))

	return boxCoords

# Make vertices for all boxes along line	
def makeAllBoxVertices(xInt,yInt,boxWidth):
	verts=[] # intitalize
	slopes=calculateSlopes(xInt,yInt) # calculate "slope" (azimuth) of line segments
	perpSlopes=calculatePerpendicularSlopes(slopes) # calculate perpendicular "slope" (azimuth)
	
	for i in range(0,len(xInt)-1):
		centersX = np.array((xInt[i],xInt[i+1]))
		centersY = np.array((yInt[i],yInt[i+1]))		
		vertsNow = makeOneBoxVertices(centersX,centersY,perpSlopes[i],boxWidth) # get verticies of box
		verts.append(vertsNow) # append to list
		
	return verts

# Function to make a multipolygon vector file from corner vertices of boxes
# Input: List of list containing x,y coords for each box corner. boxVertices[0] = 4x2 array containing x,y coords of one box
# srs = spatial reference (osr.SpatialReference()), shpFnOut = filename for output shapefile
def makeMultipolygonFromBoxVertices(boxVertices,srs,shpFnOut):
	driver = ogr.GetDriverByName("ESRI Shapefile") # get driver
	if os.path.exists(shpFnOut):
		driver.DeleteDataSource(shpFnOut) # error if data source already exists
	newDataSource=driver.CreateDataSource(shpFnOut) # create datasaource
	newLayer=newDataSource.CreateLayer("samplingPolygons",srs,geom_type=ogr.wkbMultiPolygon) # create layer
	fieldDefn=ogr.FieldDefn('id',ogr.OFTInteger)  # create field for polygon id
	newLayer.CreateField(fieldDefn)	
	
	numBoxes = len(boxVertices) # number of boxes to make

	# Iterating to create polygons and set fields
	for i in range(numBoxes):

		# create the feature
		feature = ogr.Feature(newLayer.GetLayerDefn())

		# Set Geometry
		ring=ogr.Geometry(ogr.wkbLinearRing)
		ring.AddPoint(boxVertices[i][0,0],boxVertices[i][0,1]) # adding easting/northing for each vertex
		ring.AddPoint(boxVertices[i][1,0],boxVertices[i][1,1])
		ring.AddPoint(boxVertices[i][2,0],boxVertices[i][2,1])
		ring.AddPoint(boxVertices[i][3,0],boxVertices[i][3,1])
		ring.CloseRings()

		polygon = ogr.Geometry(ogr.wkbPolygon)
		polygon.AddGeometry(ring)
		feature.SetGeometryDirectly(polygon)

		# Setting fields
		feature.SetField('id',i) # id number
#		featDist=lineDist[i]
#		feature.SetField('dist',featDist) # distance along line

		# Create the feature in the layer (shapefile)
		newLayer.CreateFeature(feature)

		# Destroy the feature to free resources
		feature.Destroy()

	# Destroy the data source to free resources
	newDataSource.Destroy()
	
	print "Sucessfully wrote shapefile: " + str(shpFnOut)
	

def makePolylineShpFromCoords(x,y,epsgNum,shpFnOut):
	'''
	# Function uses shapely and fiona to make a polyline form coords
	# Inputs = x coordinate list, y coordinate list, epsgNum = epsg coordinate system code, schema = dict defining shapefile properties
	# Modified from Mike T answer on http://gis.stackexchange.com/questions/52705/how-to-write-shapely-geometries-to-shapefiles
	# As of 08 sep 2016 doesn't encode crs because missing fiona.crs package? Once have that, finish following  http://gis.stackexchange.com/questions/100010/writing-shapefile-with-projection-defined-crashes-fiona
	'''
	coordPairs = []
	for i in range(0,len(x)):
		coordPairs.append( (x[i],y[i]) )

	line = geometry.LineString(coordPairs)

	# Define schema
	schema = {
		'geometry': 'LineString',
		'properties': {'id': 'int'},
	}	

	# Write a new Shapefile
	with fiona.open(shpFnOut, 'w', 'ESRI Shapefile', schema, crs=from_epsg(epsgNum)) as c:
		## If there are multiple geometries, put the "for" loop here
		c.write({
			'geometry': geometry.mapping(line),
			'properties': {'id':0},
		})
		
	print "Sucessfully wrote shapefile: " + shpFnOut
	
	
# Bounding box to pixel offsets for gdal.ReadAsArray()
# from https://gist.github.com/perrygeo/5667173
def bbox_to_pixel_offsets(gt, bbox):
    originX = gt[0]
    originY = gt[3]
    pixel_width = gt[1]
    pixel_height = gt[5]
    x1 = int((bbox[0] - originX) / pixel_width)
    x2 = int((bbox[1] - originX) / pixel_width) + 1

    y1 = int((bbox[3] - originY) / pixel_height)
    y2 = int((bbox[2] - originY) / pixel_height) + 1

    xsize = x2 - x1
    ysize = y2 - y1
    return (x1, y1, xsize, ysize)

# Python-implemented zonal stats
# from https://gist.github.com/perrygeo/5667173
def zonal_stats(vector_path, raster_path, nodata_value=None, global_src_extent=False):
    rds = gdal.Open(raster_path, GA_ReadOnly)
    assert(rds)
    rb = rds.GetRasterBand(1)
    rgt = rds.GetGeoTransform()

    if nodata_value:
        nodata_value = float(nodata_value)
        rb.SetNoDataValue(nodata_value)

    vds = ogr.Open(vector_path, GA_ReadOnly)  # TODO maybe open update if we want to write stats
    assert(vds)
    vlyr = vds.GetLayer(0)

    # create an in-memory numpy array of the source raster data
    # covering the whole extent of the vector layer
    if global_src_extent:
        # use global source extent
        # useful only when disk IO or raster scanning inefficiencies are your limiting factor
        # advantage: reads raster data in one pass
        # disadvantage: large vector extents may have big memory requirements
        src_offset = bbox_to_pixel_offsets(rgt, vlyr.GetExtent())
        src_array = rb.ReadAsArray(*src_offset)

        # calculate new geotransform of the layer subset
        new_gt = (
            (rgt[0] + (src_offset[0] * rgt[1])),
            rgt[1],
            0.0,
            (rgt[3] + (src_offset[1] * rgt[5])),
            0.0,
            rgt[5]
        )

    mem_drv = ogr.GetDriverByName('Memory')
    driver = gdal.GetDriverByName('MEM')

    # Loop through vectors
    stats = []
    feat = vlyr.GetNextFeature()
    while feat is not None:

        if not global_src_extent:
            # use local source extent
            # fastest option when you have fast disks and well indexed raster (ie tiled Geotiff)
            # advantage: each feature uses the smallest raster chunk
            # disadvantage: lots of reads on the source raster
            src_offset = bbox_to_pixel_offsets(rgt, feat.geometry().GetEnvelope())
            src_array = rb.ReadAsArray(*src_offset)

            # calculate new geotransform of the feature subset
            new_gt = (
                (rgt[0] + (src_offset[0] * rgt[1])),
                rgt[1],
                0.0,
                (rgt[3] + (src_offset[1] * rgt[5])),
                0.0,
                rgt[5]
            )

        # Create a temporary vector layer in memory
        mem_ds = mem_drv.CreateDataSource('out')
        mem_layer = mem_ds.CreateLayer('poly', None, ogr.wkbPolygon)
        mem_layer.CreateFeature(feat.Clone())

        # Rasterize it
        rvds = driver.Create('', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
        rvds.SetGeoTransform(new_gt)
        gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1])
        rv_array = rvds.ReadAsArray()

        # Mask the source data array with our current feature
        # we take the logical_not to flip 0<->1 to get the correct mask effect
        # we also mask out nodata values explictly
        masked = np.ma.MaskedArray(
            src_array,
            mask=np.logical_or(
                src_array == nodata_value,
                np.logical_not(rv_array)
            )
        )

        feature_stats = {
            'min': float(masked.min()),
            'mean': float(masked.mean()),
            'max': float(masked.max()),
            'std': float(masked.std()),
            'sum': float(masked.sum()),
            'count': int(masked.count()),
            'fid': int(feat.GetFID())}

        stats.append(feature_stats)

        rvds = None
        mem_ds = None
        feat = vlyr.GetNextFeature()

    vds = None
    rds = None
    return stats


def readZonalStatsOutput(zonalStatsOut):
	'''
	# Function to put zonal stats output into an array
	# Ouputs a dictionary of statistics
	'''
	
	boxNum = len(zonalStatsOut) # number of boxes over which statistics were calculated

	# intialize
	mins = []
	maxs = []
	means = []
	counts = []
	stds = []
	fids = []

	# Read stats and append to list
	for i in range(0,boxNum):
		boxNow = zonalStatsOut[i]
		fids.append(boxNow['fid'])	
		mins.append(boxNow['min'])
		means.append(boxNow['mean'])
		maxs.append(boxNow['max'])
		counts.append(boxNow['count'])	
		stds.append(boxNow['std'])
		
#	statsOut = np.array( (fids,mins,means,maxs,stds,counts) ).transpose() # stick together to write out
	statsOut = {'fids':fids,'means':means,'mins':mins,'maxs':maxs,'stds':stds,'counts':counts}	
	return statsOut

def readVertexCoordsFromPolyline(shapefileFn):
	'''
	Read vertex coordinates of polyline
	Input = polyline shapefile full filepath
	Outout = dataOut = {'xLocal','yLocal','lat','lon','samplePtsXy','samplePtsLatLon'}
	'''
	with fiona.open(shapefileFn,'r') as c:
	
		# Initialize
		sample_pts_local=[]
		sample_pts_lon_lat=[]
		xList = []
		yList = []
		lonList = []
		latList = []		
	
		# Get coordinate reference system from original shapefile
		original = Proj(c.crs)
		destination = Proj(init='EPSG:4326') # dest is WGS84 = EPSG 4326

		f = c[0] # open feature
		
		for j in f['geometry']['coordinates']:
			x = j[0]
			y = j[1]
			xList.append(j[0])
			yList.append(j[1])
			sample_pts_local.append( (x,y) )
			lon,lat = transform(original, destination,x,y)
			lonList.append(lon)
			latList.append(lat)
			sample_pts_lon_lat.append( (lon,lat) )
		
		# Stick it together in a dict for storage
		dataOut = {'xLocal':xList,'yLocal':yList,'lat':latList,'lon':lonList,'samplePtsXy':sample_pts_local,'samplePtsLatLon':sample_pts_lon_lat}
		
		return dataOut


	
def readLineCoordsFromMultiLineShp(shapefileFn,transNum):
	'''
	# TAKEN FROM M. FAHNESTOCK'S L8_sample_frames... code. 08 sep 2016
	# Unfinished as of 08 sep 2016 wha
	# 22 sep 2016 think it's finished wha

	# Mark's notes	
	# use fiona to read in points along centerline in local projection (called PS here)
	# save the PS original points in sample_pts_PS, then reproject to 4326, push back into
	# coordinates of input Point object (f) and make a list of these (out_pts_ll) and when done save
	# them as a Shapely MultiPoint object

	'''
	with fiona.open(shapefileFn,'r') as c:
	
		# Initialize
		sample_pts_local=[]
		sample_pts_lon_lat=[]
		xList = []
		yList = []
		lonList = []
		latList = []		
	
		# Get coordinate reference system from original shapefile
		original = Proj(c.crs)
		destination = Proj(init='EPSG:4326') # dest is WGS84 = EPSG 4326

		f = c[transNum] # open specified feature (transect)
		
		for j in f['geometry']['coordinates']:
			x = j[0]
			y = j[1]
			xList.append(j[0])
			yList.append(j[1])
			sample_pts_local.append( (x,y) )
			lon,lat = transform(original, destination,x,y)
			lonList.append(lon)
			latList.append(lat)
			sample_pts_lon_lat.append( (lon,lat) )
		
		# Stick it together in a dict for storage
		dataOut = {'xLocal':xList,'yLocal':yList,'lat':latList,'lon':lonList,'samplePtsXy':sample_pts_local,'samplePtsLatLon':sample_pts_lon_lat}
		
		return dataOut

# Make a shapely geometry from a list = ( (x,y) ) of coordinates
def makeGeomFromLineCoords(xyList):
	line = geometry.LineString(xyList)
	return line
	


def getPathRowList(profGeom,wrs2Fn):	
	'''
	Function to search WRS shapefile to see which landsat path/rows should be sampled for a certain geometry
	Modified from Mark Fahnestock code
	Input: profGeom = shapely geometry containing object for sampling, 
	Input: wrs2Fn = full filepath for wrs2 (worldwide reference system) shapefile that shows image locations
	Output: pr_list = list of tuples containing (path, row) of intersecting image footprints
	'''
	pr_list=[] # initialize for later storage
	wrs2_fc = fiona.open(wrs2Fn) # open wrs shapefile in fiona
	for f in wrs2_fc: # iterate over path/row features
		shape=geometry.asShape( f['geometry'] ) # make feature into shapely geometry
		if shape.intersects(profGeom): # true if path/row overlaps profile geometry to be sampled
			#print f['properties']['PATH'],f['properties']['ROW'] # print overlapping path/rows
			pr_list.append((f['properties']['PATH'],f['properties']['ROW'])) # add to list
	return pr_list


def makePolygonShpFromCoords(x,y,crs,shpFnOut):
	'''
	# Function uses shapely and fiona to make a polygon form coords
	# Inputs = x coordinate list, y coordinate list, crs = coordinate reference wkt, schema = dict defining shapefile properties
	# Modified from Mike T answer on http://gis.stackexchange.com/questions/52705/how-to-write-shapely-geometries-to-shapefiles
	# As of 08 sep 2016 doesn't encode crs because missing fiona.crs package? Once have that, finish following  http://gis.stackexchange.com/questions/100010/writing-shapefile-with-projection-defined-crashes-fiona
	'''
	coordPairs = []
	for i in range(0,len(x)):
		coordPairs.append( (x[i],y[i]) )

	poly = geometry.Polygon(coordPairs)

	# Define schema
	schema = {
		'geometry': 'Polygon',
		'properties': {'id': 'int'},
	}	

	# Write a new Shapefile
	with fiona.open(shpFnOut, 'w', 'ESRI Shapefile', schema, crs=crs) as c:
		## If there are multiple geometries, put the "for" loop here
		c.write({
			'geometry': geometry.mapping(poly),
			'properties': {'id':0},
		})
		
	print "Sucessfully wrote shapefile: " + shpFnOut


def makeRasterBoundingBoxIntoShapefile(rasterFn,shapefileOutFn):
	'''
	Gives bounding box but not correct way. Accounts for dark space around image.
	'''
	rast = rasterio.open(rasterFn)
	crs = Proj(rast.get_crs()).srs
	bounds = rast.bounds
	ul = ( bounds.left, bounds.top)
	ur = ( bounds.right, bounds.top)
	lr = ( bounds.right, bounds.bottom)
	ll = ( bounds.left, bounds.bottom)
	
	msk = rast.read_masks(1)
	i,j = np.where(msk==255)
	mostRight = np.where(j==np.max(j))
	x,y = zip(ul,ur,lr,ll)
	
	makePolygonShpFromCoords(x,y,crs,shapefileOutFn)
	
	
def doyToDate(year,doy):
	'''
	Simple function to turn a year & day of year into a date
	Input: year = year (int), doy = doy (int)
	Output: dateObject = date (datetime date object)
	'''
	dateObject = dt.date.fromordinal(dt.date(year-1,12,31).toordinal()+doy)
	return dateObject


def sinceEpoch(date): # returns seconds since epoch
	'''
	marks function to get seconds in a year for use in toYearFraction()
	'''
	return time.mktime(date.timetuple())

def toYearFraction(date):
	'''
	Turns a datetime date object into a fractional year date (float)
	Input: date (datetime date object)
	Output: fractional year date (float)
	Written by Mark Fahnestock
	
	from Mark:
	########################################################
	# following from accepted answer here: http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years
	# but I import datetime package as dt, not datetime.datetime as dt, so this is modfied to reflect that
	########################################################
	'''
	s = sinceEpoch

	year = date.year
	startOfThisYear = dt.datetime(year=year, month=1, day=1)
	startOfNextYear = dt.datetime(year=year+1, month=1, day=1)

	yearElapsed = s(date) - s(startOfThisYear)
	yearDuration = s(startOfNextYear) - s(startOfThisYear)
	fraction = yearElapsed/yearDuration

	return date.year + fraction
 

def parseL8filename(l8fn):
	'''
	Parse Landsat 8 filename
	Input : landsat 8 filename (str)
	Output : data.viewkeys() = ['path','row','year','doy']
	'''
	path = l8fn[3:6]
	row = l8fn[6:9]
	year = l8fn[9:13]
	doy = l8fn[13:16]
	
	return {'path':path,'row':row,'year':year,'doy':doy}


def getCorrelationInfoFromFn(corrFnFullFilepath):
	'''
	Get information on image cross-correlation from correlation filename. To be used for data storage.
	Input: correlation full filepath (str) expecting format: 'LC80640172013342LGN0_LC80640172014185LGN0_2014-081_208_24_40_10_vv_filt_corr0.15_delcorr0.05.tif'
	Output: profile_dict.viewkeys() = ['delts_days','series','start_dec_year','stop_dec_year','mid_dec_year','path','row','infile','indir'
	'''
	
	filepath,part,filename=corrFnFullFilepath.rpartition('/') # parse into filename and directory
	
	# Parsing filename
	try: # if filtered
		image1,image2,yearDoyMedian,daysBtwn,sourceSize,targetSize,step,vType,filtStatus,corr,delcorr,spatialReference=filename.split('_')	
	except ValueError: # if not filtered, or if just wrong size..
		image1,image2,yearDoyMedian,daysBtwn,sourceSize,targetSize,step,vType=filename.split('_')[0:8]
		filtStatus = None
		corr = None
		delcorr = None

	series_tag=None # this just to be consistent with Mark's old scheme			
	doy1, year1 = int(parseL8filename(image1)['doy']) ,int(parseL8filename(image1)['year'])
	doy2, year2 = int(parseL8filename(image2)['doy']) ,int(parseL8filename(image2)['year'])
	path,row = int(parseL8filename(image1)['path']) ,int(parseL8filename(image1)['row'])			
# 	doy1, year1 = int(wgt.parseL8filename(image1)['doy']) ,int(wgt.parseL8filename(image1)['year'])
# 	doy2, year2 = int(wgt.parseL8filename(image2)['doy']) ,int(wgt.parseL8filename(image2)['year'])
# 	path,row = int(wgt.parseL8filename(image1)['path']) ,int(wgt.parseL8filename(image1)['row'])			


	dt_1=dt.datetime.fromordinal(dt.date(year1-1,12,31).toordinal()+doy1)
	dt_2=dt.datetime.fromordinal(dt.date(year2-1,12,31).toordinal()+doy2)
	timedelta=dt_2-dt_1
	delts_days=timedelta.days
	mid_dt=dt_1+(timedelta/2)
	start_dec_year=toYearFraction(dt_1.date())
	stop_dec_year=toYearFraction(dt_2.date())
	mid_dec_year=toYearFraction(mid_dt.date())
# 	start_dec_year=wgt.toYearFraction(dt_1.date())
# 	stop_dec_year=wgt.toYearFraction(dt_2.date())
# 	mid_dec_year=wgt.toYearFraction(mid_dt.date())	

	profile_dict={ 'delts_days':delts_days,'image1fn':image1,'image2fn':image2,\
					'series':series_tag, 'start_dec_year':start_dec_year, 'stop_dec_year':stop_dec_year, 'mid_dec_year':mid_dec_year,\
					'path':path, 'row':row, 'infile':filename, 'indir':filepath }
					
	return profile_dict



def listImagesToSample(pr_list,inDir):
	'''
	Get list of dicts containing image filenames to sample
	Input: pr_list = List of tuples containing  (path,row) of intersecting wrs footprints; inDir = directory containing images to sample
	Output: all_frames = list of dicts where dict.viewkeys() = ['infiles','path','row'] where 'infiles' is list of full path filenames
	'''
	total_files=0 # start counter
	all_frames=[] # initialize for storage
	for path,row in pr_list: # loop over scenes
		infiles=glob.glob(inDir+'LC8%03d%03d*vv**.tif'%(path,row)) # find images derived from this scene
		print "p:"+str(path) + ",r:"+str(row)+' number of frames: %d'%(len(infiles)) # say how many there are
		frame_info={'path':path, 'row':row, 'infiles':infiles} # save info for this scene
		all_frames.append(frame_info) # store info for this scene
		total_files+=len(infiles) # update counter

	print "Total files to sample: " + str(total_files)

	return all_frames


def convertProfilesListToPointTimeseries(profile_dict):
	'''
	#convert from profiles on individual days to time series at each point
	written by mark fahnestock
	'''
	delts_days=[]
	indir=[]
	infile=[]
	start_dec_year=[]
	mid_dec_year=[]
	stop_dec_year=[]
	path=[]
	row=[]
	series=[]
#	pt_speed=[[] for a in range(len(sample_pts_PS))]  # create a list of blank lists to append the speeds of each point into...
	# changed wha 26sep2016 b/c sample_pts_PS defined outside. Below line should be same length
	pt_speed=[[] for a in range(len(profile_dict[0]['speed']))]  # create a list of blank lists to append the speeds of each point into...
	for prof in profile_dict:
		delts_days.append(prof['delts_days'])
		indir.append(prof['indir'])
		infile.append(prof['infile'])
		start_dec_year.append(prof['start_dec_year'])
		mid_dec_year.append(prof['mid_dec_year'])
		stop_dec_year.append(prof['stop_dec_year'])
		path.append(prof['path'])
		row.append(prof['row'])
		series.append(prof['series'])
		# iterate over points, append timeseries of velocity at that point
		for indx,sp in enumerate(prof['speed']):
			pt_speed[indx].append(sp)

	# for storage
	pts={ 'delts_days':delts_days, 'indir':indir, 'infile':infile, 'start_dec_year':start_dec_year, 'mid_dec_year':mid_dec_year,\
			'stop_dec_year':stop_dec_year, 'path':path, 'row':row, 'series':series, 'pt_speed':pt_speed }

	return pts


def calculateSeasonalVelocityStats(uSeason):		
	'''
	Function to calculate season velocity percentiles		
	Takes in seasonal velocity (here called uSu b/c just copied from below)
	Ouputs dataOut = (x, u10, u25, u50, u75, u90)
	Something going wrong at 75th and 90th percentiles where just all nan. Think might have to do with treatment of nan
	
	Taken from old code, probably sloppy 26sep16 wha
	'''
	# initialize
	season10 = None
	season25 = None
	season50 = None
	season75 = None
	season90 = None
	
	# identify missing data and turn to nan
	nanInd = np.where(uSeason==-1)
	uSeason[nanInd] = np.nan
	
	if len(uSeason.shape)>1:
		xLen = uSeason.shape[1]
	
		for i in range(0,xLen):
			dataInd = np.isnan(uSeason[:,i]) != True
			colData = uSeason[dataInd,i]
			if len(colData) == 0:
				season10x = np.nan
				season25x = np.nan
				season50x = np.nan				
				season75x = np.nan				
				season90x = np.nan
			else:		
				season10x, season25x, season50x, season75x, season90x = np.percentile(colData,[10,25,50,75,90])
	
			if i == 0:
				season10 = season10x
				season25 = season25x
				season50 = season50x						
				season75 = season75x
				season90 = season90x
			else:
				season10 = np.column_stack((season10,season10x))	
				season25 = np.column_stack((season25,season25x))	
				season50 = np.column_stack((season50,season50x))	
				season75 = np.column_stack((season75,season75x))	
				season90 = np.column_stack((season90,season90x))																		
	else: # if not enough observations then nan
		season10 = np.empty((1,len(uSeason)))
		season10.fill(np.nan)
		season25 = np.copy(season10)
		season50 = np.copy(season10)		
		season75 = np.copy(season10)
		season90 = np.copy(season10)
		
	# Append data
	dataOut = np.row_stack((season10,season25,season50,season75,season90))
	
	
	return dataOut

def getNumDataFromJson(jsonData,startSeasons=[304, 81, 121, 265]):
	'''
	Function that records number of observations within each season
	Inputs: jsonData = json file with velocity data; startSeasons = [winter, spring, summer, fall] start season doy
	Outputs: dataTogether = [x, nWinter, nSpring, nSummer, nFall]
	# start days of each season, these based of local met data (304 = Oct 28; 81 = Mar 22; 121 = May 1; 265 = Sep 22)
	From old script, probably bad and clunky 26sep16 wha
	'''
	shapefileName = jsonData['profile_shapefile']
	#transName = shapefileName.split('_')[0]
	transName = jsonData['transect_name']
	x = np.array(jsonData['sample_pts_frontdist']) # distance from headwall [m]
	profiles = jsonData['profiles'] # getting all profiles
	numProfiles=len(profiles)

	dataNumAll = np.zeros((len(x)))
	dataNumWi = np.copy(dataNumAll)
	dataNumSp = np.copy(dataNumAll)
	dataNumSu = np.copy(dataNumAll)
	dataNumFa = np.copy(dataNumAll)

	for i in range(0,numProfiles):
		profNow = profiles[i]
		startDecYear = profNow['start_dec_year']
		midDecYear = profNow['mid_dec_year']
		endDecYear = profNow['stop_dec_year']
	
		middleDoy = (midDecYear - np.floor(midDecYear))*365.25
	
		u = np.array(profNow['speed'])
		dataInd=np.where(u!=-1)
		dataNumAll[dataInd] = dataNumAll[dataInd] + 1
	
		# Categorize by season
		if int(middleDoy) >= startSeasons[0] or int(middleDoy) < startSeasons[1]:
			seasonColor = 'c' # winter color
			midSeason = 'winter'
			dataNumWi[dataInd] = dataNumWi[dataInd] + 1
		elif int(middleDoy) >= startSeasons[1] and int(middleDoy) < startSeasons[2]:
			seasonColor = 'y' # spring color
			midSeason = 'spring'
			dataNumSp[dataInd] = dataNumSp[dataInd] + 1
		elif int(middleDoy) >= startSeasons[2] and int(middleDoy) < startSeasons[3]:
			seasonColor = 'g' # summer color
			midSeason = 'summer'
			dataNumSu[dataInd] = dataNumSu[dataInd] + 1
		elif int(middleDoy) >= startSeasons[3] and int(middleDoy) < startSeasons[0]:
			seasonColor = 'r' # fall color
			midSeason = 'fall'
			dataNumFa[dataInd] = dataNumFa[dataInd] + 1
	
	dataTogether = np.array((x,dataNumAll,dataNumWi,dataNumSp,dataNumSu,dataNumFa))
	return dataTogether.transpose()


def computeStatsAndPlot(fnIn,pathOut,startSeasons=[304, 81, 121, 265],axis='auto'):
	'''
	Script to read JSON velocity file, perform seasonal statistics, and plot results
	Input: json filename to process
	Output: data = x, nAll, uAll, nWinter, uWinter ... 
	Old and clunky
	
	Added NDSI 05oct2016 wha
	'''

	with open(fnIn) as data_file:
		data = json.load(data_file)

	# Read in JSON data to see how many observations at each x
	obsNum = getNumDataFromJson(data)

	shapefileName = data['profile_shapefile']
	transName = data['transect_name']
	x = np.array(data['sample_pts_frontdist']) # distance from headwall [m]
	profiles = data['profiles'] # getting all profiles
	numProfiles=len(profiles)

	uMatx = np.empty(len(x))
	uWi = []
	uSp = []
	uSu = []
	uFa = []

	for i in range(0,numProfiles): # iterate over profiles
		u = np.array(profiles[i]['speed'])
		if i == 0:
			uMatx = u
		else:
			uMatx = np.column_stack((uMatx,u))
		startDecYear = np.array(profiles[i]['start_dec_year'])
		startYear = np.floor(startDecYear)
		startDoy=np.round((startDecYear-startYear)*365)
		middleDecYear = np.array(profiles[i]['mid_dec_year'])
		middleYear = np.floor(middleDecYear)
		middleDoy = np.round((middleDecYear-middleYear)*365)
		endDecYear = np.array(profiles[i]['stop_dec_year'])
		endYear = np.floor(endDecYear)
		endDoy=np.round((endDecYear-endYear)*365)
		daysBtwn = np.array(profiles[i]['delts_days'])
	
		ndsi1 = np.array(profiles[i]['image1_snowIndex'])
		ndsi2 = np.array(profiles[i]['image2_snowIndex'])
		
		# Plot point
		#Start season based
		# Identify start seasons

		if int(middleDoy) >= startSeasons[0] or int(middleDoy) < startSeasons[1]:
			seasonColor = 'c' # winter color
			midSeason = 'winter'
		elif int(middleDoy) >= startSeasons[1] and int(middleDoy) < startSeasons[2]:
			seasonColor = 'y' # spring color
			midSeason = 'spring'
		elif int(middleDoy) >= startSeasons[2] and int(middleDoy) < startSeasons[3]:
			seasonColor = 'g' # summer color
			midSeason = 'summer'
		elif int(middleDoy) >= startSeasons[3] and int(middleDoy) < startSeasons[0]:
			seasonColor = 'r' # fall color
			midSeason = 'fall'
		
		# Compiling seasonal velocity
		if midSeason == 'winter':
			if len(uWi) == 0:
				uWi = u
			else:
				uWi = np.row_stack((uWi,u))
		elif midSeason == 'spring':
			if len(uSp) == 0:
				uSp = u
			else:
				uSp = np.row_stack((uSp,u))				
		elif midSeason == 'summer':
			if len(uSu) == 0:
				uSu = u
			else:
				uSu = np.row_stack((uSu,u))			
		elif midSeason == 'fall':
			if len(uFa) == 0:
				uFa = u
			else:
				uFa = np.row_stack((uFa,u))	
			

		#ax1.plot(x/1e3,u,linewidth=0,marker='.',color=seasonColor,alpha=0.5)

	# If no correlations in that season, then nan
	if len(uWi) == 0:
		uWi = np.empty((1,len(x)))
		uWiStats = np.empty((5,len(x)))
		uWi.fill(np.nan)
		uWiStats.fill(np.nan)
	else: #otherwise calculate stats
		uWiStats = calculateSeasonalVelocityStats(uWi)
	if len(uSp) == 0:
		uSp = np.empty((1,len(x)))
		uSpStats = np.empty((5,len(x)))		
		uSp.fill(np.nan)
		uSpStats.fill(np.nan)
	else:
		uSpStats = calculateSeasonalVelocityStats(uSp)			
	if len(uSu) == 0:
		uSu = np.empty((1,len(x)))
		uSuStats = np.empty((5,len(x)))			
		uSu.fill(np.nan)
		uSuStats.fill(np.nan)
	else:
		uSuStats = calculateSeasonalVelocityStats(uSu)		
	if len(uFa) == 0:
		uFa = np.empty((1,len(x)))
		uFaStats = np.empty((5,len(x)))	
		uFa.fill(np.nan)
		uFaStats.fill(np.nan)
	else:
		uFaStats = calculateSeasonalVelocityStats(uFa)			
						
	uAllStats = calculateSeasonalVelocityStats(uMatx.transpose())






	if len(uWiStats.shape)<=1: # this true if not enough data to properly calculate	
		uWiStats = np.empty((5,len(x)))
		uWiStats[:] = np.nan
	if len(uSpStats.shape)<=1: # this true if not enough data to properly calculate	
		uSpStats = np.empty((5,len(x)))
		uSpStats[:] = np.nan
	if len(uSuStats.shape)<=1: # this true if not enough data to properly calculate	
		uSuStats = np.empty((5,len(x)))
		uSuStats[:] = np.nan	
	if len(uFaStats.shape)<=1: # this true if not enough data to properly calculate	
		uFaStats = np.empty((5,len(x)))
		uFaStats[:] = np.nan			
	
	# Difference from median			
	uWiDiff50 = uWiStats[2,:] - uAllStats[2,:]
	uSpDiff50 = uSpStats[2,:] - uAllStats[2,:]
	uSuDiff50 = uSuStats[2,:] - uAllStats[2,:]
	uFaDiff50 = uFaStats[2,:] - uAllStats[2,:]				

	# Greater than 50th/75th %

	uWiGt50 = uWiDiff50 >= 0
	uSpGt50 = uSpDiff50 >= 0
	uSuGt50 = uSuDiff50 >= 0	
	uFaGt50 = uFaDiff50 >= 0	
	
	uWiGt75 = uWiStats[2,:] >= uAllStats[3,:]
	uSpGt75 = uSpStats[2,:] >= uAllStats[3,:]	
	uSuGt75 = uSuStats[2,:] >= uAllStats[3,:]	
	uFaGt75 = uFaStats[2,:] >= uAllStats[3,:]	

	# Less than 25ths
	uWiLt25 = uWiStats[2,:] <= uAllStats[1,:]
	uSpLt25 = uSpStats[2,:] <= uAllStats[1,:]
	uSuLt25 = uSuStats[2,:] <= uAllStats[1,:]
	uFaLt25 = uFaStats[2,:] <= uAllStats[1,:]

	# Filter out seasons with too little data
	obsNumThreshold = 3
	wiBadObsInd = np.where(obsNum[:,2]<obsNumThreshold)
	uWi50plot = np.copy(uWiStats[2,:])
	uWi50plot[wiBadObsInd] = np.nan

	spBadObsInd = np.where(obsNum[:,3]<obsNumThreshold)
	uSp50plot = np.copy(uSpStats[2,:])
	uSp50plot[spBadObsInd] = np.nan	

	suBadObsInd = np.where(obsNum[:,4]<obsNumThreshold)
	uSu50plot = np.copy(uSuStats[2,:])
	uSu50plot[suBadObsInd] = np.nan		

	faBadObsInd = np.where(obsNum[:,5]<obsNumThreshold)
	uFa50plot = np.copy(uFaStats[2,:])
	uFa50plot[faBadObsInd] = np.nan		

	# Set up plot
	gs = gridspec.GridSpec(4,1,height_ratios=[1,3,1,1])
	ax1=plt.subplot(gs[1])
	ax2=plt.subplot(gs[0])
	ax3=plt.subplot(gs[3])
	ax4=plt.subplot(gs[2])

	if 1:
		# all observations
		ax1.plot(x/1e3,uWi.transpose(),marker='.',color='c',lw=0,alpha=0.25)
		ax1.plot(x/1e3,uSp.transpose(),marker='.',color='y',lw=0,alpha=0.25)
		ax1.plot(x/1e3,uSu.transpose(),marker='.',color='g',lw=0,alpha=0.25)
		ax1.plot(x/1e3,uFa.transpose(),marker='.',color='r',lw=0,alpha=0.25)

	# all season IQR and median
	ax1.plot(x/1e3,uAllStats[1,:],lw=1,color='k')
	ax1.plot(x/1e3,uAllStats[2,:],lw=2,color='k',label='Median: all data')	
	ax1.plot(x/1e3,uAllStats[3,:],lw=1,color='k',label='25th-75th %: all data')
	ax1.fill_between(x/1e3,	uAllStats[1,:],uAllStats[3,:],color='gray',alpha=0.1)

	# seasonal medians
	if 1:
		ax1.plot(x/1e3,uWiStats[2,:],lw=2,color='c')
		ax1.plot(x/1e3,uSpStats[2,:],lw=2,color='y')
		ax1.plot(x/1e3,uSuStats[2,:],lw=2,color='g')	
		ax1.plot(x/1e3,uFaStats[2,:],lw=2,color='r')	
	if 0:
		ax1.plot(x/1e3,uWi50plot,lw=2,color='c')
		ax1.plot(x/1e3,uSp50plot,lw=2,color='orange')
		ax1.plot(x/1e3,uSu50plot,lw=2,color='g')
		ax1.plot(x/1e3,uFa50plot,lw=2,color='y')		

	# Seasonal difference from all data median
	ax4.plot(x/1e3,uWiDiff50,lw=2,ls='-',color='c')
	ax4.plot(x/1e3,uSpDiff50,lw=2,ls='-',color='y')
	ax4.plot(x/1e3,uSuDiff50,lw=2,ls='-',color='g')
	ax4.plot(x/1e3,uFaDiff50,lw=2,ls='-',color='r')			
	ax4.plot(x/1e3,np.zeros(len(x)),lw=1,ls='--',color='k')

	xLims4 = [-.3, .3]
	xMax4 = xLims4[1]
	xMin4 = xLims4[0]
	trueLoc1 = 0.02
	trueLoc2 = 0.04
	trueLoc3 = 0.06
	trueLoc4 = 0.08
		
	# Where seasonal velocity is gt 75th percentile
	ax4.plot(x[uWiGt75==True]/1e3,uWiGt75[uWiGt75==True]*(xMax4 - trueLoc1),lw=0,marker='o',color='c',alpha=0.5)
	ax4.plot(x[uSpGt75==True]/1e3,uSpGt75[uSpGt75==True]*(xMax4 - trueLoc2),lw=0,marker='o',color='y',alpha=0.5)
	ax4.plot(x[uSuGt75==True]/1e3,uSuGt75[uSuGt75==True]*(xMax4 - trueLoc3),lw=0,marker='o',color='g',alpha=0.5)
	ax4.plot(x[uFaGt75==True]/1e3,uFaGt75[uFaGt75==True]*(xMax4 - trueLoc4),lw=0,marker='o',color='r',alpha=0.5)		

	# Where is velocity lt 25th percentile
	ax4.plot(x[uWiLt25==True]/1e3,uWiLt25[uWiLt25==True]*(xMin4 + trueLoc1),lw=0,marker='o',color='c',alpha=0.5)
	ax4.plot(x[uSpLt25==True]/1e3,uSpLt25[uSpLt25==True]*(xMin4 + trueLoc2),lw=0,marker='o',color='y',alpha=0.5)
	ax4.plot(x[uSuLt25==True]/1e3,uSuLt25[uSuLt25==True]*(xMin4 + trueLoc3),lw=0,marker='o',color='g',alpha=0.5)
	ax4.plot(x[uFaLt25==True]/1e3,uFaLt25[uFaLt25==True]*(xMin4 + trueLoc4),lw=0,marker='o',color='r',alpha=0.5)

	saveData = np.vstack((x,obsNum[:,1],uAllStats[1,:],uAllStats[2,:],uAllStats[3,:],obsNum[:,2],uWiStats[2,:],obsNum[:,3],uSpStats[2,:],obsNum[:,4],uSuStats[2,:],obsNum[:,5],uFaStats[2,:]))
	saveDataT = saveData.transpose()
	np.savetxt(pathOut + transName +'_velocityStats.csv',saveDataT,delimiter=',')
		 
					
	## GET ELEVATION PROFILE
	#elevData = np.genfromtxt(elevationProfileFolder + shapefileName[:-4] + '_elevationProfile.csv',delimiter=',')
	if(0):
		elevData = np.genfromtxt('/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/dtm/elevationProfiles/' + transName + '_elevationProfile.csv',delimiter=',')
		ax3.plot(elevData[:,0]/1e3,elevData[:,1],color='gray')		

	ax3.plot(x/1e3,ndsi1,lw=2,alpha=0.7,c='b')
	ax3.plot(x/1e3,ndsi2,lw=2,alpha=0.7,c='b')	

	ax2.plot(obsNum[:,0]/1e3,obsNum[:,1],marker='.',color='k')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,2],marker='.',color='c')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,3],marker='.',color='y')	
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,4],marker='.',color='g')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,5],marker='.',color='r')		
				


	ax3.set_xlabel('Distance from headwall [km]',fontsize=16)
	ax1.set_ylabel('Velocity [m d$^{-1}$]',fontsize=16)
	ax2.set_title(transName,fontsize=18)

	ax1.plot(-1,-1,marker='.',color='c',linewidth=0,label='Winter')
	ax1.plot(-1,-1,marker='.',color='y',linewidth=0,label='Spring')
	ax1.plot(-1,-1,marker='.',color='g',linewidth=0,label='Summer')
	ax1.plot(-1,-1,marker='.',color='r',linewidth=0,label='Fall')


	#ax3.set_ylim([np.nanmin(elevData[:,1]),	np.nanmax(elevData[:,1])])
	ax2.set_ylabel('# Observations',fontsize=18)
	ax2.yaxis.set_label_position("right")
	ax3.set_ylabel('NDSI [-]',fontsize=18)
	#ax2.set_yticks(np.arange( round(elevData[:,1].min()/500)*500,round(elevData[:,1].max()/500)*500+500,500))
	ax3.set_yticks((np.arange(-0.25,1.25,0.25)))
	#ax1.set_ylim([0,np.nanmin((2.5,np.nanmax(uAllStats[4,:])))]) # auto-scaling limits
	ax1.set_ylim([0,2])
	ax1.set_xlim([0,x.max()/1e3])
	ax1lims = ax1.xaxis.get_view_interval()
	ax3.set_xlim([ax1lims[0],ax1lims[1]])
	ax2.set_xlim([ax1lims[0],ax1lims[1]])
	ax4.set_xlim([ax1lims[0],ax1lims[1]])
	ax4.set_ylim(xLims4)
	ax4.set_ylabel('Diff. from annual median',fontsize=18)
	ax4.yaxis.set_label_position("right")
	
	ax1.legend(loc='best',prop={'size':6},numpoints=1)
	plt.draw()
	plt.savefig(pathOut + transName + '_profilePlot_wSeasonalMedian_wElevation_wDataNum_wDiffFromMean.pdf')
	#plt.show()
	plt.close()		
		
	return saveData


def makeRGBfromBands(folderPath,outName,tr = 'auto', bs=[4,3,2]):
	'''
	Make color composite from three separate bands
	Inputs: folderPath = path containing image bands, tr = 100 100 = target resolution for 100 m pixels, bs = (4,3,2) = bands to use for R,G,B
	Outputs: tiff at outName
	'''
	# Remove output files if they exist 
	if os.path.isfile(outName):
		print "Output filename already exists; deleting"
		os.remove(outName)
	
	if os.path.isfile('tmp.tif'):
		os.remove(outName)
	
	# Get filenames for bands to be used for R,G,B
	rFn = glob.glob(folderPath + '*B' + str(bs[0]) + '.TIF' or folderPath + '*B' + str(bs[0]) + '.tif')[0]
	gFn = glob.glob(folderPath + '*B' + str(bs[1]) + '.TIF' or folderPath + '*B' + str(bs[1]) + '.tif')[0]
	bFn = glob.glob(folderPath + '*B' + str(bs[2]) + '.TIF' or folderPath + '*B' + str(bs[2]) + '.tif')[0]

	print "Red : " + rFn
	print "Green : " + gFn
	print "Blue : " + bFn
	
	if tr == 'auto': # assume output same resolution of not specified
		mergeCommandString = "gdal_merge.py -separate -co " + '"' + 'COMPRESS=LZW' + '"' + " -o " + "'" + outName  + "'" + " " + rFn + " " + gFn + " " + bFn
	else: # if target resolution specified
		mergeCommandString = "gdal_merge.py -separate -co " + '"' + 'COMPRESS=LZW' + '"' + " -ps " + str(tr) + " " + str(tr) + " -o " + "'" + 'tmp.tif'  + "'" + " " + rFn + " " + gFn + " " + bFn
	
	# run gdal_merge.py to make multiband
	os.system(mergeCommandString)
	
	# run gdal_translate to make RGB interpretation
	translateCommandString = "gdal_translate -of 'GTiff' -co " + '"' + 'PHOTOMETRIC=RGB' + '" -scale 0 100 -a_nodata 0 tmp.tif ' + outName
	os.system(translateCommandString)
	
	print "Created file: " + outName


def interpolateShapelyLineGeom(line,dDist):
	'''
	Interpolate points at specified distance along shapely line geometry
	Inputs: line = shapely line geometry, dDist = specified interpolation distance
	Outputs: x, y coordinates at specified interpolation distance
	'''

	lineLen = line.length # get line length

	interpDist = np.arange(0,lineLen,dDist) # make an array of desired distances to get coordiantes for

	# initialize
	xInt = []
	yInt = []

	# iterate over interpolation list and get line coords at specified distances
	for distNow in interpDist:
		xNow,yNow = line.interpolate(distNow).xy
		xInt.append(xNow[0])
		yInt.append(yNow[0])
	
	return xInt, yInt

def calculateNDSI(imageDirectory,saveImageToggle=0,overwrite=0):
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

		if 0: # disabling functionality for now 05oct2016 wha
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

def getUtmZones(shapelyGeometry):
	'''
	Get UTM zone numbers intersecting a geometry
	Input: shapely geometry
	Output: list of utm zone numbers
	'''
	
	utmZoneShpFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/qgis/UTM_Zone_Boundaries/UTM_Zone_Boundaries.shp'
	utmZones = fiona.open(utmZoneShpFn)

	zones = []
	for zone in utmZones:
		zoneGeom = geometry.asShape(zone['geometry'])
		if zoneGeom.intersects(shapelyGeometry):
			zones.append(zone['properties']['ZONE'])

	return zones

def getEPSGforUTMzone(utmZone):
	'''
	Very simple fn to get EPSG number for a given utm zone (assumes northern hemisphere)
	'''
	return 32600 + int(utmZone)

def yearFracToDoy(yearFrac):
	'''
	Turn a year fraction to a day of year. Could be off by one in leap years..
	Input = fractional year date (e.g., 2015.328)
	Output = day of year of input (DOY 120 for above yearFrac)
	'''
	return int(np.round((yearFrac - int(yearFrac))*365))



def flipDistance(dist):
	'''
	Change distance from headwall to distance from terminus
	'''
	
	totalDist = np.max(dist)
	termDist = np.ones(len(dist))*totalDist - np.array(dist)
	
	return termDist

def normalizeDistance(dist):
	'''
	Normalize distance by total glacier length
	'''
	totalDist = np.max(dist)
	normDist = np.ones(len(dist))*dist / totalDist
	
	return normDist
	








