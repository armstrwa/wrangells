'''
New script for swath profiler that does not require GRASS
06 Sep 2016
William Armstrong

Inputs: line geometry, raster
Outputs: raster swath profiled along line geometry

'''

## -------------- ##
## IMPORT MODULES ##
## -------------- ##

import os
import sys
import numpy as np
#from krb_vecTools import *
import scipy as sci 
from osgeo import ogr, osr, gdal
from osgeo import osr
from osgeo import gdal
from gdalconst import *
import csv
import glob
import matplotlib.pyplot as plt
import json
from shapely import geometry
import fiona
from pyproj import Proj, transform
from fiona.crs import from_epsg

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
	

# Function uses shapely and fiona to make a polyline form coords
# Inputs = x coordinate list, y coordinate list, epsgNum = epsg coordinate system code, schema = dict defining shapefile properties
# Modified from Mike T answer on http://gis.stackexchange.com/questions/52705/how-to-write-shapely-geometries-to-shapefiles
# As of 08 sep 2016 doesn't encode crs because missing fiona.crs package? Once have that, finish following  http://gis.stackexchange.com/questions/100010/writing-shapefile-with-projection-defined-crashes-fiona
def makePolylineShpFromCoords(x,y,epsgNum,shpFnOut):
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


# Function to put zonal stats output into an array
# Ouputs a dictionary of statistics
def readZonalStatsOutput(zonalStatsOut):
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
	
# TAKEN FROM M. FAHNESTOCK'S L8_sample_frames... code. 08 sep 2016
# Unfinished as of 08 sep 2016 wha

# Mark's notes	
# use fiona to read in points along centerline in local projection (called PS here)
# save the PS original points in sample_pts_PS, then reproject to 4326, push back into
# coordinates of input Point object (f) and make a list of these (out_pts_ll) and when done save
# them as a Shapely MultiPoint object

def readLineCoordsFromMultiLineShp(shapefileFn,transNum):
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
	
			
#		print f['geometry']['coordinates'][0]
	
		#x,y =  f['geometry']['coordinates']
#		sample_pts_PS.append((x,y))
# 		sample_pts_frontdist.append(f['properties']['cng_(Meter'])   #'cng_(Meter' a holdover from chainage QGIS plugin - will need to fix
#		sample_pts_frontdist.append(float(f['properties'][key_for_distance_field]))
#		if (args.key_for_pointID):
#			sample_pts_point_id.append(f['properties'][args.key_for_pointID])
#		if (args.key_for_lineID):
#			sample_pts_line_id.append(f['properties'][args.key_for_lineID])
		
		#sample_pts_lon_lat.append((lon,lat))
		# 		print (long,lat,x,y)
		# change only the coordinates of the feature
		#f['geometry']['coordinates'] = (lon,lat)
		#out_pts_ll.append(f)


	

## ----------- ##
## USER INPUTS ##
## ----------- ##



dDist = 500 # interpolation distance along transect line
boxWidth = 1000 # swath profiler width across transect line
epsgNum = 32607 # spatial reference EPSG code
		
#shapefileIn='/Users/wiar9509/git/pySwath/stEliasTransectsUtm7n.shp'		
shapefileIn='/Users/wiar9509/git/pySwath/wrangellsTransectsUtm7n.shp'
rasterFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/AK_BA_62_17/akv2ba_S8_062_017_016_2014_043_2014_059_hp_filt_3.0_vv.tif'
rasterDirectory = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/test/'
jsonDataFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/kaskWB_evenlySpace_100m_profiles_sampled_2015-12-04.json'
wrs2Fn = '/Users/wiar9509/Documents/generalScripts/fromMarkFahnestock/wrs2_descending/wrs2_descending.shp'
				
## -------- ##
## ANALYSIS ##
## -------- ##


# Bringing in a test centerline geometry
#shapefileIn='/Users/wiar9509/Documents/generalScripts/swath/shapefileInputs/newWStETransects_19nov15_fromQGIS_EGPSG102006.shp'

shp = ogr.Open(shapefileIn)
lyr = shp.GetLayer(0)

if 0: # to get names and associated index
	for i in range(0,lyr.GetFeatureCount()):
		feat = lyr.GetFeature(i)
		print i, feat.GetField(1)
	
#feat = lyr.GetFeature(59) # kaskWB if using the epsg102006 transects
#feat = lyr.GetFeature(19) # kaskWB if using the st elias utm7n transects
feat = lyr.GetFeature(8) # nabeEB if using the wrangells utm7n transects
lineGeom = feat.geometry()

# Make shapefile for swath profiling
xInt,yInt = interpolateVerticesOfPolyline(lineGeom,dDist) # interpolate vertices
slopes=calculateSlopes(xInt,yInt) # calculate "slope" (azimuth) of line segments
perpSlopes=calculatePerpendicularSlopes(slopes) # calculate perpendicular "slope" (azimuth)
allBoxVerts = makeAllBoxVertices(xInt,yInt,boxWidth) # get box coordinates for each interpolated point, with a specified box width (cross-transect coord)

# Create the spatial reference
srs = osr.SpatialReference()
srs.ImportFromEPSG(epsgNum) # This hard-coded to UTM 7N w/ WGS84 datum. Would be nice to have this defined based of input shapefile

makeMultipolygonFromBoxVertices(allBoxVerts,srs,'test.shp') # Make shapefile for swath profiling


profileCoords = readLineCoordsFromMultiLineShp(shapefileIn,14)
geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']}
profGeom = geometry.shape(geomToPass)

#profGeomLatLon = makeGeomFromLineCoords(profileCoords['samplePtsLatLon'])

#profGeom = geometry.shape(tuple(profileCoords['samplePtsLatLon']))


#points_ll = geometry.multipoint.MultiPoint([geometry.shape(f['geometry']) for f in out_pts_ll])

#front_dist=[x['properties'][key_for_distance_field] for x in out_pts_ll]

# now run through the WRS2 outlines and find path,rows that intersect with the MultiPoint object

pr_list=[]
wrs2_fc = fiona.open(wrs2Fn)
for f in wrs2_fc:
    shape=geometry.asShape( f['geometry'] )
    if shape.intersects(profGeomLatLon):
        print f['properties']['PATH'],f['properties']['ROW']
        pr_list.append((f['properties']['PATH'],f['properties']['ROW']))


total_files=0
all_frames=[]
for path,row in pr_list:
	infiles=glob.glob(indir+'/'+fn_start+'S8_%03d_%03d_*vv.tif'%(path,row))
	print path,row,' number of frames: %d'%(len(infiles))
	frame_info={'path':path, 'row':row, 'infiles':infiles}
	all_frames.append(frame_info)
	total_files+=len(infiles)





# Brute force, just to every file in directory
fileList = glob.glob(rasterDirectory+'*064017*')
for rasterFn in fileList:
	print "Profiling file: " + fileList[0].split('/')[-1]
	raster = gdal.Open(rasterFn)
	xyTfm = createTransformation(rasterFn) # transformation parameters for converting raster pixels to ground coordinates

	stats=zonal_stats('test.shp', rasterFn, nodata_value=-1)
	statsDict = readZonalStatsOutput(stats)
	dist = distanceAlongLine(xInt,yInt)
	dist = np.delete(dist,0) # delete the first entry, which we don't have a value for b/c don't have a slope for the first entry (b/c need to points), so interpolator doesn't make a 0,0 point .. this quick and dirty
	statsDict['dist'] = dist

	plt.plot(statsDict['dist'],statsDict['means'],lw=2)
	plt.plot(statsDict['dist'],statsDict['counts'],c='gray',ls='--')

plt.show()
plt.close()






		
		
	


