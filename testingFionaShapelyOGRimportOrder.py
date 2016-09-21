'''
Running tests to see what orders of module import don't cause python to crash when performing
simple geometric operations such as shape.centroid. Need to import shapely first, then fiona/OGR
20 Sep 2016
'''

# testing why keeps crashing
if 1:
	if 0: # this WORKS
		from shapely import geometry
		import fiona
		import json
	
		with open('profileCoords.json') as data_file:
			profileCoords = json.load(data_file)
		geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']}
		profGeom = geometry.shape(geomToPass)
		
		print "centoid coords:"
		print profGeom.centroid.coords.xy
		
	if 0: # this WORKS	
		from shapely import geometry
		from osgeo import ogr, osr, gdal
		import fiona
		import json
	
		with open('profileCoords.json') as data_file:
			profileCoords = json.load(data_file)
		geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']}
		profGeom = geometry.shape(geomToPass)
		
		print "centoid coords:"
		print profGeom.centroid.coords.xy

	if 0: # this WORKS	
		from shapely import geometry
		import fiona
		from osgeo import ogr, osr, gdal		
		import json
	
		with open('profileCoords.json') as data_file:
			profileCoords = json.load(data_file)
		geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']}
		profGeom = geometry.shape(geomToPass)
		
		print "centoid coords:"
		print profGeom.centroid.coords.xy
	
	if 0: # this WORKS	
		from osgeo import gdal
		from shapely import geometry
		import fiona
		import json
	
		with open('profileCoords.json') as data_file:
			profileCoords = json.load(data_file)
		geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']}
		profGeom = geometry.shape(geomToPass)
		
		print "centoid coords:"
		print profGeom.centroid.coords.xy		

	if 0: # this fails	
		from osgeo import ogr
		from shapely import geometry
		import fiona
		import json
	
		with open('profileCoords.json') as data_file:
			profileCoords = json.load(data_file)
		geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']}
		profGeom = geometry.shape(geomToPass)
		
		print "centoid coords:"
		print profGeom.centroid.coords.xy		

	if 0: # this fails	
		from osgeo import osr
		from shapely import geometry
		import fiona
		import json
	
		with open('profileCoords.json') as data_file:
			profileCoords = json.load(data_file)
		geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']}
		profGeom = geometry.shape(geomToPass)
		
		print "centoid coords:"
		print profGeom.centroid.coords.xy		
		
	if 0: # this fails
		import fiona	
		from shapely import geometry
		import json
	
		with open('profileCoords.json') as data_file:
			profileCoords = json.load(data_file)
		geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']}
		profGeom = geometry.shape(geomToPass)
		
		print "centoid coords:"
		print profGeom.centroid.coords.xy
		
	if 1: # THIS WORKS -- NEED GDALCONST TO BE AFTER SHAPELY TOO
		shapefileIn='/Users/wiar9509/git/pySwath/wrangellsTransectsUtm7n.shp'
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
				
		import os
		import sys
		import numpy as np
		#from krb_vecTools import *
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
		
		profileCoords = readLineCoordsFromMultiLineShp(shapefileIn,14)
		geomToPass = {'type':'LineString','coordinates':profileCoords['samplePtsLatLon']}
		profGeom = geometry.shape(geomToPass)
		
		print "coords:"
		print profGeom.centroid.coords.xy
		
		
		
