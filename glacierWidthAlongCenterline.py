## SCRIPT TO FIND THE CROSS-GLACIER WIDTH ALONG KIENHOLZ CENTERLINE PROFILES
#  01 July 2016
#  William Armstrong

# IMPORT MODULES
import numpy as np
from osgeo import ogr, gdal, osr
import matplotlib.pyplot as plt
import sys
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection, PolyCollection
import time

# FUNCTIONS
# Calculate the slope of a line
# Stolen from krb_vectools
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

# Calculate the slope and perpendicular of a line
# Stolen from krb_vectools
def calcSlopeAndPerpSlope(ptArray, hwin, units="radians"):

    slope=[]
    perpSlope=[]
    
    for i in range(hwin): 
        firstPt=ptArray[0]
        lastPt=ptArray[i+hwin]
        slope.append(calculateSlope(firstPt, lastPt,units))
        perpSlope.append(((slope[i]+2*np.pi)-(np.pi/2))%(2*np.pi))
        
    for i in range(hwin,len(ptArray)-hwin):
        firstPt=ptArray[i-hwin]
        lastPt=ptArray[i+hwin]
        slope.append(calculateSlope(firstPt, lastPt,units))
        perpSlope.append(((slope[i]+2*np.pi)-(np.pi/2))%(2*np.pi))
        
    for i in range(len(ptArray)-hwin, len(ptArray)): 
        firstPt=ptArray[i-hwin]
        lastPt=ptArray[-1]
        slope.append(calculateSlope(firstPt, lastPt,units))
        perpSlope.append(((slope[i]+2*np.pi)-(np.pi/2))%(2*np.pi))    
    return slope, perpSlope

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

# Calculate distahce along a line
# Stolen from krb_vectools
def distanceAlongLine(x,y):   
    xOut=np.asarray(x)
    yOut=np.asarray(y)    
    segments=(((xOut[1:]-xOut[0:-1])**2)+(yOut[1:]-yOut[0:-1])**2)**(0.5)
    dist=np.cumsum(segments)
    dist=np.insert(dist, 0, 0)
    return dist
 
# Function to transform ogr geometry into different spatial reference system
def transformGeometry(geom,desiredSrsEpsg):
	# Get spatial reference system of original geometry
	sourceSrs = osr.SpatialReference()
	geomSrsProj4 = geom.GetSpatialReference().ExportToProj4()
	sourceSrs.ImportFromProj4(geomSrsProj4) #

	desiredSrs = osr.SpatialReference()
	desiredSrs.ImportFromEPSG(desiredSrsEpsg) # Creating new spatial reference

	srsTfm = osr.CoordinateTransformation(sourceSrs, desiredSrs) # Create coordinate transformation

	# Convert geometry to desired srs
	transformedGeom = geom.Clone() # make a copy
	transformedGeom.Transform(srsTfm) # transform
	
	return transformedGeom

# Calculate UTM zone for a geometry
def getUtmEpsgFromGeometry(geom):
	tfmGeom = transformGeometry(geom,4326) # convert to wgs 84
	
	geomCentroid = tfmGeom.Centroid() # Get geometry centroid
	
	# Find utm zone
	geomCentroidLon = geomCentroid.GetX()	
	utmZone = np.floor((geomCentroidLon + 180) / 6.).astype('int') + 1
	epsgForUtmZone1 = 32601 - 1 # EPSG code for WGS 84 UTM Zone 1 N, -1 b/c going to do + zone
	epsgForUtmZone = epsgForUtmZone1 + utmZone # EPSG for desired UTM zone
	
	return epsgForUtmZone

# This not especially useful b/c GLIMS id's in centerlines shapefile is different from RGI for same glaciers..
def returnGlacierInd(rgiLyr,glimsIdDesired):
	glId = None
	nGlaciers = rgiLyr.GetFeatureCount()
	for j in range(0,nGlaciers):
		gl = rgiLyr.GetFeature(j)
		glimsIdNow = gl.GetField(1)
		
		if glimsIdNow == glimsIdDesired:
			glId = j # index for glacier in which centerline resides
	
	# True if no match
	if glId is None:
		glId = 'nan'
		
	return glId

# Function to calculate glacier width along a centerline profile. Could definitely pare this down by adding functions within.	
# Inputs: centerline geometry, outputs: width at each centerline vertex
def calculateWidthAlongCenterline(clGeom,rgiLyr):

	nPoints = clGeom.GetPointCount()
	clGeomUtm = clGeom.Clone() # Make copy for later	
	middleX, middleY, middleZ = clGeom.GetPoint(np.floor(nPoints/2).astype('int')) # coordinates of point half way down centerline profile
	middleGeom = ogr.Geometry(ogr.wkbPoint)
	middleGeom.AddPoint(middleX,middleY,0)
	# Making copy to transform to WGS later
	middleGeomWgs = middleGeom.Clone()

	# Centroid location (not esp. useful b/c curved centerlines)
	cent = clGeom.Centroid()
	centX = cent.GetX(0)
	centY = cent.GetY(0)	
	t1 = time.time()
	rgiInd = None
	glimsId = None
	middleContained = False
	# Loop over RGI to find containing glacier
	nGlaciers = rgiLyr.GetFeatureCount()
	for j in range(0,nGlaciers):
	#while middleContained == False: # tried this with a break and it actually made the process slower		
		gl = rgiLyr.GetFeature(j)
		glGeom = gl.geometry()

		middleContained = middleGeom.Intersects(glGeom)
		
		if middleContained == True:
			glimsId = gl.GetField(1)
			rgiInd = j # index for glacier in which centerline resides
			break


	if middleContained == False:
		sys.exit('Iteration: ' + str(i) + '; No glacier found to contain middle point on centerline')

	containingFeat = rgiLyr.GetFeature(rgiInd)
	containingGeom = containingFeat.geometry()

	# Find other glaciers that touch and merge them into one outline
	polyGeom = ogr.Geometry(ogr.wkbMultiPolygon)
	for j in range(0,nGlaciers):
		otherGl = rgiLyr.GetFeature(j)
		otherGlGeom = otherGl.geometry()
	
		# 1 if glaciers touch
		touchingCurrentGl = containingGeom.Intersects(otherGlGeom)
	
		if touchingCurrentGl == 1: # add geometry to multipolygon if touching
			otherGlGeomRef = otherGl.GetGeometryRef()
			polyGeom.AddGeometry(otherGlGeomRef)

	# "Union" (merge) touching polygons
	glacierUnionGeom = polyGeom.UnionCascaded()

	# Convert spatial reference to get into useful system for measurement
	sourceSrs = osr.SpatialReference()
	sourceSrs.ImportFromEPSG(3338) # Alaska Albers

	wgsSrs = osr.SpatialReference()
	wgsSrs.ImportFromEPSG(4326) # WGS

	tfmAlbersToWgs = osr.CoordinateTransformation(sourceSrs, wgsSrs) # Create coordinate transformation

	# Convert glacier geometry and current centerline middle to WGS
	containingGeomWgs = containingGeom.Clone()
	containingGeomWgs.Transform(tfmAlbersToWgs)
	middleGeomWgs.Transform(tfmAlbersToWgs)

	# Find utm zone
	middleGeomLon = middleGeomWgs.GetX()	
	utmZone = np.floor((middleGeomLon + 180) / 6.).astype('int') + 1
	epsgForUtmZone1 = 32601 - 1 # EPSG code for WGS 84 UTM Zone 1 N, -1 b/c going to do + zone
	epsgForUtmZone = epsgForUtmZone1 + utmZone # EPSG for desired UTM zone

	# Create transformation for WGS to UTM
	utmSrs = osr.SpatialReference()
	utmSrs.ImportFromEPSG(epsgForUtmZone)
	tfmWgsToUtm = osr.CoordinateTransformation(wgsSrs,utmSrs)

	# And tfm from Albers to Utm
	tfmAlbersToUtm = osr.CoordinateTransformation(sourceSrs, utmSrs)

	# Convert to UTM
	# First make copies
	middleGeomUtm = middleGeomWgs.Clone()
	containingGeomUtm = containingGeomWgs.Clone()
	containingGlacierUnionGeomUtm=glacierUnionGeom.Clone()	
	middleGeomUtm.Transform(tfmWgsToUtm)
	containingGeomUtm.Transform(tfmWgsToUtm)
	clGeomUtm.Transform(tfmAlbersToUtm)
	containingGlacierUnionGeomUtm.Transform(tfmAlbersToUtm)

	middleXUtm = middleGeomUtm.GetX()
	middleYUtm = middleGeomUtm.GetY()

	# Utm coord list
	# Centerline coordinates
	# Initialize
	clX = []
	clY = []	
	for k in range(0,nPoints):
		clX.append(clGeomUtm.GetX(k))
		clY.append(clGeomUtm.GetY(k))

	# Find angle between 90 bearing (east, really 0 rad) and centerline segments
	dE = np.diff(clX) # change in easting between points
	dN = np.diff(clY) # change in northing between points
	theta = np.arctan((dN/dE)) # angle between (1,0) and dN/dE, lying between -pi/2 and pi/2. Not caring about fact that arctan(-1/-1) is same as arctan(1/1) because drawing line in both directions
	thetaPerp = ((theta+2*np.pi)-(np.pi/2))%(2*np.pi) # angle perpendicular to theta
	thetaPerp = np.insert(thetaPerp,0,thetaPerp[0]) # assigning first point value of 2nd point (not there b/c diff)

	width = []

	for k in range(0,nPoints): # calculate width along centerline profile
		# Make verticies and center coords for a line perpendicular to flow direction 
		perpLineCoords = makeLine((clX[k],clY[k]),thetaPerp[k],15000)

		# Make transverse line a geometry
		perpLineGeom = ogr.Geometry(ogr.wkbLineString)
		perpLineGeom.AddPoint(perpLineCoords[0,0],perpLineCoords[0,1],0)	
		perpLineGeom.AddPoint(perpLineCoords[2,0],perpLineCoords[2,1],0)		

		# Find crossing points between transverse line and glacier geometry
		crossingGeom = perpLineGeom.Intersection(containingGlacierUnionGeomUtm)

		# Iterate over points where perpendicular line intersects glacier geometry
		crossX = []
		crossY = []
	
		# Not sure why there are two ways the crossing geometry can come..
		if crossingGeom.GetGeometryCount()>0:
			for n in range(0,crossingGeom.GetGeometryCount()):
				crossGeom = crossingGeom.GetGeometryRef(n)
				for m in range(0,crossGeom.GetPointCount()):
					crossX.append(crossGeom.GetX(m))
					crossY.append(crossGeom.GetY(m))
		elif crossingGeom.GetPointCount()>0:
			for n in range(0,crossingGeom.GetPointCount()):
				crossX.append(crossingGeom.GetX(n))
				crossY.append(crossingGeom.GetY(n))				
	
		dX = np.array(crossX) - clX[k]
		dY = np.array(crossY) - clY[k]
		dist = np.hypot(dX,dY)	
	
		minusInd = np.where(dX < 0)[0]
		plusInd = np.where(dX >= 0)[0]
		closestPlus = np.where(dist[plusInd] == np.min(dist[plusInd]))[0]
		closestMinus = np.where(dist[minusInd] == np.min(dist[minusInd]))[0]	
		closeX = np.array((crossX[plusInd[closestPlus]],crossX[minusInd[closestMinus]]))
		closeY = np.array((crossY[plusInd[closestPlus]],crossY[minusInd[closestMinus]]))
	
		if k == 0: # special treatment for first point
			width = np.array(dist[plusInd[closestPlus]] + dist[minusInd[closestMinus]])
		else:
			width = np.append(width,dist[plusInd[closestPlus]] + dist[minusInd[closestMinus]])

	t2 = time.time()
	print "Took: " + str(t2-t1) + " seconds to calculate width for " + str(k) + " points"


	# Just for debugging
	if 1:
		if 0: # for old, un-unioned shapefiles
			# Just for debugging, points to plot
			nGeom = containingGeomUtm.GetGeometryCount()
			bX = []
			bY = []
			for q in range(0,nGeom):
				boundNow = containingGeomUtm.GetGeometryRef(q)
				for m in range(0,boundNow.GetPointCount()):
					bX.append(boundNow.GetX(m))
					bY.append(boundNow.GetY(m))
		if 1: # for unioned shapefiles
			# Just for debugging, points to plot
			nGeom = containingGlacierUnionGeomUtm.GetGeometryCount()
			bX = []
			bY = []
			for q in range(0,nGeom):
				boundNow = containingGlacierUnionGeomUtm.GetGeometryRef(q)
				for m in range(0,boundNow.GetPointCount()):
					bX.append(boundNow.GetX(m))
					bY.append(boundNow.GetY(m))
				
			
			plt.scatter(bX,bY,lw=0,marker='.',facecolor='k')
			plt.scatter(clX,clY,c=width,s=width/50,lw=0)
			cbar = plt.colorbar()
			cbar.set_label('Glacier width [m]',fontsize=16)
			#plt.plot(clX[k],clY[k],marker='v',markersize=15)
			#plt.plot(perpLineCoords[:,0],perpLineCoords[:,1])
			#plt.plot(crossX,crossY,marker='o',lw=0,markersize=15,markerfacecolor=None)
			#plt.plot(closeX,closeY,marker='.',lw=0,markersize=10)		
			#plt.savefig('kennGlacierWidthExample.png',r=300)
			plt.show()
			plt.close()			

	return width	


## FILENAMES
#centerlineFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/qgis/Alaska_centerline_from_ChristianKleinholz/Centerlines.shp'
centerlineFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/qgis/keinholz_centerlines_paredDownToStudyArea.shp'
rgiFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/qgis/rgi_50_epsg3338.shp'

## ANALYSIS

# Open shapefiles
clShp = ogr.Open(centerlineFn)
rgiShp = ogr.Open(rgiFn)
clLyr = clShp.GetLayer(0)
rgiLyr = rgiShp.GetLayer(0)
nCl = clLyr.GetFeatureCount()

# Iterate through centerlines
for i in range(0,nCl):
		
	cl = clLyr.GetFeature(i)
	objectId = cl.GetField(0)
	clGeom = cl.geometry()

	# Get x and y coordinates along centerline
	nPointsCl = clGeom.GetPointCount()
	clX = []
	clY = []
	utmEpsg = getUtmEpsgFromGeometry(clGeom) # get epsg code for utm zone
	clGeomUtm = transformGeometry(clGeom,utmEpsg) # convert geometry to utm	
	for j in range(0,nPointsCl):
		clX.append(clGeomUtm.GetX(j))
		clY.append(clGeomUtm.GetY(j))	
	
	dist = distanceAlongLine(clX,clY)
	clWidth = calculateWidthAlongCenterline(clGeom,rgiLyr)
	
	
	
		

			
		
		
			
		
	
	
	
	
		


