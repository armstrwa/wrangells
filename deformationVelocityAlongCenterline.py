## SCRIPT TO FIND THE CROSS-GLACIER WIDTH ALONG KIENHOLZ CENTERLINE PROFILES
#  01 July 2016
#  William Armstrong

# IMPORT MODULES
import numpy as np
from osgeo import ogr, gdal, osr
import matplotlib.pyplot as plt
import sys, os
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection, PolyCollection
import time
import scipy.signal as sig
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import UnivariateSpline
import json

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

# Function to sample raster at specified coordinates
# Returns values of raster at x,y
def sampleRasterAtXY(rasterInFn,x,y,xyTfm):
	
	z = []
	
	imageIj = xyToImageIj(x,y,xyTfm)
	
	rast = gdal.Open(rasterInFn)
	band = rast.GetRasterBand(1) # assumes single band raster
	bandArr = band.ReadAsArray() # convert band into array
	
	samplePtNum = imageIj.shape[0] # number of sample points
	
	for i in range(0,samplePtNum):
		z.append( bandArr[ np.round(imageIj[i,1]), np.round(imageIj[i,0] ) ] )
	
	return z
	 
 
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

# Get coordinates of line geometry vertices
def getLineGeometryVerticesXY(lineGeom):
	nPoints = lineGeom.GetPointCount()
	clX = []
	clY = []
	
	for j in range(0,nPoints):
		clX.append(lineGeom.GetX(j))
		clY.append(lineGeom.GetY(j))	
		
	return clX, clY

def reprojectRaster(rasterFnIn,epsgIn,epsgOut):

	## SCRIPT TO ITERATE RASTER REPROJECTION
	#  Written by William Armstrong
	#  25 November 2015
	#
	#  With help from: http://stackoverflow.com/questions/8715419/passing-a-variable-set-by-the-user-as-a-parameter-for-a-os-system-command
	#  and http://gis.stackexchange.com/questions/92907/re-project-raster-image-from-mercator-to-equirectangular

	outDir = os.path.dirname(rasterFnIn) + '/'
	rastFn = rasterFnIn.split('/')[-1]
	fnOut = outDir + rastFn[0:-4] + '_EPSG' + str(epsgOut) + '.tif'

	if not os.path.isfile(fnOut):
		os.system("gdalwarp -s_srs EPSG:" +  str(epsgIn) + " -t_srs EPSG:" + str(epsgOut) + " '" + rasterFnIn + "' " + " '" + fnOut + "' ")

# Stolen from KRB_vectools
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]] # This is doing a reflection at the boundaries
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    yTrim=y[window_len/2:-window_len/2+1] # this added by wha on 08jul2016. Was outputting longer array than input because of edge effects
    return yTrim

# Estimate ice thickness based on constant basal shear stress
def estimateIceThicknessYieldStress(x,z,yieldStress,sigma):
	rho = 917 # column-integrated glacier dentisty [kg/m3]
	g = 9.81 # accel due to gravity [m/s2]
	dx = np.diff(x)
	dz = np.diff(z)
	alpha = dz/dx
	alphaSmooth = gaussian_filter1d(alpha, sigma)
	H = yieldStress/(rho*g*np.abs(alpha))
	Hsmooth = yieldStress/(rho*g*np.abs(alphaSmooth))
	
	return H, Hsmooth, alpha, alphaSmooth

# Estimate ice surface velocity by shallow ice approximation
def calculateSIAspeed(H,alpha,f=1,A=2.4e-24,dt=86400):
	rho = 850
	g = 9.81
	n = 3
	taub = rho*g*H*f*np.abs(np.sin(np.arctan(alpha)))
	us = (2*A)/(n+1)*(taub**n)*H*dt
	return us

# Fn to read data from json file
def readJsonFile(jsonFn):
	with open(jsonFn) as data_file:
		data = json.load(data_file)
	return data

# Fn to interpolate evenly spaced vertices on a polyline
# As of 13 Jul 2016, these not evenly spaced, but I'm not super concerned about that right now.
def interpolateVerticesOfPolyline(lineGeom,dDist):
	clX,clY = getLineGeometryVerticesXY(lineGeom) # line vertices in utm
	dist = distanceAlongLine(clX,clY) # distance along centerline
	nP = len(clX) # number of vertices
	
	# Initialize with first point
	xInt = np.array(clX[0])
	yInt = np.array(clY[0])
	
	for i in range(1,nP): # iterate over vertices
		nPointSeg = np.round(np.array(dist[i]-dist[i-1])/dDist) # number of points that will be on line segment
		#print nPointSeg
		dE = clX[i] - clX[i-1] # change in easting across line seg
		dN = clY[i] - clY[i-1] # change in northing across line seg
		dx = np.array(dE)/nPointSeg # change in easting between interpolated points
		dy = np.array(dN)/nPointSeg # change in easting between interpolated points
		xSeg = clX[i-1] + np.arange(0,dE,dx) # x-coordinates for interpolated points on line
		ySeg = clY[i-1] + np.arange(0,dN,dy) # y-coordinates for interpolated points on line	
		
		xInt = np.insert(xInt,-1,xSeg)
		yInt = np.insert(yInt,-1,ySeg)
		
	return xInt, yInt
		

#%%	
## FILENAMES
#centerlineFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/qgis/Alaska_centerline_from_ChristianKleinholz/Centerlines.shp'
centerlineFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/qgis/keinholz_centerlines_paredDownToStudyArea.shp'
rgiFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/qgis/rgi_50_epsg3338.shp'
#demFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/dtm/ned_allTogether/img_files/wrangellsStElias_ned_merge.tif' # think this needs to be in diff coord sys?
demFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/dtm/wrangellStElias_gdem_EPSG102006.tif' # testing w old DEM that worked..
oldCenterlineFn = '/Users/wiar9509/Documents/generalScripts/swath/shapefileInputs/newWStETransects_19nov15_fromQGIS_EGPSG102006.shp'
jsonWithVelocitiesFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/nabeEB_evenlySpace_100m_profiles_sampled_2015-12-04.json'

# Contants
g = 9.81 # acceleration due to gravity [m/s2]
rhoi = 850 # column-integrated glacier density [kg/m3]
A = 2.4e-24 # flow law rate factor [pa^-3 s^-1]
n = 3 # flow late exponent
taub = 100e3 # assumed basal shear stress [Pa]
sigmaH = 2 # smoothing width [stdev] for gaussian [units of dx]

#%%
## ANALYSIS
# Open shapefiles
#clShp = ogr.Open(centerlineFn) # keinholz centerlins
clShp = ogr.Open(oldCenterlineFn) # my old hand-traced centerlines
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
	utmEpsg = getUtmEpsgFromGeometry(clGeom) # get epsg code for utm zone
	clGeomUtm = transformGeometry(clGeom,utmEpsg) # convert geometry to utm	
#	clX,clY = getLineGeometryVerticesXY(clGeomUtm) # line vertices in utm
#  	clXalbers,clYalbers=getLineGeometryVerticesXY(clGeom) # line vertices in alaska albers
 	
	if(1): # modifications to work with old centerline profiles where densifying points
		clX,clY = interpolateVerticesOfPolyline(clGeomUtm,100) # line vertices in utm, with line segments densified	
		interpClGeom = ogr.Geometry(ogr.wkbMultiPoint)
		srs = osr.SpatialReference()
		srs.ImportFromEPSG(utmEpsg)
		interpClGeom.AssignSpatialReference(srs)	
		for i in range(0,len(clX)):
			point = ogr.Geometry(ogr.wkbPoint)
			point.AddPoint(clX[i],clY[i])
			interpClGeom.AddGeometry(point)

		interpAlbers = transformGeometry(interpClGeom,3338) # interpolated centerline geometry in albers coordinates
		clXalbers,clYalbers=getLineGeometryVerticesXY(interpAlbers)
		
	dist = distanceAlongLine(clX,clY) # distance along centerline
	clWidth = calculateWidthAlongCenterline(clGeom,rgiLyr) # width along centerline
	
	# Get elevation along centerline
	xyTfm = createTransformation(demFn)
	elev = sampleRasterAtXY(demFn,clXalbers,clYalbers,xyTfm)

	splineFn = UnivariateSpline(dist,elev,s=2e5)
	elevSmooth = splineFn(dist)


	# Estimate ice thickness
	sigmaH = 5 # gaussian smoothing window to use for smooth surface slope in ice thickness from taub=constant estimate
	Hsharp50, Hsmooth50, alphaSharp, alphaSmooth = estimateIceThicknessYieldStress(dist,elevSmooth,.5e5,sigmaH)
	Hsharp100, Hsmooth100, alphaSharp, alphaSmooth = estimateIceThicknessYieldStress(dist,elevSmooth,1e5,sigmaH)
	Hsharp150, Hsmooth150, alphaSharp, alphaSmooth = estimateIceThicknessYieldStress(dist,elevSmooth,1.5e5,sigmaH)	
	
	# Clean bad data associated with 0 surface slope
	cleanInd = np.abs(np.diff(Hsmooth100))>100
	Hsmooth50[cleanInd] = np.nan
	Hsmooth100[cleanInd] = np.nan
	Hsmooth150[cleanInd] = np.nan
	
	# Estimate ice surface velocity by shallow ice approximation
	udef50 = calculateSIAspeed(Hsmooth50,alphaSmooth,f=1,A=2.4e-24,dt=86400)
	udef100 = calculateSIAspeed(Hsmooth100,alphaSmooth,f=1,A=2.4e-24,dt=86400)
	udef150 = calculateSIAspeed(Hsmooth150,alphaSmooth,f=1,A=2.4e-24,dt=86400)
	
	# plot it up
	plt.figure(0)
	ax1 = plt.subplot(311)
	ax1.plot(dist/1e3,elev,c='k')
	ax1.plot(dist/1e3,elevSmooth,c='b',lw=2)	
	ax2 = ax1.twinx()
	ax2.plot(dist[1:]/1e3,alphaSmooth,lw=2,c='r')
	for tl in ax1.get_yticklabels():
		tl.set_color('b')
	for t2 in ax2.get_yticklabels():
		t2.set_color('r')
	ax1.set_ylabel('Elevation [m]', color='b',fontsize=16)
	ax2.set_ylabel('Surface slope [%]', color='r',fontsize=16)		
	plt.subplot(312)
	plt.fill_between(dist[1:]/1e3,Hsmooth50,Hsmooth150,lw=1,facecolor='c',alpha=0.5)
	plt.plot(dist[1:]/1e3,Hsmooth100,lw=2,label='$\tau_b=1.0e5 Pa')	
	plt.ylim(0,1500)
	plt.ylabel('Ice thickness [m]',fontsize=16)
	
	
	plt.subplot(313)
	plt.fill_between(dist[1:]/1e3,udef50,udef150,lw=1,facecolor='c',alpha=0.5)
	plt.plot(dist[1:]/1e3,udef100,lw=2,label='$\tau_b=1.0e5 Pa')	
	plt.ylim(0,0.5)
	plt.ylabel('Ice velocity [m/d]',fontsize=16)
	plt.xlabel('Distance [km]',fontsize=16)	
	#plt.savefig('defm_velocity.pdf')
	plt.show()
	plt.close()
		
	if 0: # turned this off when found old fn
		# Calculate ice thickness
		dDist = np.diff(dist) # change in distance along profile [m]		
		elevSmooth = smooth(np.array(elev),41,window='flat') # smooth elevation data with a moving average
		dz = np.diff(elevSmooth) # change in elevation along profile [m]
		slope = dz/dDist # ice surface slope
		slope = np.insert(slope,0,0) # add a leading 0 to slope. No flux upstream at ice divide.
		dDist = np.insert(dDist,0,0) # add a leading 0 to dDist to treat spacing in smoothing
		Htaub = taub/(rhoi*g*np.sin(slope))
		udef = (2*A)/(n+1)*((rhoi*g*np.sin(slope))**n)*Htaub**(n+1)*86400 # expected deformation velocity [m/d]
	
			
fig, ax1 = plt.subplots()
ax1.plot(dist/1e3,np.array(elev), 'b-',linewidth=2)
ax1.set_xlabel('Down-glacier distance [km]',fontsize=16)
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Elevation [m]', color='b',fontsize=16)
ax1.set_ylim((800,4200))
ax1.set_yticks(np.arange(1000,4500,500))
for tl in ax1.get_yticklabels():
    tl.set_color('b')


ax2 = ax1.twinx()
ax2.plot(dist/1e3,clWidth/1e3, 'r-',linewidth=2)
ax2.set_ylabel('Width [km]', color='r',fontsize=16)
ax2.set_ylim((-5,30))
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.savefig('nabesna_widthAndSlope_07jul2016.pdf')
plt.show()		
		
			
		
if 0: # testing affect of smoothing parameter on elevation profile	
	s = np.arange(1e5,3.5e5,0.5e5)

	plt.figure(0)
	for i in range(0,len(s)):
		splineFn = UnivariateSpline(dist,elev,s=s[i])
		elevSmooth = splineFn(dist)
		plt.plot(dist,elevSmooth,lw=2,label=s[i])
	
	plt.plot(dist,elev,c='k',lw=1,label='0')
	plt.legend()
	plt.xlabel('Distance [m]',fontsize=16)
	plt.ylabel('Elevation [m]',fontsize=16)
	plt.title('Effect of spline smoothing length, s, on elevation profile')
	plt.savefig('effectOfSplineSmoothingParamOnElevationProfile.pdf')
	plt.show()
	plt.close()
	
	
	
		


