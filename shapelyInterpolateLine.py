'''
shapely line interpolation
'''

import os
from shapely import geometry
import fiona
import numpy as np

shpFn = '/Users/wiar9509/Documents/generalScripts/swath/shapefileInputs/newWStETransectsUtm7n_24sep15.shp'
dDist = 100

shp = fiona.open(shpFn)

i = 0

feat = shp[i]
transName = feat['properties']['name']
geom = feat['geometry']
line = geometry.shape(geom)
x,y=line.xy

def interpolateLineGeom(line,dDist):
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
	
xInt,yInt = interpolateLineGeom(line,dDist)

plt.plot(xInt,yInt,'bv',lw=0)
plt.plot(x.tolist(),y.tolist(),'ro',lw=0)
	
	
	
	