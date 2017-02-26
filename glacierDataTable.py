# make table with glacier data
## SCRIPT TO READ SHAPEFILE ATTRIBUTE TABLE AND WRITE TO JSON
#  24 February 2016
#  William Armstrong

from osgeo import ogr
import numpy as np
import matplotlib.pyplot as plt
import json

shpFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/qgis/studyGlaciers_onlyGlaciersWithTransects_13jan2017.shp'

writeToggle = 0 # 1 to write json file


def readRipRGI(shpFn,writeToggle=1):
	'''
	Function to read in the randolph glacier inventory and turn it into a json file
	Inputs: shpFn = full filepath to RGI shapefile. writeToggle = 1 to write json file, 0 to only create variable
	
	29 February 2016
	William Armstrong
	'''

	shp = ogr.Open(shpFn)

	lyr = shp.GetLayer(0)

	featCount = lyr.GetFeatureCount()

	# initializing
	gIterIds = []
	gNames = []
	gAreas = []
	gCenLons = []
	gCenLats = []
	gZmins = []
	gZmeds = []
	gZmaxs = []
	gSlopes = []
	gAspects = []
	gLens = []
	gGlacTypes = []
	rgiIds = []
	
	# Initialize output JSON file
	glacierData = []

	for i in range(0,featCount):
		feat = lyr.GetFeature(i)
		gIterIds.append(i)
		if feat.GetField("Name") == '\x0e': # no data value for some reason
			gNames.append(np.nan)
			nameNow = np.nan
		else:
			nameNow = feat.GetField("Name")
			gNames.append(nameNow)
		rgiId = feat.GetField("RGIId")
		rgiIds.append(rgiId)
		areaNow = feat.GetField("Area")
		gAreas.append(areaNow)
		lonNow = feat.GetField("CenLon")
		gCenLons.append(lonNow)
		latNow = feat.GetField("CenLat")
		gCenLats.append(latNow)
		minNow = feat.GetField("Zmin")
		gZmins.append(minNow)
		medNow = feat.GetField("Zmed")
		gZmeds.append(medNow)
		maxNow = feat.GetField("Zmax")
		gZmaxs.append(maxNow)
		slopeNow = feat.GetField("Slope")
		gSlopes.append(slopeNow)
		aspNow = feat.GetField("Aspect")
		gAspects.append(aspNow)
		lenNow = feat.GetField("Lmax")
		gLens.append(lenNow)
	
		# Write dict for saving to JSON
		gDict = { 'iterId': i, 'glacierName':nameNow, 'area':areaNow, 'centerLon':lonNow, 'centerLat':latNow, 'elev_min':minNow, \
				'elev_med':medNow, 'elev_max':maxNow, 'slope':slopeNow, 'aspect':aspNow, 'length':lenNow,'rgiId':rgiId}
		# Append to list of other glacier data
		glacierData.append(gDict)
	
	shp.Destroy()
	glacierDict = {'glaciers':glacierData,'lengths':gLens,'areas':gAreas,'zMaxs':gZmaxs,'zMins':gZmins,'zMeds':gZmeds,'lats':gCenLats,'lons':gCenLons,'rgiIds':rgiIds}

	if writeToggle == 1:
		print "Saving data to: " + outFilename
		with open(outFilename, 'w') as f: 
			f.write(json.dumps(glacierDict))

	return glacierDict

studyGlaciersDict = readRipRGI(shpFn,0)

speedupTypeFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/excel/glacierSpeedupTypeList_withRGIids.csv'
speedupTypes = np.genfromtxt(speedupTypeFn,delimiter=',',skip_header=1,dtype=None)

test = np.vstack((studyGlaciersDict['rgiIds'],studyGlaciersDict['lats'],studyGlaciersDict['lons'],studyGlaciersDict['lengths'],studyGlaciersDict['areas'],studyGlaciersDict['zMins'],studyGlaciersDict['zMeds'],studyGlaciersDict['zMaxs']))

	for i in range(0,len(studyGlaciersDict['rgiIds'])):
		idNow = studyGlaciersDict['rgiIds'][i]
		for j in range(0,len(speedupTypes)):
			if idNow == speedupTypes[j][0]:
				speedupType = speedupTypes[j][1]
				if speedupType == 'yes':
					plt.scatter(np.array(studyGlaciersDict['lengths'])[i]/1e3,studyGlaciersDict['areas'][i],marker='.',c='darkred',alpha=0.9,lw=1,s=400)
				elif speedupType == 'no':
					plt.scatter(np.array(studyGlaciersDict['lengths'])[i]/1e3,studyGlaciersDict['areas'][i],marker='.',c='deepskyblue',alpha=0.9,lw=1,s=400)
				elif speedupType == 'mixed':
					plt.scatter(np.array(studyGlaciersDict['lengths'])[i]/1e3,studyGlaciersDict['areas'][i],marker='.',c='yellow',alpha=0.9,lw=1,s=400)
				elif speedupType == 'exclude':
					plt.scatter(np.array(studyGlaciersDict['lengths'])[i]/1e3,studyGlaciersDict['areas'][i],marker='.',c='dimgray',alpha=0.9,lw=1,s=400)