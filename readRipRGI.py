## SCRIPT TO READ SHAPEFILE ATTRIBUTE TABLE AND WRITE TO JSON
#  24 February 2016
#  William Armstrong

from osgeo import ogr
import numpy as np
import matplotlib.pyplot as plt
import json

shpFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/qgis/01_rgi50_Alaska/01_rgi50_Alaska.shp'
outFilename = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/qgis/rgi50_alaska.json'

writeToggle = 1 # 1 to write json file


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
				'elev_med':medNow, 'elev_max':maxNow, 'slope':slopeNow, 'aspect':aspNow, 'length':lenNow}
		# Append to list of other glacier data
		glacierData.append(gDict)
	
	shp.Destroy()
	glacierDict = {'glaciers':glacierData}

	if writeToggle == 1:
		print "Saving data to: " + outFilename
		with open(outFilename, 'w') as f: 
			f.write(json.dumps(glacierDict))



## PLOTTING
ind = np.isfinite(np.log10(gLens)) & np.isfinite(np.log10(gLens))
p = np.polyfit(np.log10(np.array(gLens)[ind]),np.log10(np.array(gAreas)[ind]),1)
yFit = np.polyval(p,np.log10(np.array(gAreas)[ind]))

fig = plt.figure(0)
ax = fig.add_subplot(1,1,1)
plt.plot(np.array(gLens)/1e3,gAreas,marker='.',lw=0)
#plt.plot(np.array(gLens)[ind]/1e3,yFit,lw=2,color='r')
ax.set_yscale('linear')
ax.set_xscale('linear')
plt.xlabel('Glacier length [km]',fontsize=16)
plt.ylabel('Glacier area [km$^2$]',fontsize=16)
plt.xlim([0, np.max(np.array(gLens)/1e3)])
plt.grid()
plt.title('Glacier length-area relation from RGI 5.0',fontsize=20)
#plt.savefig('alaskaGlacierAreaLength.pdf',dpi=300)
#plt.show()
#plt.close()

dA = .5
#nBins = np.ceil(np.max(gAreas)/dA)
nBins = 200
#bins = np.linspace(0,np.max(gAreas),nBins)
bins = np.logspace(-2,np.log10(np.max(gAreas)),nBins)
nArea = np.zeros((len(bins),1))
pArea = np.zeros((len(bins),1))

for i in range(0,len(bins)):
	binNow = bins[i]
	boolArr = gAreas>=binNow
	nArea[i] = np.sum(boolArr)

pArea = nArea/len(gAreas)

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1)
plt.plot(bins,pArea,marker='.',lw=0)
ax1.set_xscale('log')
plt.xlabel('Glacier area [km$^2$]',fontsize=16)
plt.ylabel('P(Area)',fontsize=16)
plt.title('Alaska Glaciers Area Probability, RGI 5.0',fontsize=20)
plt.grid()
#plt.savefig('alaskaGlaciersAreaProbability.pdf',dpi=300)
#plt.show()
#plt.close()

fig.show()
fig1.show()





	
	