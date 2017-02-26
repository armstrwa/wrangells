'''
Script to analyze best fully automated results, now that I've done a lot of cherry picking..
20 nov 2016
william armstrong
'''

import json
import matplotlib.pyplot as plt
import numpy as np
import wrangellsGeospatialTools as wgt 
import matplotlib.gridspec as gridspec
import glob


## USER DEFINED
#jsonFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/kennCL_swathVelocitySampling_2016-10-03.json'
#jsonFolder =  '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/'
jsonFolder = '/Volumes/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/'
#jsonFolder = '/Users/wiar9509/Google Drive/wrangells/json/'
countMin = 10 # minimum pixels within swath profiler box to be considered good data
stdMax = 1 # maximum standard devuation of velocity within swath profiler box to be considered good data
summerStartEnd=((150,250)) # start and end DOYs to define the summer season for seasonal velocity split
elevationProfileFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/new_centerlines_withChugach_epsg102006_elevationProfiles.json'

## PROCESSING
## Initializing
velocityChangeDict = {}
velocityChangeDict2 = {}
allSeasonalVelocities = {}
jsonList = glob.glob(jsonFolder + '*2016-10-05.json')

for jsonFn in jsonList:
	print jsonFn


	profile_data = wgt.readJsonFile(jsonFn) # read in json velocity data
	transName = profile_data['transect_name'] # get transect name	
	print "Processing " + transName
	nProfiles = len(profile_data['profiles']) # number of profiles (time slices) within json data
	dist = profile_data['sample_pts_frontdist'] # distance from headwall [m]
	termDist = wgt.flipDistance(dist) # Flipping distance axis, so zero at terminus
	normTermDist = wgt.normalizeDistance(termDist) # Normalizing distance by glacier length
	profile_data  = wgt.swathFilter(profile_data,countMin,stdMax) # filter velocity results based on counts and standard deviation within swath profiler boxes
	
	
	# Temporal filtering
	uMatxRaw =wgt.makeVelocityMatrix(profile_data,'speed') # make velocity matrix (nCols = len(space), nRows = len(time) )
	uMatxSwath = wgt.makeVelocityMatrix(profile_data,'speedMeanClean') # make velocity matrix (nCols = len(space), nRows = len(time) )
	uMatxRawOut = wgt.removeOutliers(uMatxRaw) # remove outliers
	uMatxSwathOut = wgt.removeOutliers(uMatxSwath) # remove outliers
	uMatxSwathOut2 = wgt.removeOutliers(uMatxSwathOut) # second iteration removing outliers
	# All 2 suffixes below indicate 2 iterations of outerlier removal

	# Put temporally filtered data back into profile_data dict form
	profileDataSwathOutliers = wgt.putTempFiltDataBack(profile_data,uMatxSwathOut)
	profileDataSwathOutliers2 = wgt.putTempFiltDataBack(profile_data,uMatxSwathOut2)
	
	# Break velocities into summer and winter categories
	seasonalVelocityDict = wgt.splitSummerWinter(profileDataSwathOutliers,summerStartEnd)
	seasonalVelocityDict2 = wgt.splitSummerWinter(profileDataSwathOutliers2,summerStartEnd)
	speedupDictNow = wgt.calculateSeasonalSpeedup(seasonalVelocityDict,termDist)
	speedupDictNow2 = wgt.calculateSeasonalSpeedup(seasonalVelocityDict2,termDist)
	velocityChangeDict[transName] = speedupDictNow
	velocityChangeDict2[transName] = speedupDictNow2
	allSeasonalVelocities[transName]=seasonalVelocityDict2
	# Plot seasonal velocity data
	#wgt.plotSeasonalVelocity(seasonalVelocityDict2,speedupDictNow2,termDist,1,1)
	#wgt.plotSeasonalSpeedup(speedupDictNow2,termDist,1,1)	

# Write  speedup plots out
jsonFnOut1 = 'speedupDict_fullyAutomated_28nov2016.json'
jsonFnOut2 = 'speedupDict_2iterTempFilter_fullyAutomated_28nov2016.json'
wgt.writeJsonFile(velocityChangeDict,jsonFnOut1)
wgt.writeJsonFile(velocityChangeDict2,jsonFnOut2)


# Write seasonal velocity out
jsonFnOut3 = 'seasonalVelocityProfiles_2iterTempFilter_fullyAutomated_24jan2017.json'
wgt.writeJsonFile(allSeasonalVelocities,jsonFnOut3)

# Define speedup distance
from scipy.stats import ranksums
speedupDict = wgt.readJsonFile(jsonFnOut1)
speedupDict2 = wgt.readJsonFile(jsonFnOut2)

# Identify 'good' (i.e., well-behaved) transects
goodTransectsFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/excel/summerSpeedupTransectList_yesSpeedup.csv'
#goodTransectsFn = 'summerSpeedupTransectList_yesSpeedup.csv'
transList = np.genfromtxt(goodTransectsFn,delimiter=',',skip_header=0,dtype=None)
transNames = transList[:,0]
ranksumDict = {} # for storage

	
# Iterate over transects in speedup dict
keyIter = speedupDict2.iterkeys() # key iterator
for i in range(0,len(speedupDict2)):
	transNow = keyIter.next() # current transect nane
	if transNow in transNames: # only do good transects
		print transNow
		
		speedupNow = speedupDict2[transNow]['tempFilter']
		termDist = speedupDict2[transNow]['termDist']
		minDiff = speedupDict2[transNow]['tempFilter']['minDiff']
		medDiff = speedupDict2[transNow]['tempFilter']['medianDiff']
		maxDiff = speedupDict2[transNow]['tempFilter']['maxDiff']

		# Iterative rank sum
		wList = []
		pList = []
		for i in range(0,len(termDist)):
			# Split into up and down glacier regions
			cutoffDist = termDist[i] # distance for separating up and down glacier
			upGlacierInd = np.arange(0,i+1)
			downGlacierInd = np.arange(i+1,len(termDist))
			medDiffUp = np.array(medDiff)[upGlacierInd]
			medDiffDown = np.array(medDiff)[downGlacierInd]	
			# ranksum test for whether upglacier is statistically different from up-glacier
			w,p = ranksums(medDiffUp,medDiffDown) # w = test statistic, p = significance level
			# Append stats to list
			wList.append(w)
			pList.append(p)

		breakPoint = np.where(np.array(wList)==np.nanmin(wList))
		speedupLength = np.array(termDist)[breakPoint]	
		
		ranksumDict[transNow] = {'termDist':termDist,'medSpeedup':medDiff,'rankSumStat':wList,'rankSumPval':pList,'speedupLength':speedupLength.tolist()}
		
		if 1:
			plt.subplot(3,1,1)
			plt.plot(np.array(termDist)/1e3,medDiff,lw=1,c='gray')
			plt.plot(np.array((speedupLength,speedupLength))/1e3,(-0.5,0.5),ls='--',c='k')
			plt.plot((0,np.max(termDist)/1e3),(0,0),ls='--',c='k')
			upGlacierInd = np.arange(0,breakPoint[0]+1)
			downGlacierInd = np.arange(breakPoint[0]+1,len(termDist))
			plt.scatter(np.array(termDist)[upGlacierInd]/1e3,np.array(medDiff)[upGlacierInd],c='b',s=20)
			plt.scatter(np.array(termDist)[downGlacierInd]/1e3,np.array(medDiff)[downGlacierInd],c='r',s=20)
			plt.ylabel('speedup')
			plt.xlim((0,np.max(termDist)/1e3))
			plt.ylim((-.4,.4))
			plt.title(transNow + ' rank sum for speedup length')
			plt.gca().invert_xaxis()
			plt.subplot(3,1,2)
			plt.plot(np.array(termDist)/1e3,wList,lw=2)
			plt.ylabel('wilcoxon test stat.')
			plt.xlim((0,np.max(termDist)/1e3))
			plt.gca().invert_xaxis()			
			plt.subplot(3,1,3)
			plt.plot(np.array(termDist)/1e3,pList,lw=2)
			plt.ylabel('p val')
			plt.xlabel('distance from terminus [km]')
			plt.xlim((0,np.max(termDist)/1e3))
			plt.gca().invert_xaxis()			
			plt.savefig(transNow + 'ranksumForSpeedupLength.pdf')
			plt.show()
			plt.close()

ranksumFnOut='/Users/wiar9509/Google Drive/wrangells/json/ranksumSpeedupDict_30nov2016.json'
ranksumFnOut='ranksumSpeedupDict_30nov2016.json'
wgt.writeJsonFile(ranksumDict,ranksumFnOut)	

# Load in rank sum dict
ranksumDict = wgt.readJsonFile(ranksumFnOut)
keyIter = ranksumDict.iterkeys()
# Coloring for lines
transNums = np.linspace(0, 1, len(ranksumDict))
transCols = plt.cm.Set1(transNums)
colIter = 0
if 0:
	for i in range(0,len(ranksumDict)):
		transNow=keyIter.next()
		print transNow
		ranksumNow = ranksumDict[transNow]
		speedupLength = np.array(ranksumNow['speedupLength'])/np.max(ranksumNow['termDist'])
		glacierLength = np.max(ranksumNow['termDist'])/1e3
		
		# manually correcting kluvesna, that has bad RS results b/c ice fall
		if transNow == 'kluvCL':
			speedupLength = 12000./np.max(ranksumNow['termDist'])
			print transNow, glacierLength, speedupLength
					
		plt.scatter(glacierLength,speedupLength,lw=1,color='gray',facecolor=transCols[colIter],s=60,label=transNow)
		colIter+=1

	#plt.axis('equal')
	plt.xlim((0,80))
	plt.xlabel('Glacier length [km]',fontsize=18)
	plt.ylabel('Speedup length ratio [-]',fontsize=18)
	plt.ylim((0,1))
	plt.legend(prop={'size':8},frameon=False,ncol=3,loc=3,numpoints=1,scatterpoints=1)
	plt.savefig('glacierLength_vs_speedupLength_normalize_fixKluvCL_bigDots.pdf')
	plt.show()
	plt.close()	

if 1:
	fig = plt.figure()
	ax = fig.add_subplot(111)
	xList = []
	yList = []
	for i in range(0,len(ranksumDict)):
		transNow=keyIter.next()
	
		ranksumNow = ranksumDict[transNow]
		speedupLength = np.array(ranksumNow['speedupLength'])/1e3
		glacierLength = np.max(ranksumNow['termDist'])/1e3
		
		print transNow, speedupLength, glacierLength
		xList.append(glacierLength)
		yList.append(speedupLength)
		plt.scatter(glacierLength,speedupLength,lw=0.5,color='k',facecolor=transCols[colIter],s=60,label=transNow)
		xy = ((glacierLength+np.random.randint(0,3),speedupLength-np.random.randint(-2,2)))
		#ax.annotate(transNow, xy=xy, textcoords='data',fontsize=8) 
		colIter+=1
	
# 	iterCount = 0
# 	for xy in zip(xList,yList):
# 		xy=np.array(xy)
# 		xy[0]=xy[0]+np.array(np.random.randint(0,3)).astype('float')
# 		xy[1]=xy[1]+np.array(np.random.randint(-3,0)).astype('float')	
# 		nameNow = transNames[iterCount]                                      # <--
# 		#ax.annotate(nameNow, xy=xy, textcoords='data',fontsize=8) 
# 		iterCount+=1
	plt.plot((0,80),(0,80),ls='--',c='k')
	#plt.axis('equal')
	plt.xlim((0,80))
	plt.xlabel('Glacier length [km]',fontsize=18)
	plt.ylabel('Speedup length [km]',fontsize=18)
	plt.ylim((0,80))
	plt.legend(prop={'size':8},frameon=False,ncol=3,loc=2,numpoints=1,scatterpoints=1)
	plt.savefig('glacierLength_vs_speedupLength_bigDots.pdf')
	plt.show()
	plt.close()	

## SQUIGGLE PLOT (seasonal speedup collective profiles)
# Make squiggle plot (speedup)
speedupDict = wgt.readJsonFile(jsonFnOut1)
speedupDict2 = wgt.readJsonFile(jsonFnOut2)

# Identify 'good' (i.e., well-behaved) transects
goodTransectsFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/excel/summerSpeedupTransectList_yesSpeedup.csv'
transList = np.genfromtxt(goodTransectsFn,delimiter=',',skip_header=0,dtype=None)
transNames = transList[:,0]

# Count cutoff for plotting velocity change
# I.e., minimum number of observations at a point
minCountToPlot = 4

# Coloring for lines
#transNums = np.linspace(0.25, 1, len(transNames))
#transCols = plt.cm.Accent(transNums)
transNums = np.linspace(0, 1, len(transNames))
transCols = plt.cm.Set1(transNums)

# Iterate over transects in speedup dict
keyIter = speedupDict2.iterkeys() # key iterator
colIter = 0 # for line color
for i in range(0,len(speedupDict2)):
	transNow = keyIter.next() # current transect nane
	if transNow in transNames: # only do good transects
		# Pulling data
		speedupNow = speedupDict2[transNow]['tempFilter']
		termDist = speedupDict2[transNow]['termDist']
		minDiff = speedupDict2[transNow]['tempFilter']['minDiff']
		medDiff = speedupDict2[transNow]['tempFilter']['medianDiff']
		maxDiff = speedupDict2[transNow]['tempFilter']['maxDiff']

		# Hiding data based on observation count
		winterDataNum = speedupNow['winterDataNum']
		summerDataNum = speedupNow['summerDataNum']
		
		lowCountInd = np.logical_or(np.array(winterDataNum)<minCountToPlot,np.array(summerDataNum)<minCountToPlot) # where is there not enough data in one season

		medDiffMa = np.ma.masked_array(medDiff,lowCountInd)
		minDiffMa = np.ma.masked_array(minDiff,lowCountInd)
		maxDiffMa = np.ma.masked_array(maxDiff,lowCountInd)

		if transNow in {'rootCL'}:
			icefallInd = np.logical_and(np.array(termDist)<23.*1e3,np.array(termDist)>17.8*1e3) # where is there bad data associated with ice falls
			bothInds = np.logical_or(icefallInd == True,lowCountInd == True) # where is there either too few data or icefall
			#medDiffMa = np.ma.masked_array(medDiff,bothInds)
			medDiffMa = np.copy(medDiff)
			medDiffMa[bothInds]=np.nan
			maxDiffMa = np.copy(maxDiff)
			maxDiffMa[bothInds] = np.nan
			minDiffMa = np.copy(minDiff)
			minDiffMa[bothInds] = np.nan			
# 			minDiffMa = np.ma.masked_array(minDiff,bothInds)
# 			maxDiffMa = np.ma.masked_array(maxDiff,bothInds)			
		medDiffSmMa = wgt.smooth(medDiffMa,9)
		
								
		plt.plot(np.array(termDist)/np.max(termDist),medDiffSmMa,c=transCols[colIter],lw=2,label=transNow,alpha=0.7)
		plt.fill_between(np.array(termDist)/np.max(termDist),minDiff,maxDiff,color=transCols[colIter],alpha=0.1,lw=0)
		colIter+=1
		
plt.plot((0,80),(0,0),lw=2,c='k',ls='--')
plt.ylabel('Speedup [m d$^{-1}$]',fontsize=18)
plt.xlabel('Distance from terminus [km]',fontsize=18)
plt.legend(loc=3,frameon=False,ncol=3,prop={'size':8})
plt.ylim((-.4,.4))
plt.xlim((0,1))
plt.gca().invert_xaxis()
plt.savefig('speedupSquiggles_withFill_fullAuto_exludeIfLt4Obs_onlyWithYesSpeedup_normDist_icefallNans.pdf')
plt.show()
plt.close()


## ELEVATION PROFILES

elevData = wgt.readJsonFile(elevationProfileFn)
keyIter = elevData.iterkeys()
subsetTransFn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/excel/summerSpeedupTransectList_yesAndNoSpeedup.csv'
transList = np.genfromtxt(subsetTransFn,delimiter=',',skip_header=0,dtype=None)
transNames = transList[:,0]
speedupType = transList[:,1] # 1 for yes speedup
alphaVal = 0.8
speedupElevList = [] # for storage
speedupDistList = []
allElevList = 1e3
allDistList = 0

for i in range(0,len(elevData)):
	transNow = keyIter.next()
	# Just plot transects with speedups
	if transNow in transList[:,0]:
		transInd = np.where(transList[:,0]==transNow)
		speedType = transList[transInd,1]
		elevNow = elevData[transNow]
		zMean = elevNow['elevMean']
		slope = elevNow['slope']
		distUp = wgt.flipDistance(elevNow['dist']) # not exact same as ranksumNow['termDist']
		distDown = elevNow['dist'] # distance from headwall
		
		
		try:
			ranksumNow = ranksumDict[transNow]
			bp = ranksumNow['speedupLength']
		
			downElevInd = np.where(np.array(distUp)<=bp)
			upElevInd = np.where(np.array(distUp)>=bp)
			downElevInd2 = np.where(np.array(distDown)<=bp)
			upElevInd2 = np.where(np.array(distDown)>=bp)
			
			upGlacierDist = distUp[upElevInd]
			downGlacierDist = distUp[downElevInd]			
			upGlacierElev = np.array(zMean)[upElevInd]
			downGlacierElev = np.array(zMean)[downElevInd]
			
			# Add a point in between to connect the two lines
			distBetween = np.mean( (np.min(upGlacierDist),np.max(downGlacierDist)) )			
			elevBetween = np.mean( (np.min(upGlacierElev),np.max(downGlacierElev)) )
			
			downGlacierDist2 = np.insert(downGlacierDist,0,distBetween)
			downGlacierElev2 = np.insert(downGlacierElev,0,elevBetween)
			upGlacierDist2 = np.append(upGlacierDist,distBetween)
			upGlacierElev2 = np.append(upGlacierElev,elevBetween)
			
			# store speedup maximum elevation
			speedupElevList.append(elevBetween)
			speedupDistList.append(distBetween)			
			
			allElevList = np.append(allElevList,zMean)
			allDistList = np.append(allDistList,distUp)
			
			plt.plot(np.array(downGlacierDist2)/1e3,downGlacierElev2,alpha=alphaVal,c='r',lw=2)
			plt.plot(np.array(upGlacierDist2)/1e3,upGlacierElev2,alpha=alphaVal,c='b',lw=2)
		except KeyError: # this happens for non-speedup types
			plt.plot(np.array(distUp)/1e3,np.array(zMean),alpha=alphaVal,c='gray',lw=2)

		#plt.annotate(transNow,(np.max(np.array(distUp)/1e3),np.max(np.array(zMean))+np.random.randint(-50,50)),fontsize=8)


s25,s50,s75 = np.percentile(speedupElevList,(25,50,75))
all25,all50,all75 = np.percentile(allElevList,(25,50,75))

#plt.fill_between((-1,80),(all25,all25),(all75,all75),facecolor='gray',alpha=0.2,linewidth=0)
plt.fill_between((-1,80),(s25,s25),(s75,s75),facecolor='darkred',alpha=0.2,linewidth=0)
#plt.plot((-1,80),(all50,all50),c='gray',ls='--',lw=1)		
plt.plot((-1,80),(s50,s50),c='darkred',ls='--',lw=1)		

plt.plot(-10,-10,c='r',alpha=alphaVal,lw=2,label='Summer speedup reach')
plt.plot(-10,-10,c='b',alpha=alphaVal,lw=2,label='No speedup reach')
plt.plot(-10,-10,c='gray',alpha=alphaVal,lw=2,label='No speedup glacier')
plt.xlim((-1,80))
plt.ylim((500,4700))
plt.xlabel('Distance from terminus [km]',fontsize=18)
plt.ylabel('Elevation [m]',fontsize=18)
plt.title('Elevation profiles and speedup',fontsize=20)
plt.legend(frameon=False,loc=4)
plt.savefig('elevAndSpeedup_bottomUp_noNames_wIQR.pdf')
plt.show()
plt.close()



### SPEEDUP BOX AND WHISKER

# Elevation
if 1:
	locations = [1,2]
	names = ['Speedup','All']

	plt.figure()
	plot2 = plt.boxplot([speedupElevList,allElevList], 
						widths=0.7,
						notch=False,             # adds median notch
						positions=locations,    # boxes locations
						patch_artist=True,
						)
					
	plt.grid(axis='y',          # set y-axis grid lines
			linestyle='--',     # use dashed lines
			which='major',      # only major ticks
			color='lightgrey',  # line colour
			alpha=0.7)          # make lines semi-translucent

	#for box, colour in zip(plot2['boxes'], box_colours):
	plt.setp(plot2['boxes'][0], 
			 linewidth=1.5, 
			 facecolor='darkred')
		 
	plt.setp(plot2['boxes'][1], 
		 linewidth=1.5, 
		 facecolor='lightgray')
		 
	plt.setp(plot2['whiskers'], color='gray', linewidth=1.5)
	plt.setp(plot2['caps'], color='gray', linewidth=1.5)
	plt.setp(plot2['fliers'], color='gray', marker='o', markersize=10)
	plt.setp(plot2['medians'], color='salmon', linewidth=1.5)


	plt.xticks(locations,               # tick marks
			   names,                   # labels
			   rotation='horizontal',fontsize=16)     # rotate the labels
	plt.yticks(fontsize=16)

	plt.ylabel('Elevation [m]',fontsize=18)                  # y-axis label
	plt.ylim((500,4700))
	plt.savefig('speedupElev_vs_allElev_boxPlot.pdf')
	plt.show()
	plt.close()

# Distance
if 1:										 
	locations = [1,2]
	names = ['Speedup','All']

	plt.figure()
	plot2 = plt.boxplot([np.array(speedupDistList)/1e3,np.array(allDistList)/1e3], 
						widths=0.7,
						notch=False,             # adds median notch
						positions=locations,    # boxes locations
						patch_artist=True,
						)
					
	plt.grid(axis='y',          # set y-axis grid lines
			linestyle='--',     # use dashed lines
			which='major',      # only major ticks
			color='lightgrey',  # line colour
			alpha=0.7)          # make lines semi-translucent

	#for box, colour in zip(plot2['boxes'], box_colours):
	plt.setp(plot2['boxes'][0], 
			 linewidth=1.5, 
			 facecolor='darkred')
		 
	plt.setp(plot2['boxes'][1], 
		 linewidth=1.5, 
		 facecolor='lightgray')
		 
	plt.setp(plot2['whiskers'], color='gray', linewidth=1.5)
	plt.setp(plot2['caps'], color='gray', linewidth=1.5)
	plt.setp(plot2['fliers'], color='gray', marker='o', markersize=10)
	plt.setp(plot2['medians'], color='salmon', linewidth=1.5)


	plt.xticks(locations,               # tick marks
			   names,                   # labels
			   rotation='horizontal',fontsize=16)     # rotate the labels
	plt.yticks(fontsize=16)

	plt.ylabel('Distance from terminus [km]',fontsize=18)                  # y-axis label

	plt.savefig('speedupDist_vs_allDist_boxPlot.pdf')
	plt.show()
	plt.close()
                                         


## PETTIT TEST
# 12 jan 2017

def sign(num):
	'''
	get the sign of a number
	input = num (a number)
	output = sgnOut = -1/0/1
	'''
	if num > 0:
		sgnOut = 1
	elif num == 0:
		sgnOut = 0 
	elif num < 0:
		sgnOut = -1
	return sgnOut
	
def pettitTest(data):
	'''
	Pettit change point test
	Input: data = numeric data list or array
	Output: dataOut = {'testStat','changePointValue','changePointInd'}
	Output: testStat = array of test statistics (len = data), changePointValue = max(abs(testStat)), changePointInd = index value for where change point occurs
	'''
	nPoints = len(data) # number of data points
	uAllList = [] # initialize list of test statistics
	# iterate over data points
	for i in range(0,nPoints):
		uNowList = [] # initialize
		# iterate over following data points
		for j in range(i+1,nPoints):
			uNow = sign(data[i]-data[j]) # take sign of ith data point minus jth data point
			uNowList.append(uNow) # add sign to a list for this i
	
		uAllList.append(np.sum(uNowList)) # sum the signs for this i

	# Calculate p value
	maxU = np.max(uAllList)
	# pVal doesn't seem to be right.. should there not be a 2 in front? Is <1 then.
	pVal = 2 * np.exp( (-6*(maxU**2)) / ( nPoints**3 + nPoints**2 ) )
	
	# Find change point
	changePointInd = np.where(uAllList==maxU)[0]
	changePointValue = data[changePointInd]
	dataOut = {'testStats':uAllList,'changePointValue':changePointValue,'changePointInd':int(changePointInd),'pValue':pVal}

	return dataOut


rsSpeedupLength = ranksumDict['kennCL']['speedupLength']
ptSpeedupLength = termDist[dataOut['changePointInd']]
plt.subplot(2,1,1)
plt.plot(termDist,data)	
plt.plot((ptSpeedupLength,ptSpeedupLength),(-0.5,0.5),ls='--',c='r',label='Pettit')
plt.plot((rsSpeedupLength,rsSpeedupLength),(-0.5,0.5),ls='--',c='b',label='Ranksum')
plt.legend(loc=3,frameon=False)
plt.ylabel('Speedup [m/d]',fontsize=18)
ax1 = plt.subplot(2,1,2)
plt.plot(termDist,dataOut['testStats'],c='r',lw=2)
#ax1.yaxis.label.set_label('Pettit test stat')
ax1.set_ylabel('Pettit test stat', color='r',fontsize=18)
ax1.yaxis.label.set_color('r')
ax2 = ax1.twinx()
ax2.plot(termDist,ranksumDict['kennCL']['rankSumStat'],c='b',lw=2)
#ax2.yaxis.label.set_label('Ranksum test stat')
ax2.set_ylabel('Ranksum test stat', color='b',fontsize=18)
plt.xlabel('Distance from terminus [m]',fontsize=18)
plt.savefig('ranksumVsPettitChangePoints_kennCL.pdf')
plt.show()
plt.close()


