### SCRIPT TO COMPUT STATS FROM JSON FILE
#   19 Nov 2015
#   William Armstrong

import json
import matplotlib.pyplot as plt
import numpy as np
from osgeo import ogr
import matplotlib.gridspec as gridspec
import scipy.stats as stats

annualToggle=0
myDataToggle=1
#startSeasons=[306, 92, 153, 242] # start days of each season (306 = Nov 1; 92 = Apr 1; 153 = Jun 1; 242 = Sep 1)
startSeasons=[304, 81, 121, 265] # start days of each season, these based of local met data (304 = Oct 28; 81 = Mar 22; 121 = May 1; 265 = Sep 22)
shapefileIn='/Users/wiar9509/Documents/generalScripts/swath/shapefileInputs/newWStETransects_19nov15_fromQGIS_EGPSG102006.shp'
elevationProfileFolder = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/dtm/elevationProfiles/'

def getNumDataFromJson(jsonData):
	shapefileName = jsonData['profile_shapefile']
	transName = shapefileName.split('_')[0]
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

# Function that takes a matrix of velocities size = (len(x), numProfiles).
def findWhenMedianVelocityStabilizesAllData(uMatx,iterMax=10):

	uMatx = uMatx.transpose() # starting writing this wrong, thought new profiles were on rows, but actually columns. Flipping to have profiles on rows
	
	numProfiles = uMatx.shape[0]
	lenX = uMatx.shape[1]
	xInd = np.arange(0,lenX)
	
	# Initialize
	sampleNow = np.empty((numProfiles,lenX,iterMax))
	medianNow = np.empty((numProfiles,lenX,iterMax))
	diffFromLast = np.empty((numProfiles,lenX,iterMax))
	diffFromFinal = np.empty((numProfiles,lenX,iterMax))	
	medianNow.fill(np.nan)
	sampleNow.fill(np.nan)
	diffFromLast.fill(np.nan)	
	diffFromFinal.fill(np.nan)
	
	# Median velocity at x using all observations
	allObsMedian = stats.nanmedian(uMatx,axis=0)

	# Monte carlo loop
	for j in range(0,iterMax):
		
		rowList = np.arange(0,numProfiles)
		
		# Random selection of timeslice loop
		for i in range(0,numProfiles):
			rowChoice = np.random.choice(rowList) # pick a random row
			choiceInd = np.where(rowList == rowChoice) # where is that random row?
			rowList = np.delete(rowList,choiceInd) # delete that row number so can't pick again
		
			addRow = uMatx[rowChoice,:] # pull row from umatx
			sampleNow[i,:,j] = addRow # add to median now
			medianNow[i,:,j] = stats.nanmedian(sampleNow[:,:,j],axis=0) # calculate column median with current sample
			
			diffFromFinal[i,:,j] = medianNow[i,:,j] - allObsMedian
					
			if i > 0:
				diffFromLast[i,:,j] = medianNow[i,:,j] - medianNow[i-1,:,j]
	

	
	#diff10MC, diff25MC, diff75MC, diff90MC = np.percentile(diffFromLast,[10,25, 75,90],axis=2)
	diff10MC, diff25MC, diff75MC, diff90MC = np.percentile(diffFromFinal,[10,25, 75,90],axis=2)
	
	#diff10MCmeanX = np.nanmean(diff10MC,axis=1)
	#diff25MCmeanX = np.nanmean(diff25MC,axis=1)
	#diff75MCmeanX = np.nanmean(diff75MC,axis=1)
	#diff90MCmeanX = np.nanmean(diff90MC,axis=1)
	diff10MCmeanX = stats.nanmedian(diff10MC,axis=1)
	diff25MCmeanX = stats.nanmedian(diff25MC,axis=1)
	diff75MCmeanX = stats.nanmedian(diff75MC,axis=1)
	diff90MCmeanX = stats.nanmedian(diff90MC,axis=1)	
	
	diffMedianAfterMC = stats.nanmedian(diffFromLast,axis=2)

	plt.fill_between(range(0,numProfiles),diff25MCmeanX,diff75MCmeanX,alpha=0.5,color='gray')
	
	for k in range(0,numProfiles):
		plt.scatter(k*np.ones(lenX),diffMedianAfterMC[k,:],linewidth=0,marker='.',edgecolor='',s=30, c=xInd, cmap=plt.cm.coolwarm)
		#plt.plot(k*np.ones(lenX),diff25MC[k,:],linewidth=2)
		#plt.plot(k*np.ones(lenX),diff75MC[k,:],linewidth=2)		
		
	plt.plot(range(0,numProfiles),diff10MCmeanX,linewidth=1,color='k')	
	plt.plot(range(0,numProfiles),diff25MCmeanX,linewidth=1,color='k')

	plt.plot(range(0,numProfiles),diff75MCmeanX,linewidth=1,color='k')
	plt.plot(range(0,numProfiles),diff90MCmeanX,linewidth=1,color='k')
	plt.plot(range(0,numProfiles),0.05*np.ones(numProfiles),linestyle='--',color='k')
	plt.plot(range(0,numProfiles),-0.05*np.ones(numProfiles),linestyle='--',color='k')
		
	plt.xlabel('Number of observations',fontsize=18)
	plt.ylabel('Velocity change relative to end [m d$^{-1}$]',fontsize=18)
	plt.xlim([0,numProfiles])
	plt.ylim([-.25,.25])
	plt.title('kennCL velocity stability, red = terminus',fontsize=20)
	plt.savefig('kennCLvelocityStabilityMarkData.pdf')
	plt.show()
	plt.close()
	
# Function that tests whether seasonal velocities are statistically different
def ranksumSeasonSpeed(transName):	
	if myDataToggle == 1:
		if annualToggle == 0:
			# <=96 day corr
			fnIn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/' + transName + '_evenlySpace_100m_profiles_sampled_2015-12-04.json'
			pathOut = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/'
			titleTag = ' <= 96 d'
		elif annualToggle == 1:
			# annual corr
			fnIn='/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/annualCorrelations/L8SamplerOutput/' + transName + '_evenlySpace_500m_profiles_sampled_2015-11-19.json'
			pathOut = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/annualCorrelations/L8SamplerOutput/'
			titleTag = ' annual'
	elif myDataToggle == 0:
		if annualToggle == 0:
			# corr from Mark's data
			fnIn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/markL8correlations/L8SamplerOutput/'  + transName + '_evenlySpace_500m_profiles_sampled_2015-11-17.json'
			pathOut = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/markL8correlations/L8SamplerOutput/'
			titleTag = ' Mark'
			
	
	#def readFromJson(fnIn):
	with open(fnIn) as data_file:
		data = json.load(data_file)

	shapefileName = data['profile_shapefile']
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
	
		# Plot point
		#Start season based
		# Identify start seasons
		if annualToggle == 0:
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
								
	# Replace -1 with nan
	nanInd = np.where(uMatx==-1)
	uWiNanInd = np.where(uWi==-1)
	uSpNanInd = np.where(uSp==-1)	
	uSuNanInd = np.where(uSu==-1)	
	uFaNanInd = np.where(uFa==-1)	
	
	uMatx[nanInd] = np.nan	
	uWi[uWiNanInd] = np.nan
	uSp[uSpNanInd] = np.nan	
	uSu[uSuNanInd] = np.nan	
	uFa[uFaNanInd] = np.nan			
	
	# Read in JSON data to see how many observations at each x
	obsNum = getNumDataFromJson(data)
	wiObs = obsNum[:,2]
	spObs = obsNum[:,3]	
	suObs = obsNum[:,4]
	faObs = obsNum[:,5]			
	
	for i in range(0,len(x)):		
		# Winter versus other seasons speeds
		if wiObs[i] >1 and spObs[i] >1:
			#s1, p1 = stats.ranksums(uWi[:,i],uSp[:,i])
			s1, p1 = stats.mannwhitneyu(uWi[:,i],uSp[:,i])
		else:
			s1 = np.nan
			p1 = np.nan			
		if wiObs[i] >1 and suObs[i] >1:
			s2, p2 = stats.mannwhitneyu(uWi[:,i],uSu[:,i])
		else:
			s2 = np.nan
			p2 = np.nan			
		if wiObs[i] >1 and faObs[i] >1:		
			s3, p3 = stats.mannwhitneyu(uWi[:,i],uFa[:,i])
		else:
			s3 = np.nan
			p3 = np.nan

		# Spring versus other seasons speeds
		if spObs[i] >1 and suObs[i] >1:		
			s4, p4 = stats.mannwhitneyu(uSp[:,i],uSu[:,i])
		else:
			s4 = np.nan
			p4  = np.nan
		if spObs[i] >1 and faObs[i] >1:			
			s5, p5 = stats.mannwhitneyu(uSp[:,i],uFa[:,i])
		else:
			s5 = np.nan
			p5 = np.nan
		# Summer versus other seasons speeds
		if suObs[i] >1 and faObs[i] >1:			
			s6, p6 = stats.mannwhitneyu(uSu[:,i],uFa[:,i])	
		else:
			s6 = np.nan
			p6 = np.nan		

		# Save out
		if i == 0:
			rsWiSp = np.array((s1, p1))
			rsWiSu = np.array((s2, p2))			
			rsWiFa = np.array((s3, p3))						
			rsSpSu = np.array((s4, p4))		
			rsSpFa = np.array((s5, p5))			
			rsSuFa = np.array((s6, p6))	
		if i > 0:
			rsWiSp = np.row_stack((rsWiSp,(s1,p1)))
			rsWiSu = np.row_stack((rsWiSu,(s2,p2)))			
			rsWiFa = np.row_stack((rsWiFa,(s3,p3)))			
			rsSpSu = np.row_stack((rsSpSu,(s4,p4)))			
			rsSpFa = np.row_stack((rsSpFa,(s5,p5)))			
			rsSuFa = np.row_stack((rsSuFa,(s6,p6)))						
	
	plt.subplot(311)						
	plt.plot(x/1e3,uWi.transpose(),color='c',marker='.',ms=10,lw=0,alpha=0.5)
	plt.plot(x/1e3,uSp.transpose(),color='y',marker='.',ms=10,lw=0,alpha=0.5)
	plt.plot(x/1e3,uSu.transpose(),color='g',marker='.',ms=10,lw=0,alpha=0.5)	
	plt.plot(x/1e3,uFa.transpose(),color='r',marker='.',ms=10,lw=0,alpha=0.5)	
	plt.ylim([0,2])
	plt.xlim([0,x.max()/1e3])
	plt.ylabel('Velocity [m d$^{-1}$]',fontsize=14)
	
	# make fake points for legend
	plt.plot(-1,-1,color='c',marker='.',ms=10,lw=0,alpha=0.5,label='winter')
	plt.plot(-1,-1,color='y',marker='.',ms=10,lw=0,alpha=0.5,label='spring')
	plt.plot(-1,-1,color='g',marker='.',ms=10,lw=0,alpha=0.5,label='summer')
	plt.plot(-1,-1,color='r',marker='.',ms=10,lw=0,alpha=0.5,label='fall')			
	
	plt.legend(loc='best',numpoints=1,prop={'size':8})
	
	plt.subplot(312)						
	plt.plot(x/1e3,wiObs,color='c',marker='.',ms=10)
	plt.plot(x/1e3,spObs,color='y',marker='.',ms=10)	
	plt.plot(x/1e3,suObs,color='g',marker='.',ms=10)	
	plt.plot(x/1e3,faObs,color='r',marker='.',ms=10)	
	plt.ylabel('# of observations',fontsize=14)
	plt.xlim([0,x.max()/1e3])
	plt.ylim([0,30])	
			
	plt.subplot(313)						
	plt.plot(x/1e3,rsWiSp[:,1],color='c',marker='v',ms=10,lw=0,alpha=0.5)
	plt.plot(x/1e3,rsWiSu[:,1],color='c',marker='o',ms=10,lw=0,alpha=0.5)
	plt.plot(x/1e3,rsWiFa[:,1],color='c',marker='s',ms=10,lw=0,alpha=0.5)	
	plt.plot(x/1e3,rsSpSu[:,1],color='y',marker='o',ms=10,lw=0,alpha=0.5)	
	plt.plot(x/1e3,rsSpFa[:,1],color='y',marker='s',ms=10,lw=0,alpha=0.5)		
	plt.plot(x/1e3,rsSuFa[:,1],color='g',marker='s',ms=10,lw=0,alpha=0.5)		
	plt.ylim([-0.01,0.1])
	plt.ylabel('mann-whitney p-val',fontsize=14)	
	plt.xlim([0,x.max()/1e3])
	
	# make fake points for legend
	plt.plot(-1,-1,color='c',marker='v',ms=10,lw=0,alpha=0.5,label='wi-sp')
	plt.plot(-1,-1,color='c',marker='o',ms=10,lw=0,alpha=0.5,label='wi-su')
	plt.plot(-1,-1,color='c',marker='s',ms=10,lw=0,alpha=0.5,label='wi-fa')
	plt.plot(-1,-1,color='y',marker='o',ms=10,lw=0,alpha=0.5,label='sp-su')
	plt.plot(-1,-1,color='y',marker='s',ms=10,lw=0,alpha=0.5,label='sp-fa')
	plt.plot(-1,-1,color='g',marker='s',ms=10,lw=0,alpha=0.5,label='su-fa')						
	plt.legend(loc='best',numpoints=1,prop={'size':8})
	plt.xlabel('Distance from headwall [km]',fontsize=14)
	
	plt.savefig('rega_statSigSpeedChange.pdf')	
	plt.show()
	plt.close()

def computeStatsAndPlot_retired_06dec2015(transName):

	gs = gridspec.GridSpec(3,1,height_ratios=[3,1,1])
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[1])
	ax3=plt.subplot(gs[2])

	if myDataToggle == 1:
		if annualToggle == 0:
			# <=96 day corr
			fnIn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/' + transName + '_evenlySpace_100m_profiles_sampled_2015-12-04.json'
			pathOut = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/'
			titleTag = ' <= 96 d'
		elif annualToggle == 1:
			# annual corr
			fnIn='/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/annualCorrelations/L8SamplerOutput/' + transName + '_evenlySpace_500m_profiles_sampled_2015-11-19.json'
			pathOut = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/annualCorrelations/L8SamplerOutput/'
			titleTag = ' annual'
	elif myDataToggle == 0:
		if annualToggle == 0:
			# corr from Mark's data
			fnIn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/markL8correlations/L8SamplerOutput/'  + transName + '_evenlySpace_500m_profiles_sampled_2015-11-17.json'
			pathOut = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/markL8correlations/L8SamplerOutput/'
			titleTag = ' Mark'
			
	#uAnnl = np.genfromtxt('/Users/wiar9509/Documents/generalScripts/database/regaMedianAnnualVelocity.csv',delimiter=',')
	
	#def readFromJson(fnIn):
	with open(fnIn) as data_file:
		data = json.load(data_file)

	# Read in JSON data to see how many observations at each x
	obsNum = getNumDataFromJson(data)

	shapefileName = data['profile_shapefile']
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
	
		# Plot point
		#Start season based
		# Identify start seasons
		if annualToggle == 0:
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
					
		else: # true if annual toggle == 1 (only plot `annual' profiles)
			seasonColor = 'grey'
	
		ax1.plot(x/1e3,u,linewidth=0,marker='.',color=seasonColor,alpha=0.5)
	
	# Replace -1 with nan
	nanInd = np.where(uMatx==-1)
	uMatx[nanInd] = np.nan	
	
	uMin = []
	uMed = []
	uMean = []
	uMax = []
	uRange = []
	uStd = []
	uSort = []
	
	for j in range(0,uMatx.shape[0]):
	
		# Now seasonal statistics. Think we could probably streamline and do these parts more together
	
		# WINTER CORRELATIONS
		if len(uWi) > 0:
			# Make -1 into nan
			nanIndWi = np.where(uWi==-1)
			uWi[nanIndWi] = np.nan
			numDimensions = uWi.shape
			if len(numDimensions) > 1: # was breaking for 0-d arrays before. This == 1 if 0d
				uWiSort=np.sort(uWi[:,j])
				dataIndWi = np.isnan(uWiSort) != True
				goodDataWi = uWiSort[dataIndWi]
				dataNumWi = len(goodDataWi)
		
				if dataNumWi > 0:
					lower10ind = np.round(dataNumWi*.1)
					lower25ind = np.round(dataNumWi*.25)
					medianInd = np.round(dataNumWi*.5)
					upper75ind = np.round(dataNumWi*.75)
					upper90ind = np.round(dataNumWi*.9)
					lower10now = goodDataWi[lower10ind]
					lower25now = goodDataWi[lower25ind]
					medianNow = goodDataWi[medianInd]
					if upper75ind > dataNumWi-1: # was indexing out of bounds
						#upper75now = goodDataWi[-1]
						upper75now = np.nan
					else:			
						upper75now = goodDataWi[upper75ind]
					if upper90ind > dataNumWi-1: # was indexing out of bounds
						#upper90now = goodDataWi[-1]
						upper90now = np.nan
					else:
						upper90now = goodDataWi[upper90ind]
				else:
					lower10now = np.nan 
					lower25now = np.nan
					medianNow = np.nan
					upper75now = np.nan
					upper90now = np.nan
		
				if j == 0: # if first iteration need to create array
					uWi10 = lower10now
					uWi25 = lower25now
					uWiMedian = medianNow
					uWi75 = upper75now
					uWi90 = upper90now
				else: # otherwise append		
					uWi10 = np.row_stack((uWi10,lower10now))
					uWi25 = np.row_stack((uWi25,lower25now))
					uWiMedian = np.row_stack((uWiMedian,medianNow))			
					uWi75 = np.row_stack((uWi75,upper75now))			
					uWi90 = np.row_stack((uWi90,upper90now))
				
			else: # if 0-d array then just one correlation, so median doesn't mean anything
				uWiMedian = np.empty(len(x))
				uWiMedian[:] = np.nan
		else:
			uWiMedian = np.empty(len(x))
			uWiMedian[:] = np.nan				

		# SPRING CORRELATIONS
		if len(uSp) > 0:
			# Make -1 into nan
			nanIndSp = np.where(uSp==-1)
			uSp[nanIndSp] = np.nan
			numDimensions = uSp.shape
			if len(numDimensions) > 1: # was breaking for 0-d arrays before. This == 1 if 0d
				uSpSort=np.sort(uSp[:,j])
				dataIndSp = np.isnan(uSpSort) != True
				goodDataSp = uSpSort[dataIndSp]
				dataNumSp = len(goodDataSp)
		
				if dataNumSp > 0:
					lower10ind = np.round(dataNumSp*.1)
					lower25ind = np.round(dataNumSp*.25)
					medianInd = np.round(dataNumSp*.5)
					upper75ind = np.round(dataNumSp*.75)
					upper90ind = np.round(dataNumSp*.9)
					lower10now = goodDataSp[lower10ind]
					lower25now = goodDataSp[lower25ind]
					medianNow = goodDataSp[medianInd]
					if upper75ind > dataNumSp-1: # was indexing out of bounds
						#upper75now = goodDataSp[-1]
						upper75now = np.nan
					else:			
						upper75now = goodDataSp[upper75ind]
					if upper90ind > dataNumSp-1: # was indexing out of bounds
						upper90now = goodDataSp[-1]
						upper90now = np.nan
					else:
						upper90now = goodDataSp[upper90ind]
				else:
					lower10now = np.nan 
					lower25now = np.nan
					medianNow = np.nan
					upper75now = np.nan
					upper90now = np.nan
		
				if j == 0: # if first iteration need to create array
					uSp10 = lower10now
					uSp25 = lower25now
					uSpMedian = medianNow
					uSp75 = upper75now
					uSp90 = upper90now
				else: # otherwise append		
					uSp10 = np.row_stack((uSp10,lower10now))
					uSp25 = np.row_stack((uSp25,lower25now))
					uSpMedian = np.row_stack((uSpMedian,medianNow))			
					uSp75 = np.row_stack((uSp75,upper75now))			
					uSp90 = np.row_stack((uSp90,upper90now))	

			else: # if 0-d array then just one correlation, so median doesn't mean anything
				uSpMedian = np.empty(len(x))
				uSpMedian[:] = np.nan
		else:
			uSpMedian = np.empty(len(x))
			uSpMedian[:] = np.nan
			
		# SUMMER CORRELATIONS
		if len(uSu) > 0:
			# Make -1 into nan
			nanIndSu = np.where(uSu==-1)
			uSu[nanIndSu] = np.nan
			numDimensions = uSu.shape
			if len(numDimensions) > 1: # was breaking for 0-d arrays before. This == 1 if 0d		
				uSuSort=np.sort(uSu[:,j])
				dataIndSu = np.isnan(uSuSort) != True
				goodDataSu = uSuSort[dataIndSu]
				dataNumSu = len(goodDataSu)
		
				if dataNumSu > 0:
					lower10ind = np.round(dataNumSu*.1)
					lower25ind = np.round(dataNumSu*.25)
					medianInd = np.round(dataNumSu*.5)
					upper75ind = np.round(dataNumSu*.75)
					upper90ind = np.round(dataNumSu*.9)
					lower10now = goodDataSu[lower10ind]
					lower25now = goodDataSu[lower25ind]
					medianNow = goodDataSu[medianInd]
					if upper75ind > dataNumSu-1: # was indexing out of bounds
						#upper75now = goodDataSu[-1]
						upper75now = np.nan
					else:			
						upper75now = goodDataSu[upper75ind]
					if upper90ind > dataNumSu-1: # was indexing out of bounds
						upper90now = goodDataSu[-1]
						upper90now = np.nan
					else:
						upper90now = goodDataSu[upper90ind]
				else:
					lower10now = np.nan 
					lower25now = np.nan
					medianNow = np.nan
					upper75now = np.nan
					upper90now = np.nan
		
				if j == 0: # if first iteration need to create array
					uSu10 = lower10now
					uSu25 = lower25now
					uSuMedian = medianNow
					uSu75 = upper75now
					uSu90 = upper90now
				else: # otherwise append		
					uSu10 = np.row_stack((uSu10,lower10now))
					uSu25 = np.row_stack((uSu25,lower25now))
					uSuMedian = np.row_stack((uSuMedian,medianNow))			
					uSu75 = np.row_stack((uSu75,upper75now))			
					uSu90 = np.row_stack((uSu90,upper90now))	

			else: # if 0-d array then just one correlation, so median doesn't mean anything
				uSuMedian = np.empty(len(x))
				uSuMedian[:] = np.nan
		else:
			uSuMedian = np.empty(len(x))
			uSuMedian[:] = np.nan
										
		# FALL CORRELATIONS
		if len(uFa) > 0:
			# Make -1 into nan
			nanIndFa = np.where(uFa==-1)
			uFa[nanIndFa] = np.nan
			numDimensions = uFa.shape
			if len(numDimensions) > 1: # was breaking for 0-d arrays before. This == 1 if 0d	
				uFaSort=np.sort(uFa[:,j])
				dataIndFa = np.isnan(uFaSort) != True
				goodDataFa = uFaSort[dataIndFa]
				dataNumFa = len(goodDataFa)
		
				if dataNumFa > 0:
					lower10ind = np.round(dataNumFa*.1)
					lower25ind = np.round(dataNumFa*.25)
					medianInd = np.round(dataNumFa*.5)
					upper75ind = np.round(dataNumFa*.75)
					upper90ind = np.round(dataNumFa*.9)
					lower10now = goodDataFa[lower10ind]
					lower25now = goodDataFa[lower25ind]
					medianNow = goodDataFa[medianInd]
					if upper75ind > dataNumFa-1: # was indexing out of bounds
						#upper75now = goodDataFa[-1]
						upper75now = np.nan
					else:			
						upper75now = goodDataFa[upper75ind]
					if upper90ind > dataNumFa-1: # was indexing out of bounds
						#upper90now = goodDataFa[-1]
						upper90now = np.nan
					else:
						upper90now = goodDataFa[upper90ind]
				else:
					lower10now = np.nan 
					lower25now = np.nan
					medianNow = np.nan
					upper75now = np.nan
					upper90now = np.nan
				if j == 0: # if first iteration need to create array
					uFa10 = lower10now
					uFa25 = lower25now
					uFaMedian = medianNow
					uFa75 = upper75now
					uFa90 = upper90now
				else: # otherwise append		
					uFa10 = np.row_stack((uFa10,lower10now))
					uFa25 = np.row_stack((uFa25,lower25now))
					uFaMedian = np.row_stack((uFaMedian,medianNow))			
					uFa75 = np.row_stack((uFa75,upper75now))			
					uFa90 = np.row_stack((uFa90,upper90now))	

			else: # if 0-d array then just one correlation, so median doesn't mean anything
				uFaMedian = np.empty(len(x))
				uFaMedian[:] = np.nan
		else:
			uFaMedian = np.empty(len(x))
			uFaMedian[:] = np.nan
							
		# Calculating statistics regardless of season
		if j == 0: # need to create if first iteration
			uMin=np.nanmin(uMatx[j,:])
			uMean=np.nanmean(uMatx[j,:])
			uMax=np.nanmax(uMatx[j,:])
			uRange=uMax-uMin
			uStd=np.nanstd(uMatx[j,:])
			uSort=np.sort(uMatx[j,:])
			dataInd = np.isnan(uSort) != True
			goodData = uSort[dataInd]
			dataNum = len(goodData)
			if dataNum > 0: 
				lower10ind = np.round(dataNum*.1)
				lower25ind = np.round(dataNum*.25)
				medianInd = np.round(dataNum*.5)
				upper75ind = np.round(dataNum*.75)
				upper90ind = np.round(dataNum*.9)
				lower10 = goodData[lower10ind]
				lower25 = goodData[lower25ind]
				median = goodData[medianInd]
				if upper75ind > dataNum-1: # was indexing out of bounds
					upper75 = goodData[-1]
				else:			
					upper75 = goodData[upper75ind]
				if upper90ind > dataNum-1: # was indexing out of bounds
					upper90 = goodData[-1]
				else:
					upper90 = goodData[upper90ind]
			else:
				lower10 = np.nan
				lower25 = np.nan
				median = np.nan
				upper75 = np.nan
				upper90 = np.nan
		
		else:
			uMin=np.row_stack((uMin,np.nanmin(uMatx[j,:])))
			uMean=np.row_stack((uMean,np.nanmean(uMatx[j,:])))		
			uMax=np.row_stack((uMax,np.nanmax(uMatx[j,:])))		
			uRange=np.row_stack((uRange, uMax[j]-uMin[j] ))		
			uStd=np.row_stack((uStd, np.nanstd(uMatx[j,:]) ))	
			uSort=np.sort(uMatx[j,:])
			dataInd = np.isnan(uSort) != True
			goodData = uSort[dataInd]
			dataNum = len(goodData)
			if dataNum > 0:		
				lower10ind = np.round(dataNum*.1)
				lower25ind = np.round(dataNum*.25)
				medianInd = np.round(dataNum*.5)
				upper75ind = np.round(dataNum*.75)
				upper90ind = np.round(dataNum*.9)	
				lower10=np.row_stack((lower10,goodData[lower10ind]))
				lower25=np.row_stack((lower25,goodData[lower25ind]))
				median=np.row_stack((median,goodData[medianInd]))
				if upper75ind > dataNum-1:
					upper75	=np.row_stack((upper75,goodData[-1]))
				else:	
					upper75=np.row_stack((upper75,goodData[upper75ind]))
				if upper90ind > dataNum-1:
					upper90=np.row_stack((upper90,goodData[-1]))											
				else:
					upper90=np.row_stack((upper90,goodData[upper90ind]))	
			else:
				lower10=np.row_stack((lower10,np.nan))
				lower25=np.row_stack((lower25,np.nan))					
				median=np.row_stack((median,np.nan))					
				upper75=np.row_stack((upper75,np.nan))					
				upper90=np.row_stack((upper90,np.nan))				# append if not first iteration

	saveData = np.vstack((x,obsNum[:,1],lower25.flatten(),median.flatten(),upper75.flatten(),obsNum[:,2],uWiMedian.flatten(),obsNum[:,3],uSpMedian.flatten(),obsNum[:,4],uSuMedian.flatten(),obsNum[:,5],uFaMedian.flatten()))
	saveDataT = saveData.transpose()
	np.savetxt(pathOut + shapefileName[:-4] +'_velocityStats.csv',saveDataT,delimiter=',')
			 
					

	
	## GET ELEVATION PROFILE
	#elevData = np.genfromtxt(elevationProfileFolder + shapefileName[:-4] + '_elevationProfile.csv',delimiter=',')
	elevData = np.genfromtxt(elevationProfileFolder + transName + '_elevationProfile.csv',delimiter=',')
		
	
		
	## PLOTTING			
	# These axes generated way up high at start of script	

	
	if annualToggle == 0:
		# plot seasonal medians
		ax1.plot(x/1e3,uWiMedian,linewidth=2,color='c')		
		ax1.plot(x/1e3,uSpMedian,linewidth=2,color='y')		
		ax1.plot(x/1e3,uSuMedian,linewidth=2,color='g')
		ax1.plot(x/1e3,uFaMedian,linewidth=2,color='r')		

	ax2.plot(obsNum[:,0]/1e3,obsNum[:,1],marker='.',color='k')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,2],marker='.',color='c')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,3],marker='.',color='y')	
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,4],marker='.',color='g')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,5],marker='.',color='r')		
					
	ax1.plot(x/1e3,lower25,linewidth=1,color='k')
	ax1.plot(x/1e3,median,linewidth=3,color='k',label='Median: all data')
	ax1.plot(x/1e3,upper75,linewidth=1,color='k',label='25th-75th %: all data')
	
	#ax2 = ax1.twinx()

	ax3.plot(elevData[:,0]/1e3,elevData[:,1],color='gray')

	ax3.set_xlabel('Distance from headwall [km]',fontsize=16)
	ax1.set_ylabel('Velocity [m d$^{-1}$]',fontsize=16)
	ax1.set_title(shapefileName + titleTag,fontsize=18)
	
	if annualToggle == 0:
		ax1.plot(-1,-1,marker='.',color='c',linewidth=0,label='Winter')
		ax1.plot(-1,-1,marker='.',color='y',linewidth=0,label='Spring')
		ax1.plot(-1,-1,marker='.',color='g',linewidth=0,label='Summer')
		ax1.plot(-1,-1,marker='.',color='r',linewidth=0,label='Fall')

		
	#ax3.set_ylim([np.nanmin(elevData[:,1]),	np.nanmax(elevData[:,1])])
	ax2.set_ylabel('# Observations',fontsize=18)
	ax2.yaxis.set_label_position("right")
	ax3.set_ylabel('Elevation [m]',fontsize=18)
	#ax2.set_yticks(np.arange( round(elevData[:,1].min()/500)*500,round(elevData[:,1].max()/500)*500+500,500))
	ax1.set_ylim([0,np.nanmin((2.5,np.nanmax(upper75)))])
	ax1.set_xlim([0,x.max()/1e3])
	ax1lims = ax1.xaxis.get_view_interval()
	ax3.set_xlim([ax1lims[0],ax1lims[1]])
	ax2.set_xlim([ax1lims[0],ax1lims[1]])
	ax1.legend(loc='best',prop={'size':10},numpoints=1)
	plt.draw()
	plt.savefig(pathOut + shapefileName[0:-4] + '_profilePlot_wSeasonalMedian_wElevation_wDataNum.pdf')
	#plt.close()
	#plt.show()
	plt.close()		


# Function to modify json file to add keys for data separated annually
def separateDataInJson(jsonData):
	



# Function to calculate season velocity percentiles		
# Takes in seasonal velocity (here called uSu b/c just copied from below)
# Ouputs dataOut = (x, u10, u25, u50, u75, u90)
# Something going wrong at 75th and 90th percentiles where just all nan. Think might have to do with treatment of nan
def calculateSeasonalVelocityStats(uSeason):		

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


def computeStatsAndPlot_v2(transName):

	if myDataToggle == 1:
		if annualToggle == 0:
			# <=96 day corr
			fnIn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/' + transName + '_evenlySpace_100m_profiles_sampled_2015-12-04.json'
			pathOut = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/'
			titleTag = ' <= 96 d'
		elif annualToggle == 1:
			# annual corr
			fnIn='/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/annualCorrelations/L8SamplerOutput/' + transName + '_evenlySpace_500m_profiles_sampled_2015-11-19.json'
			pathOut = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/annualCorrelations/L8SamplerOutput/'
			titleTag = ' annual'
	elif myDataToggle == 0:
		if annualToggle == 0:
			# corr from Mark's data
			fnIn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/markL8correlations/L8SamplerOutput/'  + transName + '_evenlySpace_500m_profiles_sampled_2015-11-17.json'
			pathOut = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/markL8correlations/L8SamplerOutput/'
			titleTag = ' Mark'
			
	#uAnnl = np.genfromtxt('/Users/wiar9509/Documents/generalScripts/database/regaMedianAnnualVelocity.csv',delimiter=',')
	
	#def readFromJson(fnIn):
	with open(fnIn) as data_file:
		data = json.load(data_file)

	# Read in JSON data to see how many observations at each x
	obsNum = getNumDataFromJson(data)

	shapefileName = data['profile_shapefile']
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
	
		# Plot point
		#Start season based
		# Identify start seasons
		if annualToggle == 0:
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
					
		else: # true if annual toggle == 1 (only plot `annual' profiles)
			seasonColor = 'grey'
	
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

	if 0:
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
	np.savetxt(pathOut + shapefileName[:-4] +'_velocityStats.csv',saveDataT,delimiter=',')
			 
						
	## GET ELEVATION PROFILE
	#elevData = np.genfromtxt(elevationProfileFolder + shapefileName[:-4] + '_elevationProfile.csv',delimiter=',')
	elevData = np.genfromtxt(elevationProfileFolder + transName + '_elevationProfile.csv',delimiter=',')
		
	
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,1],marker='.',color='k')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,2],marker='.',color='c')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,3],marker='.',color='y')	
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,4],marker='.',color='g')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,5],marker='.',color='r')		
					
	ax3.plot(elevData[:,0]/1e3,elevData[:,1],color='gray')

	ax3.set_xlabel('Distance from headwall [km]',fontsize=16)
	ax1.set_ylabel('Velocity [m d$^{-1}$]',fontsize=16)
	ax2.set_title(shapefileName + titleTag,fontsize=18)
	
	ax1.plot(-1,-1,marker='.',color='c',linewidth=0,label='Winter')
	ax1.plot(-1,-1,marker='.',color='y',linewidth=0,label='Spring')
	ax1.plot(-1,-1,marker='.',color='g',linewidth=0,label='Summer')
	ax1.plot(-1,-1,marker='.',color='r',linewidth=0,label='Fall')

	
	#ax3.set_ylim([np.nanmin(elevData[:,1]),	np.nanmax(elevData[:,1])])
	ax2.set_ylabel('# Observations',fontsize=18)
	ax2.yaxis.set_label_position("right")
	ax3.set_ylabel('Elevation [m]',fontsize=18)
	#ax2.set_yticks(np.arange( round(elevData[:,1].min()/500)*500,round(elevData[:,1].max()/500)*500+500,500))
	ax1.set_ylim([0,np.nanmin((2.5,np.nanmax(uAllStats[4,:])))])
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
#	plt.savefig(pathOut + shapefileName[0:-4] + '_profilePlot_wSeasonalMedian_wElevation_wDataNum_wDiffFromMean.pdf')
	plt.show()
	plt.close()		
		
### FOR ITERATING OVER EACH TRANSECT
ds = ogr.Open(shapefileIn)
lyr = ds.GetLayer(0)
numLines = lyr.GetFeatureCount()


### FOR ITERATING OVER EACH TRANSECT
ds = ogr.Open(shapefileIn)
lyr = ds.GetLayer(0)
numLines = lyr.GetFeatureCount()

for i in range(0,numLines):
	line = lyr.GetFeature(i)
	transName = line.GetField(1)
	computeStatsAndPlot_v2(transName)
#transName='rootCL'
		
		