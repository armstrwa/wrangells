'''
Need to run the first part of bestFullyAutomated... to get seasonal velocity dict
'''

def splitSummerWinterToJson(jsonData,summerStartEnd=((150,250))):
	'''
	Split up correlations into winter/summer categories.
	Inputs: jsonData = mark L8 sampler-like velocity dict
	Inputs: summerStartEnd = ((150,250)) = DOYs for start/end of summer season for categorizing
	Outputs: seaonsalVelocityDict.viewkeys() = {'summerCorr','summerSpeeds','winterCorr','winterSpeeds'}
	Outputs: where summer/winter speeds are velocity matrices (len(space) cols and len(time) rows), and sumer/winter corr are lists of profile dicts.
	'''
	# Classify as summer/winter correlations and clean bad data
	winterCorr = [] # intializing
	summerCorr = []
	winterSwathSpeeds = np.nan * np.ones((len(jsonData['profiles'][0]['speed']))) # intializing arrays for stats
	winterRawSpeeds = np.nan * np.ones((len(jsonData['profiles'][0]['speed']))) # intializing arrays for stats
	winterTempFilterSpeeds = np.nan * np.ones((len(jsonData['profiles'][0]['speed']))) # intializing arrays for stats
	summerSwathSpeeds = np.nan * np.ones((len(jsonData['profiles'][0]['speed'])))
	summerTempFilterSpeeds = np.nan * np.ones((len(jsonData['profiles'][0]['speed'])))
	summerRawSpeeds = np.nan * np.ones((len(jsonData['profiles'][0]['speed']))) # intializing arrays for stats
	
	nProfiles =  len(jsonData['profiles'])
	# Iterate over profiles
	for i in range(0,nProfiles):
		profileNow = jsonData['profiles'][i] # get current profile
		midYearFrac = profileNow['mid_dec_year'] # get date (fractional year form)
		midDoy = wgt.yearFracToDoy(midYearFrac) # convert fractional year to DOY
		uCleanNow = profileNow['speedMeanClean'] # get swath-filtered velocity data
		uTempFilterNow = profileNow['speedTemporalFilter'] # temporally filtered velocity
		uRawNow  = profileNow['speed'] # raw velocity data
		
		# If it's winter
		if midDoy < summerStartEnd[0] or midDoy > summerStartEnd[1]:
			winterCorr.append(profileNow) # add correlation to collection of winter correlations
			winterSwathSpeeds = np.vstack((winterSwathSpeeds,uCleanNow)) # add speed to matrix of winter speeds
			winterTempFilterSpeeds = np.vstack((winterTempFilterSpeeds,uTempFilterNow)) # add temporally filtered speed to matrix of winter speeds
			winterRawSpeeds = np.vstack((winterRawSpeeds,uRawNow)) # add raw speed to matrix of winter speeds
		# If it's summer
		else:
			summerCorr.append(profileNow)
			summerSwathSpeeds = np.vstack((summerSwathSpeeds,uCleanNow))	
			summerTempFilterSpeeds = np.vstack((summerTempFilterSpeeds,uTempFilterNow))
			summerRawSpeeds = np.vstack((summerRawSpeeds,uRawNow))
			
	# Make dict for output	
	transName = jsonData['transect_name']
#	seaonsalVelocityDict = {'transName':transName,'summerCorr':summerCorr,'summerSwathSpeeds':summerSwathSpeeds.tolist(),'summerTempFilterSpeeds':summerTempFilterSpeeds.tolist(),'winterCorr':winterCorr,'winterSwathSpeeds':winterSwathSpeeds.tolist(),'winterTempFilterSpeeds':winterTempFilterSpeeds.tolist(),'summerRawSpeeds':summerRawSpeeds.tolist(),'winterRawSpeeds':winterRawSpeeds.tolist()}
	seaonsalVelocityDict = {'transName':transName,'summerSwathSpeeds':summerSwathSpeeds.tolist(),'summerTempFilterSpeeds':summerTempFilterSpeeds.tolist(),'winterSwathSpeeds':winterSwathSpeeds.tolist(),'winterTempFilterSpeeds':winterTempFilterSpeeds.tolist(),'summerRawSpeeds':summerRawSpeeds.tolist(),'winterRawSpeeds':winterRawSpeeds.tolist()}
	return seaonsalVelocityDict


keyIter = seaonsalVelocityDict.iterkeys()
for i in range(0,len(seaonsalVelocityDict)):
	now = keyIter.next()
	print now, type(seaonsalVelocityDict[now])
	
wgt.writeJsonFile(seaonsalVelocityDict,'test.json')