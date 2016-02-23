import json


testData = {
	'transName' : 'kennCL',
	'year' : 2015,
	'xDist': (0,10,20,30,40,50)
	}
	
json_str = json.dumps(testData,separators=(',',':'))

test_json = json.loads(json_str)


### SCRIPT STARTING HERE

import json
import numpy as np
import matplotlib.pyplot as plt

testJsonFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput/kennCL_evenlySpace_100m_profiles_sampled_2015-12-04.json'
yearStartDoy = 274 # doy 274 = oct 1. will put things from post oct 1 in with next year's observations.
recordYears = [2013,2014,2015,2016] # years in which data exists

with open(testJsonFn) as data_file:
	jsonData = json.load(data_file)

numProfiles = len(jsonData["pts"]["mid_dec_year"]) # how many profiles are there?


for yearNow in recordYears: # iterate over years of record
	print "Record year: " + str(yearNow)
	x = jsonData["sample_pts_frontdist"] # x-coordinate
	jsonData[str(yearNow) + 'profiles'] = None # clear for each year
	uList = [] # clear for each year
	
	for i in range(0,numProfiles): # iterate over profiles
		midDoy = jsonData["pts"]["mid_dec_year"][i] # midpoint of correlation date
		if int(np.floor(midDoy)) == yearNow: # true if this profile in year of interest
			print midDoy
			if jsonData["profiles"][i]["mid_dec_year"] == midDoy: # this should always be true, but just checking that referencing the right file
				uNow = jsonData["profiles"][i]["speed"]
				uList.append(uNow)
				
		if i == numProfiles-1: # true if gone through all profiles
			#print {str(yearNow) + 'profiles':uList}
			jsonData.update({str(yearNow) + 'profiles':uList})

					
if yearNow == recordYears[-1] and i == numProfiles-1: # true if gone through all profiles and all years
	for year in recordYears:
		samplesInYear = len(jsonData[str(year) + 'profiles']) # how many samples w/in year
		print "Current year: " + str(year) + " has # samples: " + str(samplesInYear)
		for j in range(0,samplesInYear):
			if year == 2013:
				col = 'b'
				aVal = 0.75
			elif year == 2014:
				col = 'c'
				aVal = 0.4
			elif year == 2015:
				col = 'm'
				aVal = 0.3
			elif year == 2016:
				col = 'r'
				aVal = 0.2
				
			plt.plot(x,jsonData[str(year) + 'profiles'][j],color=col,lw=0,marker='.',alpha=aVal)
	
	y1 = plt.plot(-1,-1,color='b',lw=0,marker='.',alpha=0.75,label='2013')
	y2 = plt.plot(-1,-1,color='c',lw=0,marker='.',alpha=0.75,label='2014')	
	y3 = plt.plot(-1,-1,color='m',lw=0,marker='.',alpha=0.75,label='2015')		
#	y4 = plt.plot(-1,-1,color='m',lw=0,marker='.',alpha=0.75,label='2015')			
	plt.legend(numpoints=1)
	plt.xlim(0,np.max(x))
	plt.xlabel('Down-glacier distance [m]',fontsize=16)
	plt.ylabel('Horizontal velocity [m d$^{-1}$]',fontsize=16)
	plt.ylim(0,2)
	plt.show()
	plt.close()
	
				
