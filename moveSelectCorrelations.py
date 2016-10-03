'''
Script to use a .csv list to identify select correlations to move
02 oct 2016
'''

### MODULES ###

import os
import glob
import numpy as np
import datetime as dt
import time

### FUNCTIONS ###

def sinceEpoch(date): # returns seconds since epoch
	'''
	marks function to get seconds in a year for use in toYearFraction()
	'''
	return time.mktime(date.timetuple())

def toYearFraction(date):
	'''
	Turns a datetime date object into a fractional year date (float)
	Input: date (datetime date object)
	Output: fractional year date (float)
	Written by Mark Fahnestock
	
	from Mark:
	########################################################
	# following from accepted answer here: http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years
	# but I import datetime package as dt, not datetime.datetime as dt, so this is modfied to reflect that
	########################################################
	'''
	s = sinceEpoch

	year = date.year
	startOfThisYear = dt.datetime(year=year, month=1, day=1)
	startOfNextYear = dt.datetime(year=year+1, month=1, day=1)

	yearElapsed = s(date) - s(startOfThisYear)
	yearDuration = s(startOfNextYear) - s(startOfThisYear)
	fraction = yearElapsed/yearDuration

	return date.year + fraction

def doyToDate(year,doy):
	'''
	Simple function to turn a year & day of year into a date
	Input: year = year (int), doy = doy (int)
	Output: dateObject = date (datetime date object)
	'''
	dateObject = dt.date.fromordinal(dt.date(year-1,12,31).toordinal()+doy)
	return dateObject


def parseL8filename(l8fn):
	'''
	Parse Landsat 8 filename
	Input : landsat 8 filename (str)
	Output : data.viewkeys() = ['path','row','year','doy']
	'''
	path = l8fn[3:6]
	row = l8fn[6:9]
	year = l8fn[9:13]
	doy = l8fn[13:16]
	
	return {'path':path,'row':row,'year':year,'doy':doy}

def parseCorrFilename(corrFn):
	'''
	Read a correlation filename and get out pertinent details
	Input: correlation filename (no path)
	Output: dictionary with correlation details
	'''
	image1fn,image2fn,yearMidDoy,dtDays,source,target,inc,trash,trash,hp,trash = corrFn.split('_')
	midYear,midDoy = yearMidDoy.split('-')
	midDate = doyToDate(int(midYear),int(midDoy))
	midYearFrac = toYearFraction(midDate)
	image1dict = parseL8filename(image1fn)
	image2dict = parseL8filename(image2fn)

	corrDict = {'image1info':image1dict,'image2info':image2dict,'path':image1dict['path'],'row':image1dict['row'],
				'dtDays':dtDays,'midYearFrac':midYearFrac,'sourceSize':source,'targetSize':target,'increment':inc,
				'image1name':image1fn,'image2name':image2fn,'midYearDoy':yearMidDoy,'gaussian':hp}

	return corrDict
	
def makeCorrFilepath(corrFolder,corrFn):
	'''
	Make filepath to correlation folder for a specified filename
	'''
	corrDict = parseCorrFilename(corrFn)
	path = corrDict['path']
	row = corrDict['row']
	fn1 = corrDict['image1name']
	fn2 = corrDict['image2name']
	midYearDoy = corrDict['midYearDoy']
	dtDays = corrDict['dtDays']
	source = corrDict['sourceSize']
	target = corrDict['targetSize']
	inc = corrDict['increment']
	hp = corrDict['gaussian']

	pathRowFolder = path + '_' + row 
	imagePairFolder = fn1 + '_' + fn2 + '_' + midYearDoy + '_' + dtDays + '_hp'
	correlationFolder = fn1 + '_' + fn2 + '_' + midYearDoy + '_' + dtDays + '_' + str(int(float(hp))) + '_' + source + '_' + target + '_' + inc

	corrFilepath = corrFolder + pathRowFolder + '/' + imagePairFolder + '/' + correlationFolder + '/'

	outDict = {'corrFilepath':corrFilepath,'masterFolder':corrFolder,'pathRowFolder':pathRowFolder,'imagePairFolder':imagePairFolder,'correlationFolder':correlationFolder}
	return outDict

### USER DEFINED ###

# list of good corelations
corrList = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/selectWinterCorrelationList.csv'

# where do all correlations sit?
corrFolder = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/'

# where to write to?
outputFolder = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/'

### ANALYSIS ###

# read in select correlations list
corrNames = np.genfromtxt(corrList, delimiter=',',skip_header=1,dtype='str')[:,0]

# loop over select correlations
for corrName in corrNames:
	
	# Generate filepath for correlation of interest
	filepathDict = makeCorrFilepath(corrFolder,corrName)
	filepath = filepathDict['corrFilepath']
	
	# Find filenames for files within this folder to move
	vvFn = glob.glob(filepath + "*_vv.tif")[0]
	vxFn = glob.glob(filepath + "*_vx.tif")[0]
	vyFn = glob.glob(filepath + "*_vy.tif")[0]
	log10Fn = glob.glob(filepath + "*_log10.tif")[0]
	corrStrengthFn = glob.glob(filepath + "*_corr.tif")[0]
	delcorrFn = glob.glob(filepath + "*_delcorr.tif")[0]
	
	filesToCopy = [vvFn, vxFn, vyFn, corrStrengthFn, delcorrFn, log10Fn]
	
	# Iterate over files to copy
	for copyFile in filesToCopy:
		shortFn = os.path.split(copyFile)[1]
		outName = outputFolder + shortFn # generate output name
		if not os.path.isfile(outName): # don't move if it's already there
			print "Moving files for " + copyFile
			copyCmd = "cp " + copyFile + " " + outName
			os.system(copyCmd) # copy over


### MOVE LOG10files
if 0:
	pathRowFolder = '064_017/'
	foldersInPathRow=glob.glob(corrFolder + pathRowFolder + "LC*")
	log10folder = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/log10/'

	for folder in foldersInPathRow:
		corrFolders = glob.glob(folder + '/*')

		for i in range(len(corrFolders)):
			folderNow = corrFolders[i] + '/'
			log10path = glob.glob(folderNow + '*log10.tif')[0]
			log10name = os.path.split(log10path[0])[1]
			outName = log10folder + log10name
		
			if not os.path.isfile(outName):
	#			print "Will copy file: " , log10path
				print "Copying file ", log10path
				copyCmd = "cp " + log10path + " " + outName
				os.system(copyCmd) # copy over
		
		



