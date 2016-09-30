'''
SCRIPT TO ITERATE IMAGE CORRELATION
29 SEP 2016
WILLIAM ARMSTRONG
'''

### IMPORT MODULES ###

import os
import glob
import datetime as dt
#%cd /Users/anderson/Desktop/ARMSTRONG/git/wrangells
#import wrangellsGeospatialTools as wgt
##cd /Users/anderson/

### FUNCTIONS ###

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


### USER INPUTS ###

folderPath='/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/63_17/band8/'
outputDirectory='/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/' # want to make functionality to specify output directory, but avoiding for now -- wha 29sep16

sourceChipSize = 24
targetChipSize = 40
increment = 10
maxDaysBetween = 400

### PROCESSING ###

masters = glob.glob(folderPath+"*_B8.TIF")

for master in masters:
	slaves = glob.glob(folderPath+"*_B8.TIF")
	
	for slave in slaves:
		if master != slave:
			
			trash,trash,mShortFn = master.rpartition('/') # parse into filename and directory
			trash,trash,sShortFn = slave.rpartition('/') # parse into filename and directory
			
			doy1, year1 = int(parseL8filename(mShortFn)['doy']) ,int(parseL8filename(mShortFn)['year'])
			doy2, year2 = int(parseL8filename(sShortFn)['doy']) ,int(parseL8filename(sShortFn)['year'])		

			dt_1=dt.datetime.fromordinal(dt.date(year1-1,12,31).toordinal()+doy1)
			dt_2=dt.datetime.fromordinal(dt.date(year2-1,12,31).toordinal()+doy2)
			timedelta=dt_2-dt_1
			delts_days=timedelta.days
			
			if delts_days > 0 and delts_days <= maxDaysBetween:
				# want to write in check to not repeat correlations that have already been done, but not doing for now -- wha 29sep16
				#outDirName = folderPath+out_name_base+'_'+str(int(args.gfilt_sigma))+'_'+str(int(args.half_source_chip))+'_'+str(int(args.half_target_chip))+'_'+str(int(args.inc))+'/'
				commandString = 'python wpycorr_v1.10_whaModified_092515.py ' + mShortFn + " " + sShortFn + " -vf -sf -v -imgdir " + folderPath + " -nsidc_out_names -use_hp -gfilt_sigma 3.0 -log10 -inc " + str(increment) + " -half_source_chip " + str(sourceChipSize) + " -half_target_chip " + str(targetChipSize)
				
				%cd $folderPath
				os.system(commandString)
				