'''
Image coverage plot
11 oct 2016
wha
'''

import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import datetime as dt
import time
import pandas as pd
import matplotlib.dates as mdates

def getCorrelationInfoFromFn(corrFnFullFilepath):
	'''
	Get information on image cross-correlation from correlation filename. To be used for data storage.
	Input: correlation full filepath (str) expecting format: 'LC80640172013342LGN0_LC80640172014185LGN0_2014-081_208_24_40_10_vv_filt_corr0.15_delcorr0.05.tif'
	Output: profile_dict.viewkeys() = ['delts_days','series','start_dec_year','stop_dec_year','mid_dec_year','path','row','infile','indir'
	'''
	
	filepath,part,filename=corrFnFullFilepath.rpartition('/') # parse into filename and directory
	
	# Parsing filename
	try: # if filtered
		image1,image2,yearDoyMedian,daysBtwn,sourceSize,targetSize,step,vType,filtStatus,corr,delcorr,spatialReference=filename.split('_')	
	except ValueError: # if not filtered, or if just wrong size..
		image1,image2,yearDoyMedian,daysBtwn,sourceSize,targetSize,step,vType=filename.split('_')[0:8]
		filtStatus = None
		corr = None
		delcorr = None

	series_tag=None # this just to be consistent with Mark's old scheme			
	doy1, year1 = int(parseL8filename(image1)['doy']) ,int(parseL8filename(image1)['year'])
	doy2, year2 = int(parseL8filename(image2)['doy']) ,int(parseL8filename(image2)['year'])
	path,row = int(parseL8filename(image1)['path']) ,int(parseL8filename(image1)['row'])			
# 	doy1, year1 = int(wgt.parseL8filename(image1)['doy']) ,int(wgt.parseL8filename(image1)['year'])
# 	doy2, year2 = int(wgt.parseL8filename(image2)['doy']) ,int(wgt.parseL8filename(image2)['year'])
# 	path,row = int(wgt.parseL8filename(image1)['path']) ,int(wgt.parseL8filename(image1)['row'])			


	dt_1=dt.datetime.fromordinal(dt.date(year1-1,12,31).toordinal()+doy1)
	dt_2=dt.datetime.fromordinal(dt.date(year2-1,12,31).toordinal()+doy2)
	timedelta=dt_2-dt_1
	delts_days=timedelta.days
	mid_dt=dt_1+(timedelta/2)
	start_dec_year=toYearFraction(dt_1.date())
	stop_dec_year=toYearFraction(dt_2.date())
	mid_dec_year=toYearFraction(mid_dt.date())
# 	start_dec_year=wgt.toYearFraction(dt_1.date())
# 	stop_dec_year=wgt.toYearFraction(dt_2.date())
# 	mid_dec_year=wgt.toYearFraction(mid_dt.date())	

	profile_dict={ 'delts_days':delts_days,'image1fn':image1,'image2fn':image2,\
					'series':series_tag, 'start_dec_year':start_dec_year, 'stop_dec_year':stop_dec_year, 'mid_dec_year':mid_dec_year,\
					'path':path, 'row':row, 'infile':filename, 'indir':filepath }
					
	return profile_dict


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




pathRowFolder = '/Users/anderson/Desktop/ARMSTRONG/wrangells/imagery/L8/'
pathRowList = glob.glob(pathRowFolder + '6*/')

allImagesDict = {}

for pr in pathRowList:
	imageList = glob.glob(pr + '/L*/')
	pathRowNum = pr.rsplit('/')[-2]

	dateList = []
	doyList = []
	imageDictList = []
	prDict = {}
	prImagesDict = {}
	
	for image in imageList:
		b1 = glob.glob(image + "*_B1.TIF")[0]
		imageFn = os.path.split(b1)[1]
		imageDict = parseL8filename(imageFn)
		imageDate = doyToDate(int(imageDict['year']),int(imageDict['doy']))
		imageDict['date']=imageDate
		
		dateList.append(imageDate)
		doyList.append(int(imageDict['doy']))
		#imageDictList.append(imageDict)
		
		yearDoy = imageDict['year'] + '-' + imageDict['doy']
		
		prImagesDict[yearDoy] = imageDict
		
	prDict = {'imageDicts':prImagesDict,'doys':doyList,'dates':dateList}
	
	allImagesDict[pathRowNum] = prDict	

### TIME SERIES PLOT ###
#prIter = sorted(allImagesDict).iterkeys()
sortedKeyList =  sorted(allImagesDict)
nums = np.arange(0,len(allImagesDict))
labelList = ['']
dateTicks = [ dt.date(2013,5,1),dt.date(2013,10,1),dt.date(2014,5,1),dt.date(2014,10,1),dt.date(2015,5,1),dt.date(2015,10,1),dt.date(2016,5,1),dt.date(2016,10,1) ]
fig, ax = plt.subplots(figsize=(22,6))
iterNum = 0
colorList = ['skyblue','indianred','lightseagreen','lawngreen','crimson','darkorange','cornflowerblue','darkviolet','darksalmon','khaki']
for i in range(0,len(allImagesDict)):
	prNow = sortedKeyList[i]
	path = int(prNow[0:2])
	row = int(prNow[3:5])
	if row == 18:
		continue 
	
	prLabel = 'Path: %03d, Row: %03d'%(path,row)
	labelList.append(prLabel)
	dates = allImagesDict[prNow]['dates'] 	
	X = pd.to_datetime(dates)
	
	ax.scatter(X, [iterNum]*len(X), c=colorList[iterNum],
			   marker='s', s=100,zorder=1)

	iterNum += 1
				   

fig.autofmt_xdate()

ax.fill_between((dateTicks[0],dateTicks[1]),-1,7,color='lemonchiffon',alpha=0.6,zorder=-1)
ax.fill_between((dateTicks[2],dateTicks[3]),-1,7,color='lemonchiffon',alpha=0.6,zorder=-1)
ax.fill_between((dateTicks[4],dateTicks[5]),-1,7,color='lemonchiffon',alpha=0.6,zorder=-1)
ax.fill_between((dateTicks[6],dateTicks[7]),-1,7,color='lemonchiffon',alpha=0.6,zorder=-1)
#ax.set_yticklabels('')
# everything after this is turning off stuff that's plotted by default
dateFmt = mdates.DateFormatter('%b %Y')
#ax.yaxis.set_visible(False)
ax.tick_params(axis='y', which='both',length=0)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_yticklabels(labelList,fontsize=18)
ax.set_xticks(dateTicks)
ax.set_xticklabels(dateTicks,fontsize=18)
ax.xaxis.set_major_formatter(dateFmt)
ax.xaxis.set_ticks_position('bottom')
#ax.get_yaxis().set_ticklabels([])
#day = pd.to_timedelta("1", unit='D')
#plt.xlim(X[0] - day, X[-1] + day)
plt.savefig('imageFrequency.pdf')
plt.show()	
plt.close()


### DOY PLOT ###
doyAll = np.nan
for pr in allImagesDict:
	#doyAllList.append(np.array(allImagesDict[pr]['doys']))
	doyAll = np.append(doyAll,np.array(allImagesDict[pr]['doys']))
	
hist,edges = np.histogram(doyAll,bins=20,range=(0,366))
plt.plot(edges[1:],hist,lw=2,c='k')
plt.xlabel('DOY',fontsize=16)
plt.ylabel('Number of images',fontsize=16)
plt.title('All scenes - image timing',fontsize=20)
plt.savefig('doy_histogram_allScenes.pdf')
plt.show()
plt.close()

	
### CORRELATION WINDOWS ###

log10folder = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/'

corrList = glob.glob(log10folder + '*_log10.tif')
midDates = []
dtDays = []
paths = []

for corrFn in corrList:
	corrDict = getCorrelationInfoFromFn(log10folder + corrFn)
	midDates.append(corrDict['mid_dec_year'])
	dtDays.append(corrDict['delts_days'])
	paths.append(int(corrDict['path']))


pathColors = (np.array(paths) - 62.)/(67.-62.)
randY = np.random.uniform(0,1000,len(midDates))
sumRand = np.cumsum(randY)
fig,ax=plt.subplots(figsize=(12,6))
ax.errorbar(midDates,sumRand,xerr=np.array(dtDays)/(365.25*2),fmt='o',ecolor='gray', lw=2, capsize=5, capthick=2)
sc=ax.scatter(midDates,sumRand,c=paths,s=100,zorder=100,cmap='Paired')
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticks([])
ax.spines['top'].set_visible(False)
yLims = ax.get_ylim()
ax.tick_params(axis='y', which='top',length=0)
ax.fill_between((2013.33,2013.74),yLims[0],yLims[1],color='lemonchiffon')
ax.fill_between((2014.33,2014.74),yLims[0],yLims[1],color='lemonchiffon')
ax.fill_between((2015.33,2015.74),yLims[0],yLims[1],color='lemonchiffon')
ax.fill_between((2016.33,2016.74),yLims[0],yLims[1],color='lemonchiffon')
ax.set_xlim((2013.25,2016.75))
plt.xlabel('Date',fontsize=18)
ax.ticklabel_format(axis='x',style='plain',useOffset=False,labelsize=20)
#clb = plt.colorbar(sc)
#fig.autofmt_xdate()
fig.savefig('correlationTimespans.pdf')	
fig.clf()

#plt.show()
#plt.close()

#fig.show()
#fig.clf()

	
	
	
	
	
	
	
	
	