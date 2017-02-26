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
from matplotlib import ticker
import wrangellsGeospatialTools as wgt

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

def decYearToDoy(decYear):
	'''
	Turn a decimal year date into DOY
	Not exactly correct, doesn't treat leap years, but good enough
	'''
	
	doy = (decYear-np.floor(decYear)) * 365
	
	return doy


#pathRowFolder = '/Users/anderson/Desktop/ARMSTRONG/wrangells/imagery/L8/'
pathRowFolder = '/Volumes/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/'
imageList = glob.glob(pathRowFolder + 'LC8*')

allImagesDict = {}
midDateList = []
midDoyList = []
dtList = []
pathList = []
rowList = []
pathRowList = []
iterList = []
list62 = []
list63 = []
list64 = []
list65 = []

iterNum=0
for image in imageList:

	corrDict = getCorrelationInfoFromFn(image)
	
	if corrDict['path'] <= 65 and corrDict['row'] <= 17:
		midDateList.append(corrDict['mid_dec_year'])
		dtList.append(corrDict['delts_days'])	
		midDoyList.append(decYearToDoy(corrDict['mid_dec_year']))
		pathList.append(corrDict['path'])
		rowList.append(corrDict['row'])	
		pathRowKey = str(corrDict['path']).zfill(3) + str(corrDict['row']).zfill(3)
		pathRowList.append(pathRowKey)
		
		#print pathRowKey
		
		if pathRowKey == '062017':
			list62.append(corrDict)
		elif pathRowKey == '063017':
			list63.append(corrDict)
		elif pathRowKey == '064017':
			list64.append(corrDict)
		elif pathRowKey == '065017':
			list65.append(corrDict)
								
		iterList.append(iterNum)
		iterNum+=1
	
allImagesDict['062017'] = list62
allImagesDict['063017'] = list63
allImagesDict['064017'] = list64
allImagesDict['065017'] = list65

midDoy62 = []
midDoy63 = []
midDoy64 = []
midDoy65 = []

for i in range(0,len(allImagesDict['062017'])):
	midDoy62.append(decYearToDoy(allImagesDict['062017'][i]['mid_dec_year']))
for i in range(0,len(allImagesDict['063017'])):
	midDoy63.append(decYearToDoy(allImagesDict['063017'][i]['mid_dec_year']))
for i in range(0,len(allImagesDict['064017'])):
	midDoy64.append(decYearToDoy(allImagesDict['064017'][i]['mid_dec_year']))
for i in range(0,len(allImagesDict['065017'])):
	midDoy65.append(decYearToDoy(allImagesDict['065017'][i]['mid_dec_year']))

lateInd62 = np.array(midDoy62) > 300
lateInd63 = np.array(midDoy63) > 300
lateInd64 = np.array(midDoy64) > 300
lateInd65 = np.array(midDoy65) > 300

midDoyT62 = np.copy(np.array(midDoy62))
midDoyT62[lateInd62] = midDoyT62[lateInd62] - 365	

midDoyT63 = np.copy(np.array(midDoy63))
midDoyT63[lateInd63] = midDoyT63[lateInd63] - 365

midDoyT64 = np.copy(np.array(midDoy64))
midDoyT64[lateInd64] = midDoyT64[lateInd64] - 365

midDoyT65 = np.copy(np.array(midDoy65))
midDoyT65[lateInd65] = midDoyT65[lateInd65] - 365

hist62,edges = np.histogram(midDoyT62,bins=np.round(365./10.),range=(-100,300))
hist63,edges = np.histogram(midDoyT63,bins=np.round(365./10.),range=(-100,300))
hist64,edges = np.histogram(midDoyT64,bins=np.round(365./10.),range=(-100,300))
hist65,edges = np.histogram(midDoyT65,bins=np.round(365./10.),range=(-100,300))

lateInd = np.array(midDoyList) > 300
midDoyArr = np.array(midDoyList)
midDoyTransform = np.copy(midDoyArr)
midDoyTransform[lateInd] = midDoyArr[lateInd] - 365	
randNumList = np.random.uniform(0,2,size=(len(midDoyList),1))

hist,edges=np.histogram(midDoyTransform,bins=np.round(365./10.),range=(-100,300))

colorList = []
for i in range(0,len(pathList)):
	if pathList[i] == 62:
		colorList.append('r')
	elif pathList[i] == 63:
		colorList.append('g')
	elif pathList[i] == 64:
		colorList.append('b')
	elif pathList[i] == 65:
		colorList.append('m')		

fig, ax1 = plt.subplots()
ax1.fill_between((120,273),(-1,-1),(40,40),facecolor='yellow',alpha=0.3,lw=1)
ax1.plot(edges[1:],hist,lw=2,c='k')
ax1.fill_between(edges[1:],0,hist62,lw=1,facecolor='r',alpha=0.5)
ax1.fill_between(edges[1:],hist62,hist62+hist63,lw=1,facecolor='g',alpha=0.5)
ax1.fill_between(edges[1:],hist62+hist63,hist62+hist63+hist64,lw=1,facecolor='b',alpha=0.5)
ax1.fill_between(edges[1:],hist62+hist63+hist64,hist62+hist63+hist64+hist65,lw=1,facecolor='m',alpha=0.5)
ax2 = ax1.twinx()
ax2.errorbar(midDoyTransform,randNumList,xerr=np.array(dtList)/2,fmt='o',ecolor='lightgray', lw=2, capsize=5, capthick=2)
ax2.scatter(midDoyTransform,randNumList,c=colorList,s=100,zorder=100,cmap='Set3')
ax2.set_ylim((-7,2.5))
ax2.set_yticks([])
ax1.set_ylim((-1,40))
ax1.set_yticks(np.arange(0,30,5))
ax2.set_xlim((-88,265))
ax1.set_xlim((-88,265))
ax1.set_ylabel('Number of observations',fontsize=18)
ax1.yaxis.set_label_coords(-.05,.35)
ax1.set_xlabel('Day of year',fontsize=18)

p62 = plt.scatter(-10,-10,c='r',s=100,zorder=100,cmap='Set3')
p63 = plt.scatter(-10,-10,c='g',s=100,zorder=100,cmap='Set3')
p64 = plt.scatter(-10,-10,c='b',s=100,zorder=100,cmap='Set3')
p65 = plt.scatter(-10,-10,c='m',s=100,zorder=100,cmap='Set3')

plt.legend((p62,p63,p64,p65),('062','063','064','065'),title="Path number",scatterpoints=1,loc='center left',frameon=False)
plt.savefig('doyCoverage_justSelect_24jan2017.pdf',dpi=500)
plt.show()
plt.close()


	
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

width = 14
dtCount = ( dtDays.count(16*1), dtDays.count(16*2), dtDays.count(16*3), dtDays.count(16*4), dtDays.count(16*5), dtDays.count(16*6), dtDays.count(16*7), dtDays.count(16*8), dtDays.count(16*9), dtDays.count(16*10)) 
dayLabel = np.arange(16,176,16)
plt.bar(dayLabel,dtCount,width=width,tick_label=dayLabel,color='cornflowerblue',edgecolor='darkblue',align='center')
plt.ylim((0,20))
plt.xlim((6,170))
plt.ylabel('Number of correlations [-]',fontsize=18)
plt.xlabel('Correlation timespan [days]',fontsize=18)
plt.savefig('correlationTimespanHistogram.pdf')
plt.show()
plt.close()

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



### TIMING PLOTS FOR INDIVIDUAL TRANSECTS

timingFn = '/Users/wiar9509/Google Drive/wrangells/json/timingDict_fullyAutomated_27jan2017.json'
timingData = wgt.readJsonFile(timingFn)	

keyIter = timingData.iterkeys()
	
for i in range(0,len(timingData)):	
	transNow = keyIter.next()
	timingNow = timingData[transNow]
	summerNow = timingNow['summer']
	winterNow = timingNow['winter']	

	dtList = np.append(summerNow['dDays'],winterNow['dDays'])
	midDoyArr = np.append(summerNow['doys'],winterNow['doys'])
	midDates = np.append(summerNow['midDates'],winterNow['midDates'])
	lateInd = np.array(midDoyArr) > 300
	midDoyTransform = np.copy(midDoyArr)
	midDoyTransform[lateInd] = midDoyArr[lateInd] - 365	
	randNumList = np.random.uniform(0,2,size=(len(midDates),1))

	fig, (ax1,ax2) = plt.subplots(2,1)
	ax1.fill_between(2013+np.array((120,273))/365.,(-1,-1),(3,3),facecolor='yellow',alpha=0.3,lw=1)	
	ax1.fill_between(2014+np.array((120,273))/365.,(-1,-1),(3,3),facecolor='yellow',alpha=0.3,lw=1)		
	ax1.fill_between(2015+np.array((120,273))/365.,(-1,-1),(3,3),facecolor='yellow',alpha=0.3,lw=1)		
	ax1.fill_between(2016+np.array((120,273))/365.,(-1,-1),(3,3),facecolor='yellow',alpha=0.3,lw=1)		
	ax1.errorbar(midDates,randNumList,xerr=np.array(dtList)/(365*2.),fmt='o',ecolor='lightgray', lw=2, capsize=5, capthick=2)
	ax1.scatter(midDates,randNumList,c='dimgray',s=100,zorder=100,cmap='Set3')	
	ax2.fill_between((120,273),(-1,-1),(3,3),facecolor='yellow',alpha=0.3,lw=1)	
	ax2.errorbar(midDoyTransform,randNumList,xerr=np.array(dtList)/2.,fmt='o',ecolor='lightgray', lw=2, capsize=5, capthick=2)
	ax2.scatter(midDoyTransform,randNumList,c='dimgray',s=100,zorder=100,cmap='Set3')



	x_formatter = ticker.ScalarFormatter(useOffset=False)
	ax1.xaxis.set_major_formatter(x_formatter)
	ax1.set_yticks([])
	ax2.set_yticks([])
	ax1.set_title(transNow + " observation timing",fontsize=20)
	ax1.set_ylabel('Year',fontsize=18)	
	ax2.set_ylabel('DOY',fontsize=18)	
	ax1.set_ylim((-.5,2.5))	
	ax2.set_ylim((-.5,2.5))
	ax2.set_xlim((-100,300))
	ax1.set_xlim((2013,2017))
	
	plt.savefig(transNow + '_observationTiming.pdf',dpi=500)
	#plt.show()
	plt.close()	
	
	
	