'''
Testing uniformity and magnitude of speedup
07feb2017
wha
'''

import numpy as np
import json
import matplotlib.pyplot as plt
import wrangellsGeospatialTools as wgt
import matplotlib.gridspec as gridspec
from scipy.stats import linregress

def calculateRMSE(obs,model):
	'''
	Calculate root mean squared error
	obs = observations
	model = best fit/model
	rmse = output root mean square error
	
	Could also do:
	def rmse(predictions, targets):
 	   return np.sqrt(((predictions - targets) ** 2).mean())
 	   
	'''
	resid = obs - model
	residSquare = resid**2
	residSumSquare = np.nansum(residSquare)
	meanSumSquare = residSumSquare / len(obs)
	rmse = np.sqrt(meanSumSquare)
	
	return rmse
	
		
seasonalVelocityFn = 'seasonalVelocityProfiles_2iterTempFilter_fullyAutomated_24jan2017.json'
speedupFn = 'speedupDict_2iterTempFilter_fullyAutomated_28nov2016.json'
elevationFn = '/Users/wiar9509/Google Drive/wrangells/json/new_centerlines_withChugach_epsg102006_elevationProfiles.json'
ranksumFn = '/Users/wiar9509/Google Drive/wrangells/json/ranksumSpeedupDict_30nov2016.json'
snowlineFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/excel/snowlineData.csv'
# Open json files
speedupDict = wgt.readJsonFile(speedupFn)
velocityDict = wgt.readJsonFile(seasonalVelocityFn)
elevationDict = wgt.readJsonFile(elevationFn)
ranksumDict = wgt.readJsonFile(ranksumFn)

# Get snowline data
snowlineData=np.genfromtxt(snowlineFn,delimiter=',',skip_header=1,dtype=None)
quantifySpeedupDict = {}

# Iterate over profiles
keyIter = speedupDict.iterkeys()
for i in range(0,len(speedupDict)):
	transNow = keyIter.next()
	
	speedupProfNow = speedupDict[transNow]['tempFilter']
	velocityProfNow = velocityDict[transNow]	
	elevNow = elevationDict[transNow]
	dist = speedupDict[transNow]['termDist']
			
	# Make distance relative to snowline
	# Get snowline data for this transect (v. clunky..)
	for j in range(0,len(snowlineData)):
		snowlineTrans = snowlineData[j][1]
		try:
			if snowlineTrans == transNow:
				snowlineDistFromTop = snowlineData[j][-1] # I determined this as 'distance from headwall' due to how the shapefile wasstructured
				snowlineHeight = snowlineData[j][2]			
				snowlineDist = np.nanmax(dist) - snowlineDistFromTop # this now in distance from terminus, as all other measurements
		except: # if didn't calculate snowline
			snowlineDist = 0 # these just no-data, should keep snowDist == dist
			snowlineHeight = 0
	
	# Make a distance array relative to snowline postion
	belowSnowInd = np.where(np.array(dist)<=snowlineDist)
	snowDist = np.array(dist) - snowlineDist
				
	
	# Pull median seasonal velocity and speedup profiles
	su50 = speedupProfNow['summer50']
	wi50 = speedupProfNow['winter50']
	sp25 = speedupProfNow['minDiff'] # speedup profiles
	sp50 = speedupProfNow['medianDiff']
	sp75 = speedupProfNow['maxDiff'] # speedup profiles
	
	
	# no ranksum file exists if no speedup
	try:
		ranksumNow = ranksumDict[transNow]
		speedupLength = ranksumNow['speedupLength']
		#dist = speedupDict[transNow]['termDist']
		dist = snowDist # doing relative to snowline
		speedupLength = speedupLength - snowlineDist # converting speedup length to snow line reference
		
		print "Processing " + transNow
		 
		# Ranksum
		ugInd = np.where(np.array(dist)>=speedupLength)[0]
		dgInd = np.where(np.array(dist)<speedupLength)[0]
	
		distDown = np.array(dist)[dgInd]	
		distUp = np.array(dist)[ugInd]
		spDown = np.array(sp50)[dgInd]
		spUp = np.array(sp50)[ugInd]
			
		# Find nan's and cut
		nanIndUp = np.where(np.isnan(spUp))[0]
		nanIndDown = np.where(np.isnan(spDown))[0]		
		distNoNanUp = np.delete(distUp,nanIndUp)
		spNoNanUp = np.delete(spUp,nanIndUp)				
		distNoNanDown = np.delete(distDown,nanIndDown)		
		spNoNanDown = np.delete(spDown,nanIndDown)				
				
		# Fit a first order polynomial to speedup
		polyUp = np.polyfit(distNoNanUp,spNoNanUp,deg=1) 
		polyDown = np.polyfit(distNoNanDown,spNoNanDown,deg=1) 	

		# Residuals
		# best fit speedup values
		speedMagDown = polyDown[1] # offset from 0
		speedSlopeDown = polyDown[0] # speedup slope (m/d/m)
		speedSlopeUp = polyUp[0] # speedup slope (m/d/m)
		speedMagUp = polyUp[1] + speedSlopeUp*np.nanmin(distUp)  # offset from 0
			
		# Best fit first order polynomial speedup values				
		bfUp = polyUp[0]*distUp+polyUp[1] 
		bfDown = polyDown[0]*distDown+polyDown[1] 
		# Residuals
		redsidUp = spUp - bfUp
		residDown = spDown - bfDown		
		# Mean of best fit line up/down glacier
		speedupMeanDown = np.nanmean(bfDown)
		speedupMeanUp = np.nanmean(bfUp)
		
				
		# Linear regression in scipy.stats (same thing, but gives r^2)
		lrUp = linregress(distUp,spUp)
		lrDown = linregress(distDown,spDown)
		# r-squared
		r2up = lrUp[2]
		r2down = lrDown[2]
		
		# Root mean square error
		rmseUp = calculateRMSE(spUp,bfUp)
		rmseDown = calculateRMSE(spDown,bfDown)				
		
		
		if transNow in ('barnWB','donjWB','donjEB','kluvCL','hawkCL','wForCL'):
			print 'Special case for: ' + transNow
			
			if transNow == 'barnWB':
				speedupStart = 18889			
				speedupEnd = 30035
			elif transNow == 'donjWB':
				speedupStart = 6686
				speedupEnd = 19935
			elif transNow == 'donjEB':
				speedupStart = 4774
				#speedupEnd = 23766
				speedupEnd = 19000 # doing this at later date than rest (16feb2017)
			elif transNow == 'kluvCL':
				speedupStart = 3644
				speedupEnd = 10301
			elif transNow == 'hawkCL':
				speedupStart = 9204
				speedupEnd = 25168			
			elif transNow == 'wForCL':
				speedupStart = 2131
				speedupEnd = 8438
			
			speedupStart = speedupStart - snowlineDist # converting
			speedupEnd = speedupEnd - snowlineDist # converting
						
			# Ranksum
			ugInd2 = np.where(np.array(dist)>=speedupEnd)[0]
			dgInd2 = np.where(np.logical_and(np.array(dist)>speedupStart,np.array(dist)<speedupEnd))[0]

			distDown2 = np.array(dist)[dgInd2]	
			distUp2 = np.array(dist)[ugInd2]
			spDown2 = np.array(sp50)[dgInd2]
			spUp2 = np.array(sp50)[ugInd2]
			
			# Find nan's and cut
			nanIndUp2 = np.where(np.isnan(spUp2))[0]
			nanIndDown2 = np.where(np.isnan(spDown2))[0]		
			distNoNanUp2 = np.delete(distUp2,nanIndUp2)
			spNoNanUp2 = np.delete(spUp2,nanIndUp2)				
			distNoNanDown2 = np.delete(distDown2,nanIndDown2)		
			spNoNanDown2 = np.delete(spDown2,nanIndDown2)				
				
			# Fit a first order polynomial to speedup
			polyUp2 = np.polyfit(distNoNanUp2,spNoNanUp2,deg=1) 
			polyDown2 = np.polyfit(distNoNanDown2,spNoNanDown2,deg=1) 	

			# Residuals
			# best fit speedup values
			speedSlopeDown2 = polyDown2[0] # speedup slope (m/d/m)
			speedSlopeUp2 = polyUp2[0] # speedup slope (m/d/m)
			speedMagDown2 = polyDown2[1] + speedSlopeDown2*np.nanmin(distDown2) # offset from 0
			speedMagUp2 = polyUp2[1] + speedSlopeUp2*np.nanmin(distUp2)  # offset from 0
			
			# Best fit first order polynomial speedup values				
			bfUp2 = polyUp2[0]*distUp2+polyUp2[1] 
			bfDown2 = polyDown2[0]*distDown2+polyDown2[1] 
			# Residuals
			redsidUp2 = spUp2 - bfUp2
			residDown2 = spDown2 - bfDown2		
			# Mean of best fit line up/downglacier
			speedupMeanDown2 = np.nanmean(bfDown2)
			speedupMeanUp2 = np.nanmean(bfUp2)

			# Linear regression in scipy.stats (same thing, but gives r^2)
			lrUp2 = linregress(distUp2,spUp2)
			lrDown2 = linregress(distDown2,spDown2)
			# r-squared
			r2up2 = lrUp2[2]
			r2down2 = lrDown2[2]
		
			# Root mean square error
			rmseUp2 = calculateRMSE(spUp2,bfUp2)
			rmseDown2 = calculateRMSE(spDown2,bfDown2)				

			# Write dict for saving out
			downGlacierStats2 = {'speedupMag':float(speedMagDown2),'meanSpeedup':float(speedupMeanDown2),'speedupSlope':float(speedSlopeDown2),'bestFitSpeedup':bfDown2.tolist(),'rSquared':float(r2down2),'speedup':spDown2.tolist(),'rmse':float(rmseDown2),'dist':distDown2.tolist()}
			upGlacierStats2 = {'speedupMag':float(speedMagUp2),'meanSpeedup':float(speedupMeanDown2),'speedupSlope':float(speedSlopeUp2),'bestFitSpeedup':bfUp2.tolist(),'rSquared':float(r2up2),'speedup':spUp2.tolist(),'rmse':float(rmseUp2),'dist':distUp2.tolist()}
			glacierStats2 = {'downGlacier':downGlacierStats2,'upGlacier':upGlacierStats2}				
															
		# Write dict for saving out
		downGlacierStats = {'speedupMag':float(speedMagDown),'meanSpeedup':float(speedupMeanDown),'speedupSlope':float(speedSlopeDown),'bestFitSpeedup':bfDown.tolist(),'rSquared':float(r2down),'speedup':spDown.tolist(),'rmse':float(rmseDown),'dist':distDown.tolist()}
		upGlacierStats = {'speedupMag':float(speedMagUp),'meanSpeedup':float(speedupMeanUp),'speedupSlope':float(speedSlopeUp),'bestFitSpeedup':bfUp.tolist(),'rSquared':float(r2up),'speedup':spUp.tolist(),'rmse':float(rmseUp),'dist':distUp.tolist()}
		glacierStats = {'downGlacier':downGlacierStats,'upGlacier':upGlacierStats,'snowlineHeight':int(snowlineHeight),'snowlineTermDist':float(snowlineDist),'distToSnow':snowDist.tolist()}
		
		if transNow in ('barnWB','donjWB','donjEB','kluvCL','hawkCL','wForCL'):
			quantifySpeedupDict[transNow]={'automated':glacierStats,'handpicked':glacierStats2}
		else:
			quantifySpeedupDict[transNow]={'automated':glacierStats}
		
		# plot
		if 1:
			if transNow in ('barnWB','donjWB','donjEB','kluvCL','hawkCL','wForCL'):	
				plt.subplot(2,1,1)
				plt.fill_between(np.array(distUp)/1e3,np.array(sp25)[ugInd],np.array(sp75)[ugInd],facecolor='cornflowerblue',alpha=0.5)
				plt.fill_between(np.array(distDown)/1e3,np.array(sp25)[dgInd],np.array(sp75)[dgInd],facecolor='salmon',alpha=0.5)		
				plt.plot(np.array(dist)/1e3,sp50,c='k')
				plt.plot((np.array(dist)[ugInd])/1e3,np.array(sp50)[ugInd],c='darkblue',lw=2)
				plt.plot((np.array(dist)[dgInd])/1e3,np.array(sp50)[dgInd],c='darkred',lw=2)
				plt.plot(np.array((speedupLength,speedupLength))/1e3,(-0.5,0.5),c='k',ls=':',lw=1)
				plt.plot(np.array((speedupStart,speedupStart))/1e3,(-0.5,0.5),c='k',ls='--',lw=1)				
				plt.plot(np.array((speedupEnd,speedupEnd))/1e3,(-0.5,0.5),c='k',ls='--',lw=1)								
				plt.plot((np.nanmin(dist),np.nanmax(dist)/1e3),(0,0),c='k',ls='--',lw=1)
				plt.plot(np.array(distUp)/1e3,polyUp[0]*distUp+polyUp[1],c='gray',lw=1)
				plt.plot(np.array(distDown)/1e3,polyDown[0]*distDown+polyDown[1],c='gray',lw=1)
				plt.plot(np.array(distUp2)/1e3,bfUp2,c='dimgray',lw=2)
				plt.plot(np.array(distDown2)/1e3,bfDown2,c='dimgray',lw=2)				
				
				plt.title(transNow,fontsize=20)
				plt.ylabel('Speedup [m/d]',fontsize=16)
				plt.xlabel('Distance from snowline [km]',fontsize=16)	
				plt.ylim((-0.5,0.5))
				if 0:
					plt.text(12,-0.15,'DG mag: ' + "{0:.3f}".format(speedupMeanDown2))
					plt.text(12,-0.25,'DG slope: ' + "{0:.1E}".format(speedSlopeDown2))		
					plt.text(12,-0.35,'UG mag: ' + "{0:.3f}".format(speedupMeanUp2))
					plt.text(12,-0.45,'UG slope: ' + "{0:.1E}".format(speedSlopeUp2))				
				plt.xlim((np.nanmin(dist)/1e3,np.nanmax(dist)/1e3))
				plt.gca().invert_xaxis()	
				
				plt.subplot(2,1,2)
				#plt.plot(distDown,rDown,c='r')
				#plt.plot(distUp,rUp,c='b')		
				plt.scatter(spUp2,bfUp2,facecolor='b')
				plt.scatter(spDown2,bfDown2,facecolor='r')
				plt.plot((-0.2,0.2),(-0.2,0.2),c='k',ls='--')
				plt.text(-.19,.05,'UG r2: ' + "{0:.3f}".format(r2up2) + '\n' +
				'UG rmse: ' + "{0:.3f}".format(rmseUp2) + '\nDG r2: ' + "{0:.3f}".format(r2down2) + '\n' + 
				'DG rmse: ' + "{0:.3f}".format(rmseDown2))
			else:
				plt.subplot(2,1,1)
				plt.fill_between(np.array(distUp)/1e3,np.array(sp25)[ugInd],np.array(sp75)[ugInd],facecolor='cornflowerblue',alpha=0.5)
				plt.fill_between(np.array(distDown)/1e3,np.array(sp25)[dgInd],np.array(sp75)[dgInd],facecolor='salmon',alpha=0.5)		
				plt.plot(np.array(dist)/1e3,sp50,c='k')
				plt.plot((np.array(dist)[ugInd])/1e3,np.array(sp50)[ugInd],c='darkblue',lw=2)
				plt.plot((np.array(dist)[dgInd])/1e3,np.array(sp50)[dgInd],c='darkred',lw=2)
				plt.plot(np.array((speedupLength,speedupLength))/1e3,(-0.5,0.5),c='k',ls='--',lw=1)
				plt.plot((np.nanmin(dist),np.nanmax(dist)/1e3),(0,0),c='k',ls='--',lw=1)
				plt.plot(np.array(distUp)/1e3,polyUp[0]*distUp+polyUp[1],c='dimgray',lw=2)
				plt.plot(np.array(distDown)/1e3,polyDown[0]*distDown+polyDown[1],c='dimgray',lw=2)
				plt.title(transNow,fontsize=20)
				plt.ylabel('Speedup [m/d]',fontsize=16)
				plt.xlabel('Distance from snowline [km]',fontsize=16)	
				plt.ylim((-0.5,0.5))
				if 0:
					plt.text(12,-0.15,'DG mag: ' + "{0:.3f}".format(speedupMeanDown))
					plt.text(12,-0.25,'DG slope: ' + "{0:.1E}".format(speedSlopeDown))		
					plt.text(12,-0.35,'UG mag: ' + "{0:.3f}".format(speedupMeanUp))
					plt.text(12,-0.45,'UG slope: ' + "{0:.1E}".format(speedSlopeDown))				
				plt.xlim((np.nanmin(dist)/1e3,np.nanmax(dist)/1e3))
				plt.gca().invert_xaxis()	
				
				plt.subplot(2,1,2)
				#plt.plot(distDown,rDown,c='r')
				#plt.plot(distUp,rUp,c='b')		
				plt.scatter(spUp,bfUp,facecolor='b')
				plt.scatter(spDown,bfDown,facecolor='r')
				plt.plot((-0.2,0.2),(-0.2,0.2),c='k',ls='--')
				plt.text(-.19,.05,'UG r2: ' + "{0:.3f}".format(r2up) + '\n' +
				'UG rmse: ' + "{0:.3f}".format(rmseUp) + '\nDG r2: ' + "{0:.3f}".format(r2down) + '\n' + 
				'DG rmse: ' + "{0:.3f}".format(rmseDown))
					
	# 		plt.text(-.19,-.13,'UG rmse: ' + "{0:.3f}".format(rmseUp))				
	# 		plt.text(-.19,-.09,'DG r2: ' + "{0:.3f}".format(r2down))
	# 		plt.text(-.19,-0.03,'DG rmse: ' + "{0:.3f}".format(rmseDown))						
			plt.xlim((-.2,.2))
			plt.ylim((-.2,.2))	
			plt.xlabel('Observed speedup [m/d]',fontsize=16)
			plt.ylabel('Best linear \n fit speedup [m/d]',fontsize=16)				
			plt.savefig(transNow + '_quantifySpeedupPlot_relToSnowline.pdf')
			#plt.show()
			plt.close()				
		
	except KeyError:
		print "Key error on " + transNow
		ranksumExists = 0
		
wgt.writeJsonFile(quantifySpeedupDict,'quantifySpeedupDict_withSnowline_handTreated_07feb2017.json')


# Look at all glaciers together
#data = wgt.readJsonFile('/Users/wiar9509/Google Drive/wrangells/json/quantifySpeedupDict_07feb2017.json')

qSpeed = wgt.readJsonFile('/Users/wiar9509/Google Drive/wrangells/json/quantifySpeedupDict_withSnowline_handTreated_07feb2017.json')
# initialize
transNameList = []
speedMagListUp = []
speedMagListDown = []
speedSlopeListUp = []
speedSlopeListDown = []

color_idx = np.linspace(0, 1, len(qSpeed))
cols = plt.cm.Paired(color_idx)

keyIter = qSpeed.iterkeys()
	
for i in range(0,len(qSpeed)):
	transNow = keyIter.next()
	#print transNow, i
	i = 4
	transNameList.append(str(transNow))
	if transNow in ('barnWB','donjWB','donjEB','kluvCL','hawkCL','wForCL'):		
		speedMagListDown.append(float(qSpeed[transNow]['handpicked']['downGlacier']['meanSpeedup']))
		speedSlopeListDown.append(float(qSpeed[transNow]['handpicked']['downGlacier']['speedupSlope']))
		speedMagListUp.append(float(qSpeed[transNow]['handpicked']['upGlacier']['meanSpeedup']))
		speedSlopeListUp.append(float(qSpeed[transNow]['handpicked']['upGlacier']['speedupSlope']))	
		bestFitUp = qSpeed[transNow]['handpicked']['upGlacier']['bestFitSpeedup']
		bestFitDown = qSpeed[transNow]['handpicked']['downGlacier']['bestFitSpeedup']		
		distUp = qSpeed[transNow]['handpicked']['upGlacier']['dist']
		distDown = qSpeed[transNow]['handpicked']['downGlacier']['dist']		
	else:
		speedMagListDown.append(float(qSpeed[transNow]['automated']['downGlacier']['meanSpeedup']))
		speedSlopeListDown.append(float(qSpeed[transNow]['automated']['downGlacier']['speedupSlope']))
		speedMagListUp.append(float(qSpeed[transNow]['automated']['upGlacier']['meanSpeedup']))
		speedSlopeListUp.append(float(qSpeed[transNow]['automated']['upGlacier']['speedupSlope']))		
		bestFitUp = qSpeed[transNow]['automated']['upGlacier']['bestFitSpeedup']
		bestFitDown = qSpeed[transNow]['automated']['downGlacier']['bestFitSpeedup']		
		distUp = qSpeed[transNow]['automated']['upGlacier']['dist']
		distDown = qSpeed[transNow]['automated']['downGlacier']['dist']		
					
#	print transNow, "{0:.3f}".format(speedMagListDown[i]), "{0:.2E}".format(speedSlopeListDown[i])
	
	# plot
	if 1:
		speedupNow = speedupDict[transNow]
		termDist = speedupNow['termDist']
		distToSnow = qSpeed[transNow]['automated']['distToSnow']
		sp50 = speedupNow['tempFilter']['medianDiff']
		sp25 = speedupNow['tempFilter']['minDiff']
		sp75 = speedupNow['tempFilter']['maxDiff']		
		randInd = np.random.randint(0,len(qSpeed))
		#plt.subplot(2,1,1)
		#plt.fill_between(np.array(termDist)/1e3,sp25,sp75,lw=0,facecolor=cols[i],alpha=0.1)
		plt.plot(np.array(termDist)/1e3,sp50,marker='.',markersize=4,lw=0.0,alpha=0.7,c=cols[i])
		plt.plot(np.array(distUp)/1e3,bestFitUp,lw=1,alpha=1,c=cols[i])
		plt.plot(np.array(distDown)/1e3,bestFitDown,lw=2,alpha=1,c=cols[i],label=transNow)		
		#plt.subplot(2,1,2)
		#plt.plot(np.array(distToSnow)/1e3,sp50,marker='.',lw=0.5,alpha=0.9,c=cols[i],label=transNow)
			
# continuing plotting from loop
if 1:
	#plt.subplot(2,1,1)
	plt.ylim((-0.25,0.35))
	plt.xlim((0,80))	
	plt.plot((0,80),(0,0),ls='--',c='k')
	plt.xlabel('Distance from terminus [km]',fontsize=16)
	plt.ylabel('Speedup [m d$^{-1}$]',fontsize=18)
	if 0:
		plt.subplot(2,1,2)
		plt.ylim((-0.25,0.35))
		plt.plot((-60,30),(0,0),ls='--',c='k')	
		plt.xlim((-60,30))
		plt.xlabel('Distance from snowline [km]',fontsize=16)		
		plt.ylabel('Speedup [m d$^{-1}$]',fontsize=18)
	plt.legend(loc=1,numpoints=1,frameon=False,ncol=3,prop={'size':8})
	#plt.savefig('speedupSquiggles_withFitLines_noFill_16feb2017.pdf')
	plt.show()
	plt.close()

dataOut = np.vstack((transNameList,speedMagListDown,speedSlopeListDown,speedMagListUp,speedSlopeListUp)).transpose()		
np.savetxt('test.csv',dataOut,fmt='%s,%s,%s,%s,%s',delimiter=',')		
		