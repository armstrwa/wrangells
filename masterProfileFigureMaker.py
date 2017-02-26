'''
Master profile figure maker
27jan2017
'''

'''
Plot speedup as a function of annual average velocity
25 jan 2017
'''

import numpy as np
import json
import matplotlib.pyplot as plt
import wrangellsGeospatialTools as wgt
import matplotlib.gridspec as gridspec

seasonalVelocityFn = 'seasonalVelocityProfiles_2iterTempFilter_fullyAutomated_24jan2017.json'
speedupFn = 'speedupDict_2iterTempFilter_fullyAutomated_28nov2016.json'
elevationFn = '/Users/wiar9509/Google Drive/wrangells/json/new_centerlines_withChugach_epsg102006_elevationProfiles.json'
ranksumFn = '/Users/wiar9509/Google Drive/wrangells/json/ranksumSpeedupDict_30nov2016.json'
# Open json files
speedupDict = wgt.readJsonFile(speedupFn)
velocityDict = wgt.readJsonFile(seasonalVelocityFn)
elevationDict = wgt.readJsonFile(elevationFn)
ranksumDict = wgt.readJsonFile(ranksumFn)

# Iterate over profiles
keyIter = speedupDict.iterkeys()
for i in range(0,len(speedupDict)):
	transNow = keyIter.next()
	
	speedupProfNow = speedupDict[transNow]['tempFilter']
	velocityProfNow = velocityDict[transNow]	
	elevNow = elevationDict[transNow]
	
	# no ranksum file exists if no speedup
	try:
		ranksumNow = ranksumDict[transNow]
		speedupLength = ranksumNow['speedupLength']
	except KeyError:
		ranksumExists = 0
		
	
	dist = speedupDict[transNow]['termDist']
	
	# Weighting summer velocity in annual average
	# Proportionate to fraction of year
	summerWeight1 = (250-150)/365. 
	winterWeight1 = 1 - summerWeight1
	
	# Testing uncertainty effects on annual average
	summerWeight2 = 0.40	
	winterWeight2 = 0.60	
	
	summerWeight3 = 0.20	
	winterWeight3 = 0.80	
	
	# Pull median seasonal velocity and speedup profiles
	su50 = speedupProfNow['summer50']
	wi50 = speedupProfNow['winter50']
	sp25 = speedupProfNow['minDiff'] # speedup profiles
	sp50 = speedupProfNow['medianDiff']
	sp75 = speedupProfNow['maxDiff'] # speedup profiles
	
	# Calculate average velocity
	uAvg1 = summerWeight1*np.array(su50) + winterWeight1*np.array(wi50) # based on DOYs used
	uAvg2 = summerWeight2*np.array(su50) + winterWeight2*np.array(wi50) # long summer	
	uAvg3 = summerWeight3*np.array(su50) + winterWeight3*np.array(wi50) # short summer	
	
	
	# Data numbers
	dNumSu = speedupProfNow['summerDataNum']
	dNumWi = speedupProfNow['winterDataNum']
		
	# Elevation profile
	zMean = elevNow['elevMean']
	dzdx = elevNow['slope']	
		
	minZ = np.nanmin(zMean)
	maxZ = np.nanmax(zMean)
	rangeZ = maxZ - minZ

	# Ranksum
	ugInd = np.where(np.array(dist)>=speedupLength)[0]
	dgInd = np.where(np.array(dist)<speedupLength)[0]
	
	distDown = np.array(dist)[dgInd]	
	distUp = np.array(dist)[ugInd]
	spDown = np.array(sp50)[dgInd]
	spUp = np.array(sp50)[ugInd]	
	
	# Fit a first order polynomial to speedup
	polyUp = np.polyfit(distUp,spUp,deg=1) 
	polyDown = np.polyfit(distDown,spDown,deg=1) 	

	plt.plot(dist,sp50,c='r')
	plt.plot((speedupLength,speedupLength),(-0.5,1.5),c='k',ls='--',lw=1)
	plt.plot(distUp,polyUp[0]*distUp+polyUp[1],c='gray',lw=1)
	plt.plot(distDown,polyDown[0]*distDown+polyDown[1],c='gray',lw=1)
	plt.show()
	plt.close()				

	# Plot
	if len(dist) == len(zMean):
		xLims = ((np.max(np.array(dist)/1e3),-0.5))
		#fig, (ax0, ax1) = plt.subplots(2, 1)
		gs = gridspec.GridSpec(4, 1,height_ratios=[1,3,2,1])

		ax00 = plt.subplot(gs[0])
		ax0 = plt.subplot(gs[1])
		ax1 = plt.subplot(gs[2])
		ax11 = plt.subplot(gs[3])

		ax00.plot(np.array(dist)/1e3,dNumSu,lw=2,c='darkred')
		ax00.plot(np.array(dist)/1e3,dNumWi,lw=2,c='darkblue')

		ax0.plot(np.array(dist)/1e3,su50,lw=2,c='darkred',ls='-',label='$u_{summer}$')
		ax0.plot(np.array(dist)/1e3,wi50,lw=2,c='darkblue',ls='-',label='$u_{winter}$')
		ax0.plot(np.array(dist)/1e3,uAvg3,lw=2,c='lightgray',ls='-')
		ax0.plot(np.array(dist)/1e3,uAvg1,lw=2,c='slategray',ls='-',label='$u_{annual}$')	
		ax0.plot(np.array(dist)/1e3,uAvg2,lw=2,c='dimgray',ls='-')		
		ax0.fill_between(np.array(dist)/1e3,sp25,sp75,facecolor='salmon',linestyle=':',alpha=0.5)
		ax0.plot(np.array(dist)/1e3,sp50,lw=2,c='crimson',ls=':',label='$\Delta u=u_s-u_w$')			
		ax0.plot(xLims,(0,0),lw=1,ls='--',c='k')

		if ranksumExists == 1:
			ax0.plot(np.array((speedupLength,speedupLength))/1e3,(-0.5,1.5),c='k',ls='--',lw=1)
			ax0.plot(np.array(distUp)/1e3,polyUp[0]*distUp+polyUp[1],c='gray',lw=1)
			ax0.plot(np.array(distDown)/1e3,polyDown[0]*distDown+polyDown[1],c='gray',lw=1)
			
			
		#ax0.invert_xaxis()

		ax1.plot(np.array(dist)/1e3,np.array(sp50)/np.array(uAvg3),lw=2,c='slategray',ls='-',label='$\Delta u / u_a$')				
		ax1.plot(np.array(dist)/1e3,np.array(sp50)/np.array(su50),lw=2,c='darkred',ls='-',label='$\Delta u / u_s$')
		ax1.plot(np.array(dist)/1e3,np.array(sp50)/np.array(wi50),lw=2,c='darkblue',ls='-',label='$\Delta u / u_w$')				

		ax1.plot(xLims,(0,0),lw=1,ls='--',c='k')

		ax11.plot(np.array(dist)/1e3,zMean,c='k',lw=2)

		ax00.invert_xaxis()	
		ax0.invert_xaxis()
		ax1.invert_xaxis()
		ax11.invert_xaxis()
	
		ax00.grid(True)
		ax0.grid(True)
		ax1.grid(True)
		ax11.grid(True)
								
		ax00.set_xlim(xLims)
		ax0.set_xlim(xLims)
		ax1.set_xlim(xLims)
		ax11.set_xlim(xLims)

		if 1: # standard axes
			ax00.set_yticks(np.arange(0,40,10))
			ax00.set_ylim([0,30])
			ax0.set_yticks(np.arange(-0.5,2.,0.5))
			ax0.set_ylim([-0.5,1.5])
			ax1.set_ylim([-2,2])
			ax1.set_yticks(np.arange(-1,2,2))
			ax11.set_ylim([500,3000])
			ax11.set_yticks(np.arange(500,3500,1000))
			ax00.set_xticks([])
			ax0.set_xticks([])
			ax1.set_xticks([])
		elif 0: # auto axes	
			print 'do something'
			#ax00.set_yticks([np.arange(np.round(minZ),np.round(maxZ),np.round(rangeZ/5))])

		ax00.yaxis.tick_right()
		ax00.yaxis.set_label_position("right")
		ax1.yaxis.tick_right()
		ax1.yaxis.set_label_position("right")	
			
		ax0.legend(loc=0,ncol=2,frameon=False,prop={'size':12})
		ax1.legend(loc=0,ncol=3,frameon=False,prop={'size':10})	
		ax00.set_ylabel('Num. obs. [-]',fontsize=18)
		ax0.set_ylabel('Velocity [m d$^{-1}$]',fontsize=18)
		ax1.set_ylabel('Norm. velocity [-]',fontsize=18)
		ax11.set_ylabel('Elevation [m]',fontsize=18)		
		ax00.yaxis.set_label_coords(1.06,.50)	
		ax1.yaxis.set_label_coords(1.06,.50)	
		ax00.set_title(transNow,fontsize=20)
		ax11.set_xlabel('Distance from terminus [km]',fontsize=18)

		#plt.savefig(transNow + '_standardAxes_countVeloSpeedupElevation.pdf')	
		plt.show()
		plt.close()
	
	