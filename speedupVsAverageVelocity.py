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

# Open json files
speedupDict = wgt.readJsonFile(speedupFn)
velocityDict = wgt.readJsonFile(seasonalVelocityFn)

# Iterate over profiles
keyIter = speedupDict.iterkeys()
for i in range(0,len(speedupDict)):
	transNow = keyIter.next()
	
	speedupProfNow = speedupDict[transNow]['tempFilter']
	velocityProfNow = velocityDict[transNow]	
	
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
	sp50 = speedupProfNow['medianDiff']

	# Calculate average velocity
	uAvg1 = summerWeight1*np.array(su50) + winterWeight1*np.array(wi50) # based on DOYs used
	uAvg2 = summerWeight2*np.array(su50) + winterWeight2*np.array(wi50) # long summer	
	uAvg3 = summerWeight3*np.array(su50) + winterWeight3*np.array(wi50) # short summer	
	
	
	# Plot
	xLims = ((np.max(np.array(dist)/1e3),-0.5))
	#fig, (ax0, ax1) = plt.subplots(2, 1)
	gs = gridspec.GridSpec(2, 1,height_ratios=[4,1])

	ax0 = plt.subplot(gs[0])
	ax1 = plt.subplot(gs[1])
	
	ax0.plot(np.array(dist)/1e3,su50,lw=2,c='darkred',ls='-',label='$u_{summer}$')
	ax0.plot(np.array(dist)/1e3,wi50,lw=2,c='darkblue',ls='-',label='$u_{winter}$')
	ax0.plot(np.array(dist)/1e3,uAvg3,lw=2,c='lightgray',ls='-')
	ax0.plot(np.array(dist)/1e3,uAvg1,lw=2,c='slategray',ls='-',label='$u_{annual}$')	
	ax0.plot(np.array(dist)/1e3,uAvg2,lw=2,c='dimgray',ls='-')		
	ax0.plot(np.array(dist)/1e3,sp50,lw=2,c='crimson',ls=':',label='$\Delta u=u_s-u_w$')			
	ax0.plot(xLims,(0,0),lw=1,ls='--',c='k')

	#ax0.invert_xaxis()

	ax1.plot(np.array(dist)/1e3,np.array(sp50)/np.array(uAvg3),lw=2,c='slategray',ls='-',label='$\Delta u / u_a$')				
	ax1.plot(np.array(dist)/1e3,np.array(sp50)/np.array(su50),lw=2,c='darkred',ls='-',label='$\Delta u / u_s$')
	ax1.plot(np.array(dist)/1e3,np.array(sp50)/np.array(wi50),lw=2,c='darkblue',ls='-',label='$\Delta u / u_w$')				
	ax1.plot(xLims,(0,0),lw=1,ls='--',c='k')
	
	ax0.invert_xaxis()
	ax1.invert_xaxis()
	
	ax0.set_xlim(xLims)
	ax1.set_xlim(xLims)
	ax0.set_xticks([])
	
	ax0.legend(loc=0,ncol=2,frameon=False,prop={'size':12})
	ax1.legend(loc=0,ncol=3,frameon=False,prop={'size':10})	
	ax0.set_ylabel('Velocity [m d$^{-1}$]',fontsize=18)
	ax1.set_ylabel('Normalized \n velocity [-]',fontsize=18)	
	ax0.yaxis.set_label_coords(-.06,.50)	
	ax1.yaxis.set_label_coords(-.06,.50)	
	ax0.set_title(transNow,fontsize=20)
	ax1.set_xlabel('Distance from terminus [km]',fontsize=18)
	
	plt.savefig(transNow + '_normalizedSpeedup.pdf')	
	#plt.show()
	plt.close()
	
	