'''
Average profile coverage
24 jan 2017
wha
'''

import glob
import wrangellsGeospatialTools as wgt
import numpy as np
import matplotlib.pyplot as plt
import os
import json

speedupDictFn = '/Users/wiar9509/Google Drive/wrangells/json/speedupDict_2iterTempFilter_fullyAutomated_28nov2016.json'

data = wgt.readJsonFile(speedupDictFn)
keyIter = data.iterkeys()


# Finding longest glacier for creating matrix of datnums
maxLenDist = 0
for i in range(0,len(data)):
	transName = keyIter.next()
	prof = data[transName]
	dist = prof['termDist']
	
	lenDist = len(dist)
	
	if lenDist > maxLenDist:
		maxLenDist = lenDist
		longestGlacier = transName
		longestTermDist = dist


# Glaciers with defined speedup
goodGlaciersList = ['gateCL','kennCL','wForCL','nabeEB','nabeWB','kluvCL',
					'chisEB','chisWB','kluaCL','donjEB','donjWB','russEB',
					'russWB','barnEB','barnWB','hawkCL','chitNB','chitSB',
					'rootCL']

# Initialize
suToget = np.nan*np.ones((1,maxLenDist))
wiToget = np.nan*np.ones((1,maxLenDist))
suTogetGood = np.nan*np.ones((1,maxLenDist))
wiTogetGood = np.nan*np.ones((1,maxLenDist))

keyIter = data.iterkeys()

# iterate over profiles
for i in range(0,len(data)):
	transName = keyIter.next()
	prof = data[transName]
	filt = prof['tempFilter']
	dNumSu = filt['summerDataNum']
	dNumWi = filt['winterDataNum']
	
	# Adding nan's to pad and make all arrays same length
	padSuNum = maxLenDist-len(dNumSu)
	padWiNum = maxLenDist-len(dNumWi)
	suNans = np.nan*np.ones((1,padSuNum))
	wiNans = np.nan*np.ones((1,padWiNum))	
	
	dNumSuPad = np.append(np.reshape(np.array(dNumSu),(1,len(dNumSu))),suNans,axis=1)
	dNumWiPad = np.append(np.reshape(np.array(dNumWi),(1,len(dNumWi))),wiNans,axis=1)
		
	suToget = np.vstack((suToget,dNumSuPad))
	wiToget = np.vstack((wiToget,dNumWiPad))	
	
	# if glacier has identifiable speedup behavior
	if transName in goodGlaciersList:
		suTogetGood = np.vstack((suTogetGood,dNumSuPad))
		wiTogetGood = np.vstack((wiTogetGood,dNumWiPad))	

		

w25,w50,w75 = np.nanpercentile(wiToget,(25,50,75),axis=0)
s25,s50,s75 = np.nanpercentile(suToget,(25,50,75),axis=0)

wg25,wg50,wg75 = np.nanpercentile(wiTogetGood,(25,50,75),axis=0)
sg25,sg50,sg75 = np.nanpercentile(suTogetGood,(25,50,75),axis=0)

nGlaciers = np.sum(np.isnan(suToget),axis=0)
nGlaciersGood = np.sum(np.isnan(suTogetGood),axis=0)

suTotalGood = np.nansum(suTogetGood,axis=0)
suTotalGoodNorm = suTotalGood/nGlaciersGood[::-1]
suTotal = np.nansum(suToget,axis=0)
suTotalNorm= suTotal/nGlaciers[::-1]

wiTotalGood = np.nansum(wiTogetGood,axis=0)
wiTotalGoodNorm = wiTotalGood/nGlaciersGood[::-1]
wiTotal = np.nansum(wiToget,axis=0)
wiTotalNorm= wiTotal/nGlaciers[::-1]

fig,ax1=plt.subplots()
#ax1.plot(wgt.flipDistance(longestTermDist)/1e3,suToget.transpose(),c='salmon',lw=0,marker='.',markersize=5,alpha=0.1)	
#ax1.plot(wgt.flipDistance(longestTermDist)/1e3,wiToget.transpose(),c='cornflowerblue',lw=0,marker='.',markersize=5,alpha=0.1)	
#ax1.fill_between(wgt.flipDistance(longestTermDist)/1e3,s25,s75,facecolor='salmon',alpha=0.5)
#ax1.fill_between(wgt.flipDistance(longestTermDist)/1e3,w25,w75,facecolor='cornflowerblue',alpha=0.5)
iqr = ax1.fill_between(wgt.flipDistance(longestTermDist)/1e3,sg25,sg75,facecolor='crimson',alpha=0.5,label='Summer IQR')
ax1.fill_between(wgt.flipDistance(longestTermDist)/1e3,wg25,wg75,facecolor='royalblue',alpha=0.5,label='Winter IQR')
w1 = ax1.plot(wgt.flipDistance(longestTermDist)/1e3,s50,c='lightcoral',lw=1,ls='-',label='S$_{50}$ all glaicers')	
s1 = ax1.plot(wgt.flipDistance(longestTermDist)/1e3,w50,c='skyblue',lw=1,ls='-',label='W$_{50}$ all glaicers')
w2 = ax1.plot(wgt.flipDistance(longestTermDist)/1e3,sg50,c='darkred',lw=2,label='Summer median')	
s2 = ax1.plot(wgt.flipDistance(longestTermDist)/1e3,wg50,c='darkblue',lw=2,label='Winter median')
ax2=ax1.twinx()
n2=ax2.plot((np.array(longestTermDist))/1e3,nGlaciers,c='gray',lw=2,label='All glaciers')
n1=ax2.plot((np.array(longestTermDist))/1e3,nGlaciersGood,c='k',lw=2,label='Select glaciers')
ax2.plot((0,np.max(longestTermDist)/1e3),(0,0),lw=1,ls='--',c='k')

ax1.set_yticks(np.arange(0,25,5))
ax2.set_yticks(np.arange(0,65,20))
ax1.set_ylim((4,30))
ax1.set_xlim((0,np.max(longestTermDist)/1e3))
ax2.set_xlim((0,np.max(longestTermDist)/1e3))
ax2.set_ylim((-300,70))
ax1.invert_xaxis()
ax2.invert_xaxis()
ax1.legend(frameon=True,loc='lower center',ncol=3,prop={'size':8})
ax2.legend(frameon=False,loc='upper right',prop={'size':8})
ax1.set_ylabel('Number of observations',fontsize=18)
ax1.set_xlabel('Distance from terminus [km]',fontsize=18)
ax1.yaxis.set_label_coords(-.05,.30)
ax2.set_ylabel('Number of \n glaciers',fontsize=18)
ax2.yaxis.set_label_coords(1.05,.90)
#plt.savefig('averageProfileSpatialCoverage.pdf')
plt.show()
plt.close()

fig.clf()
fig,ax1=plt.subplots()
l1 = ax1.plot(wgt.flipDistance(longestTermDist)/1e3,suTotalGoodNorm,lw=2,c='darkred',label='Summer avg. select')
l2 = ax1.plot(wgt.flipDistance(longestTermDist)/1e3,suTotalNorm,lw=1,c='salmon',label='Summer avg. all')
ax2=ax1.twinx()
l3 = ax2.plot(wgt.flipDistance(longestTermDist)/1e3,suTotalGood,lw=2,c='darkred',ls=':',label='Summer total select')
l4 = ax2.plot(wgt.flipDistance(longestTermDist)/1e3,suTotal,lw=1,c='salmon',ls=':',label='Summer total all')

# winter obs
l1 = ax1.plot(wgt.flipDistance(longestTermDist)/1e3,wiTotalGoodNorm,lw=2,c='darkblue',label='Winter avg. select')
l2 = ax1.plot(wgt.flipDistance(longestTermDist)/1e3,wiTotalNorm,lw=1,c='cornflowerblue',label='Winter avg. all')
l3 = ax2.plot(wgt.flipDistance(longestTermDist)/1e3,wiTotalGood,lw=2,c='darkblue',ls=':',label='Winter total select')
l4 = ax2.plot(wgt.flipDistance(longestTermDist)/1e3,wiTotal,lw=1,c='cornflowerblue',ls=':',label='Winter total all')


ax1.set_ylabel('Average number of observations',fontsize=18)
ax2.set_ylabel('Total number of observations',fontsize=18)
ax1.set_xlabel('Distance from terminus [km]',fontsize=18)
ax1.set_xlim((0,np.max(longestTermDist)/1e3))
ax1.invert_xaxis()
#ax2.invert_xaxis()

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, ncol=2,prop={'size':10},loc='upper center',frameon=False)

#plt.legend((l1,l2,l3,l4),(('Avg. select','Avg. all','Total select','Total all')),loc='upper center',frameon=False)
plt.savefig('averageSpatialCoverage_totalAndAverage.png',dpi=500)
plt.show()
plt.close()

