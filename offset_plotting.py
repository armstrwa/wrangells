from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt

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
					
	print transNow, "{0:.3f}".format(speedMagListDown[i]), "{0:.2E}".format(speedSlopeListDown[i])
	
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
	plt.savefig('speedupSquiggles_withFitLines_noFill_16feb2017.pdf')
	plt.show()
	plt.close()



if 1:

    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)

    par1 = host.twinx()
    par2 = host.twinx()

    offX = 60
    offY = 0
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right",
                                        axes=par2,
                                        offset=(offX,offY))

    par2.axis["right"].toggle(all=True)

    host.set_xlim(0, 2)
    host.set_ylim(0, 2)

    host.set_xlabel("Distance")
    host.set_ylabel("Density")
    par1.set_ylabel("Temperature")
    par2.set_ylabel("Velocity")

    p1, = host.plot([0, 1, 2], [0, 1, 2], label="Density")
    p2, = par1.plot([0, 1, 2], [0, 3, 2], label="Temperature")
    p3, = par2.plot([0, 1, 2], [50, 30, 15], label="Velocity")

    par1.set_ylim(0, 4)
    par2.set_ylim(1, 65)

    host.legend()