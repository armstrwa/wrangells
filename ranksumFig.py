ranksumFn = '/Users/wiar9509/Google Drive/wrangells/json/ranksumSpeedupDict_30nov2016 2.json'

ranksumData = wgt.readJsonFile(ranksumFn)

prof = ranksumData['nabeEB']
termDist = prof['termDist']
medDiff = prof['medSpeedup']
speedupLength = prof['speedupLength']
breakPoint = np.min(np.where(np.array(termDist)<=speedupLength))
wList = prof['rankSumStat']
pList = prof['rankSumPval']

plt.subplot(2,1,1)
plt.plot(np.array(termDist)/1e3,medDiff,lw=1,c='gray')
plt.plot(np.array((speedupLength,speedupLength))/1e3,(-0.5,0.5),ls='--',c='k')
plt.plot((0,np.max(termDist)/1e3),(0,0),ls='--',c='k')
upGlacierInd = np.arange(0,breakPoint+1)
downGlacierInd = np.arange(breakPoint+1,len(termDist))
plt.scatter(np.array(termDist)[upGlacierInd]/1e3,np.array(medDiff)[upGlacierInd],c='b',s=20)
plt.scatter(np.array(termDist)[downGlacierInd]/1e3,np.array(medDiff)[downGlacierInd],c='r',s=20)
plt.ylabel('Speedup [m d$^{-1}$]',fontsize=18)
plt.xlim((0,np.max(termDist)/1e3))
plt.ylim((-.4,.4))
plt.gca().invert_xaxis()
#plt.title(transNow + ' rank sum for speedup length')
plt.subplot(2,1,2)
plt.plot(np.array(termDist)/1e3,wList,lw=2)
plt.ylabel('Test statistic, W$_{rs}$ [-]',fontsize=18)
plt.xlim((0,np.max(termDist)/1e3))
plt.plot(np.array((speedupLength,speedupLength))/1e3,(-8,1),ls='--',c='k')
# plt.subplot(3,1,3)
# plt.plot(np.array(termDist)/1e3,pList,lw=2)
# plt.ylabel('p val')
plt.xlabel('Distance from terminus [km]',fontsize=18)
plt.xlim((0,np.max(termDist)/1e3))
plt.gca().invert_xaxis()
plt.savefig('ranksumForSpeedupLengthFigure.pdf')
plt.show()
plt.close()
