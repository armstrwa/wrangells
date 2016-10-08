import matplotlib.pyplot as plt
import numpy as np
import json
import matplotlib.gridspec as gridspec
from osgeo import gdal, ogr, osr
from gdalconst import *

# Bounding box to pixel offsets for gdal.ReadAsArray()
# from https://gist.github.com/perrygeo/5667173
def bbox_to_pixel_offsets(gt, bbox):
    originX = gt[0]
    originY = gt[3]
    pixel_width = gt[1]
    pixel_height = gt[5]
    x1 = int((bbox[0] - originX) / pixel_width)
    x2 = int((bbox[1] - originX) / pixel_width) + 1

    y1 = int((bbox[3] - originY) / pixel_height)
    y2 = int((bbox[2] - originY) / pixel_height) + 1

    xsize = x2 - x1
    ysize = y2 - y1
    return (x1, y1, xsize, ysize)

# Python-implemented zonal stats
# from https://gist.github.com/perrygeo/5667173
def zonal_stats(vector_path, raster_path, nodata_value=None, global_src_extent=False):
    rds = gdal.Open(raster_path, GA_ReadOnly)
    assert(rds)
    rb = rds.GetRasterBand(1)
    rgt = rds.GetGeoTransform()

    if nodata_value:
        nodata_value = float(nodata_value)
        rb.SetNoDataValue(nodata_value)

    vds = ogr.Open(vector_path, GA_ReadOnly)  # TODO maybe open update if we want to write stats
    assert(vds)
    vlyr = vds.GetLayer(0)

    # create an in-memory numpy array of the source raster data
    # covering the whole extent of the vector layer
    if global_src_extent:
        # use global source extent
        # useful only when disk IO or raster scanning inefficiencies are your limiting factor
        # advantage: reads raster data in one pass
        # disadvantage: large vector extents may have big memory requirements
        src_offset = bbox_to_pixel_offsets(rgt, vlyr.GetExtent())
        src_array = rb.ReadAsArray(*src_offset)

        # calculate new geotransform of the layer subset
        new_gt = (
            (rgt[0] + (src_offset[0] * rgt[1])),
            rgt[1],
            0.0,
            (rgt[3] + (src_offset[1] * rgt[5])),
            0.0,
            rgt[5]
        )

    mem_drv = ogr.GetDriverByName('Memory')
    driver = gdal.GetDriverByName('MEM')

    # Loop through vectors
    stats = []
    feat = vlyr.GetNextFeature()
    while feat is not None:

        if not global_src_extent:
            # use local source extent
            # fastest option when you have fast disks and well indexed raster (ie tiled Geotiff)
            # advantage: each feature uses the smallest raster chunk
            # disadvantage: lots of reads on the source raster
            src_offset = bbox_to_pixel_offsets(rgt, feat.geometry().GetEnvelope())
            src_array = rb.ReadAsArray(*src_offset)

            # calculate new geotransform of the feature subset
            new_gt = (
                (rgt[0] + (src_offset[0] * rgt[1])),
                rgt[1],
                0.0,
                (rgt[3] + (src_offset[1] * rgt[5])),
                0.0,
                rgt[5]
            )

        # Create a temporary vector layer in memory
        mem_ds = mem_drv.CreateDataSource('out')
        mem_layer = mem_ds.CreateLayer('poly', None, ogr.wkbPolygon)
        mem_layer.CreateFeature(feat.Clone())

        # Rasterize it
        rvds = driver.Create('', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
        rvds.SetGeoTransform(new_gt)
        gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1])
        rv_array = rvds.ReadAsArray()

        # Mask the source data array with our current feature
        # we take the logical_not to flip 0<->1 to get the correct mask effect
        # we also mask out nodata values explictly
        masked = np.ma.MaskedArray(
            src_array,
            mask=np.logical_or(
                src_array == nodata_value,
                np.logical_not(rv_array)
            )
        )

        feature_stats = {
            'min': float(masked.min()),
            'mean': float(masked.mean()),
            'max': float(masked.max()),
            'std': float(masked.std()),
            'sum': float(masked.sum()),
            'count': int(masked.count()),
            'fid': int(feat.GetFID())}

        stats.append(feature_stats)

        rvds = None
        mem_ds = None
        feat = vlyr.GetNextFeature()

    vds = None
    rds = None
    return stats


def readZonalStatsOutput(zonalStatsOut):
	'''
	# Function to put zonal stats output into an array
	# Ouputs a dictionary of statistics
	'''
	
	boxNum = len(zonalStatsOut) # number of boxes over which statistics were calculated

	# intialize
	mins = []
	maxs = []
	means = []
	counts = []
	stds = []
	fids = []

	# Read stats and append to list
	for i in range(0,boxNum):
		boxNow = zonalStatsOut[i]
		fids.append(boxNow['fid'])	
		mins.append(boxNow['min'])
		means.append(boxNow['mean'])
		maxs.append(boxNow['max'])
		counts.append(boxNow['count'])	
		stds.append(boxNow['std'])
		
#	statsOut = np.array( (fids,mins,means,maxs,stds,counts) ).transpose() # stick together to write out
	statsOut = {'fids':fids,'means':means,'mins':mins,'maxs':maxs,'stds':stds,'counts':counts}	
	return statsOut


def calculateSeasonalVelocityStats(uSeason):		
	'''
	Function to calculate season velocity percentiles		
	Takes in seasonal velocity (here called uSu b/c just copied from below)
	Ouputs dataOut = (x, u10, u25, u50, u75, u90)
	Something going wrong at 75th and 90th percentiles where just all nan. Think might have to do with treatment of nan
	
	Taken from old code, probably sloppy 26sep16 wha
	'''
	# initialize
	season10 = None
	season25 = None
	season50 = None
	season75 = None
	season90 = None
	
	# identify missing data and turn to nan
	nanInd = np.where(uSeason==-1)
	uSeason[nanInd] = np.nan
	
	if len(uSeason.shape)>1:
		xLen = uSeason.shape[1]
	
		for i in range(0,xLen):
			dataInd = np.isnan(uSeason[:,i]) != True
			colData = uSeason[dataInd,i]
			if len(colData) == 0:
				season10x = np.nan
				season25x = np.nan
				season50x = np.nan				
				season75x = np.nan				
				season90x = np.nan
			else:		
				season10x, season25x, season50x, season75x, season90x = np.percentile(colData,[10,25,50,75,90])
	
			if i == 0:
				season10 = season10x
				season25 = season25x
				season50 = season50x						
				season75 = season75x
				season90 = season90x
			else:
				season10 = np.column_stack((season10,season10x))	
				season25 = np.column_stack((season25,season25x))	
				season50 = np.column_stack((season50,season50x))	
				season75 = np.column_stack((season75,season75x))	
				season90 = np.column_stack((season90,season90x))																		
	else: # if not enough observations then nan
		season10 = np.empty((1,len(uSeason)))
		season10.fill(np.nan)
		season25 = np.copy(season10)
		season50 = np.copy(season10)		
		season75 = np.copy(season10)
		season90 = np.copy(season10)
		
	# Append data
	dataOut = np.row_stack((season10,season25,season50,season75,season90))
	
	
	return dataOut


def getNumDataFromJson(jsonData,startSeasons=[304, 81, 121, 265]):
	'''
	Function that records number of observations within each season
	Inputs: jsonData = json file with velocity data; startSeasons = [winter, spring, summer, fall] start season doy
	Outputs: dataTogether = [x, nWinter, nSpring, nSummer, nFall]
	# start days of each season, these based of local met data (304 = Oct 28; 81 = Mar 22; 121 = May 1; 265 = Sep 22)
	From old script, probably bad and clunky 26sep16 wha
	'''
	shapefileName = jsonData['profile_shapefile']
	#transName = shapefileName.split('_')[0]
	transName = jsonData['transect_name']
	x = np.array(jsonData['sample_pts_frontdist']) # distance from headwall [m]
	profiles = jsonData['profiles'] # getting all profiles
	numProfiles=len(profiles)

	dataNumAll = np.zeros((len(x)))
	dataNumWi = np.copy(dataNumAll)
	dataNumSp = np.copy(dataNumAll)
	dataNumSu = np.copy(dataNumAll)
	dataNumFa = np.copy(dataNumAll)

	for i in range(0,numProfiles):
		profNow = profiles[i]
		startDecYear = profNow['start_dec_year']
		midDecYear = profNow['mid_dec_year']
		endDecYear = profNow['stop_dec_year']
	
		middleDoy = (midDecYear - np.floor(midDecYear))*365.25
	
		u = np.array(profNow['speed'])
		dataInd=np.where(u!=-1)
		dataNumAll[dataInd] = dataNumAll[dataInd] + 1
	
		# Categorize by season
		if int(middleDoy) >= startSeasons[0] or int(middleDoy) < startSeasons[1]:
			seasonColor = 'c' # winter color
			midSeason = 'winter'
			dataNumWi[dataInd] = dataNumWi[dataInd] + 1
		elif int(middleDoy) >= startSeasons[1] and int(middleDoy) < startSeasons[2]:
			seasonColor = 'y' # spring color
			midSeason = 'spring'
			dataNumSp[dataInd] = dataNumSp[dataInd] + 1
		elif int(middleDoy) >= startSeasons[2] and int(middleDoy) < startSeasons[3]:
			seasonColor = 'g' # summer color
			midSeason = 'summer'
			dataNumSu[dataInd] = dataNumSu[dataInd] + 1
		elif int(middleDoy) >= startSeasons[3] and int(middleDoy) < startSeasons[0]:
			seasonColor = 'r' # fall color
			midSeason = 'fall'
			dataNumFa[dataInd] = dataNumFa[dataInd] + 1
	
	dataTogether = np.array((x,dataNumAll,dataNumWi,dataNumSp,dataNumSu,dataNumFa))
	return dataTogether.transpose()


def computeStatsAndPlot(fnIn,pathOut,startSeasons=[304, 81, 121, 265],axis='auto'):
	'''
	Script to read JSON velocity file, perform seasonal statistics, and plot results
	Input: json filename to process
	Output: data = x, nAll, uAll, nWinter, uWinter ... 
	Old and clunky
	
	Added NDSI 05oct2016 wha
	'''

	with open(fnIn) as data_file:
		data = json.load(data_file)

	# Read in JSON data to see how many observations at each x
	obsNum = getNumDataFromJson(data)

	shapefileName = data['profile_shapefile']
	transName = data['transect_name']
	x = np.array(data['sample_pts_frontdist']) # distance from headwall [m]
	profiles = data['profiles'] # getting all profiles
	numProfiles=len(profiles)

	uMatx = np.empty(len(x))
	uWi = []
	uSp = []
	uSu = []
	uFa = []

	for i in range(0,numProfiles): # iterate over profiles
		u = np.array(profiles[i]['speed'])
		if i == 0:
			uMatx = u
		else:
			uMatx = np.column_stack((uMatx,u))
		startDecYear = np.array(profiles[i]['start_dec_year'])
		startYear = np.floor(startDecYear)
		startDoy=np.round((startDecYear-startYear)*365)
		middleDecYear = np.array(profiles[i]['mid_dec_year'])
		middleYear = np.floor(middleDecYear)
		middleDoy = np.round((middleDecYear-middleYear)*365)
		endDecYear = np.array(profiles[i]['stop_dec_year'])
		endYear = np.floor(endDecYear)
		endDoy=np.round((endDecYear-endYear)*365)
		daysBtwn = np.array(profiles[i]['delts_days'])
	
		ndsi1 = np.array(profiles[i]['image1_snowIndex'])
		ndsi2 = np.array(profiles[i]['image2_snowIndex'])
		
		# Plot point
		#Start season based
		# Identify start seasons

		if int(middleDoy) >= startSeasons[0] or int(middleDoy) < startSeasons[1]:
			seasonColor = 'c' # winter color
			midSeason = 'winter'
		elif int(middleDoy) >= startSeasons[1] and int(middleDoy) < startSeasons[2]:
			seasonColor = 'y' # spring color
			midSeason = 'spring'
		elif int(middleDoy) >= startSeasons[2] and int(middleDoy) < startSeasons[3]:
			seasonColor = 'g' # summer color
			midSeason = 'summer'
		elif int(middleDoy) >= startSeasons[3] and int(middleDoy) < startSeasons[0]:
			seasonColor = 'r' # fall color
			midSeason = 'fall'
		
		# Compiling seasonal velocity
		if midSeason == 'winter':
			if len(uWi) == 0:
				uWi = u
			else:
				uWi = np.row_stack((uWi,u))
		elif midSeason == 'spring':
			if len(uSp) == 0:
				uSp = u
			else:
				uSp = np.row_stack((uSp,u))				
		elif midSeason == 'summer':
			if len(uSu) == 0:
				uSu = u
			else:
				uSu = np.row_stack((uSu,u))			
		elif midSeason == 'fall':
			if len(uFa) == 0:
				uFa = u
			else:
				uFa = np.row_stack((uFa,u))	
			

		#ax1.plot(x/1e3,u,linewidth=0,marker='.',color=seasonColor,alpha=0.5)

	# If no correlations in that season, then nan
	if len(uWi) == 0:
		uWi = np.empty((1,len(x)))
		uWiStats = np.empty((5,len(x)))
		uWi.fill(np.nan)
		uWiStats.fill(np.nan)
	else: #otherwise calculate stats
		uWiStats = calculateSeasonalVelocityStats(uWi)
	if len(uSp) == 0:
		uSp = np.empty((1,len(x)))
		uSpStats = np.empty((5,len(x)))		
		uSp.fill(np.nan)
		uSpStats.fill(np.nan)
	else:
		uSpStats = calculateSeasonalVelocityStats(uSp)			
	if len(uSu) == 0:
		uSu = np.empty((1,len(x)))
		uSuStats = np.empty((5,len(x)))			
		uSu.fill(np.nan)
		uSuStats.fill(np.nan)
	else:
		uSuStats = calculateSeasonalVelocityStats(uSu)		
	if len(uFa) == 0:
		uFa = np.empty((1,len(x)))
		uFaStats = np.empty((5,len(x)))	
		uFa.fill(np.nan)
		uFaStats.fill(np.nan)
	else:
		uFaStats = calculateSeasonalVelocityStats(uFa)			
						
	uAllStats = calculateSeasonalVelocityStats(uMatx.transpose())






	if len(uWiStats.shape)<=1: # this true if not enough data to properly calculate	
		uWiStats = np.empty((5,len(x)))
		uWiStats[:] = np.nan
	if len(uSpStats.shape)<=1: # this true if not enough data to properly calculate	
		uSpStats = np.empty((5,len(x)))
		uSpStats[:] = np.nan
	if len(uSuStats.shape)<=1: # this true if not enough data to properly calculate	
		uSuStats = np.empty((5,len(x)))
		uSuStats[:] = np.nan	
	if len(uFaStats.shape)<=1: # this true if not enough data to properly calculate	
		uFaStats = np.empty((5,len(x)))
		uFaStats[:] = np.nan			
	
	# Difference from median			
	uWiDiff50 = uWiStats[2,:] - uAllStats[2,:]
	uSpDiff50 = uSpStats[2,:] - uAllStats[2,:]
	uSuDiff50 = uSuStats[2,:] - uAllStats[2,:]
	uFaDiff50 = uFaStats[2,:] - uAllStats[2,:]				

	# Greater than 50th/75th %

	uWiGt50 = uWiDiff50 >= 0
	uSpGt50 = uSpDiff50 >= 0
	uSuGt50 = uSuDiff50 >= 0	
	uFaGt50 = uFaDiff50 >= 0	
	
	uWiGt75 = uWiStats[2,:] >= uAllStats[3,:]
	uSpGt75 = uSpStats[2,:] >= uAllStats[3,:]	
	uSuGt75 = uSuStats[2,:] >= uAllStats[3,:]	
	uFaGt75 = uFaStats[2,:] >= uAllStats[3,:]	

	# Less than 25ths
	uWiLt25 = uWiStats[2,:] <= uAllStats[1,:]
	uSpLt25 = uSpStats[2,:] <= uAllStats[1,:]
	uSuLt25 = uSuStats[2,:] <= uAllStats[1,:]
	uFaLt25 = uFaStats[2,:] <= uAllStats[1,:]

	# Filter out seasons with too little data
	obsNumThreshold = 3
	wiBadObsInd = np.where(obsNum[:,2]<obsNumThreshold)
	uWi50plot = np.copy(uWiStats[2,:])
	uWi50plot[wiBadObsInd] = np.nan

	spBadObsInd = np.where(obsNum[:,3]<obsNumThreshold)
	uSp50plot = np.copy(uSpStats[2,:])
	uSp50plot[spBadObsInd] = np.nan	

	suBadObsInd = np.where(obsNum[:,4]<obsNumThreshold)
	uSu50plot = np.copy(uSuStats[2,:])
	uSu50plot[suBadObsInd] = np.nan		

	faBadObsInd = np.where(obsNum[:,5]<obsNumThreshold)
	uFa50plot = np.copy(uFaStats[2,:])
	uFa50plot[faBadObsInd] = np.nan		

	# Set up plot
	gs = gridspec.GridSpec(4,1,height_ratios=[1,3,1,1])
	ax1=plt.subplot(gs[1])
	ax2=plt.subplot(gs[0])
	ax3=plt.subplot(gs[3])
	ax4=plt.subplot(gs[2])

	if 1:
		# all observations
		ax1.plot(x/1e3,uWi.transpose(),marker='.',color='c',lw=0,alpha=0.25)
		ax1.plot(x/1e3,uSp.transpose(),marker='.',color='y',lw=0,alpha=0.25)
		ax1.plot(x/1e3,uSu.transpose(),marker='.',color='g',lw=0,alpha=0.25)
		ax1.plot(x/1e3,uFa.transpose(),marker='.',color='r',lw=0,alpha=0.25)

	# all season IQR and median
	ax1.plot(x/1e3,uAllStats[1,:],lw=1,color='k')
	ax1.plot(x/1e3,uAllStats[2,:],lw=2,color='k',label='Median: all data')	
	ax1.plot(x/1e3,uAllStats[3,:],lw=1,color='k',label='25th-75th %: all data')
	ax1.fill_between(x/1e3,	uAllStats[1,:],uAllStats[3,:],color='gray',alpha=0.1)

	# seasonal medians
	if 1:
		ax1.plot(x/1e3,uWiStats[2,:],lw=2,color='c')
		ax1.plot(x/1e3,uSpStats[2,:],lw=2,color='y')
		ax1.plot(x/1e3,uSuStats[2,:],lw=2,color='g')	
		ax1.plot(x/1e3,uFaStats[2,:],lw=2,color='r')	
	if 0:
		ax1.plot(x/1e3,uWi50plot,lw=2,color='c')
		ax1.plot(x/1e3,uSp50plot,lw=2,color='orange')
		ax1.plot(x/1e3,uSu50plot,lw=2,color='g')
		ax1.plot(x/1e3,uFa50plot,lw=2,color='y')		

	# Seasonal difference from all data median
	ax4.plot(x/1e3,uWiDiff50,lw=2,ls='-',color='c')
	ax4.plot(x/1e3,uSpDiff50,lw=2,ls='-',color='y')
	ax4.plot(x/1e3,uSuDiff50,lw=2,ls='-',color='g')
	ax4.plot(x/1e3,uFaDiff50,lw=2,ls='-',color='r')			
	ax4.plot(x/1e3,np.zeros(len(x)),lw=1,ls='--',color='k')

	xLims4 = [-.3, .3]
	xMax4 = xLims4[1]
	xMin4 = xLims4[0]
	trueLoc1 = 0.02
	trueLoc2 = 0.04
	trueLoc3 = 0.06
	trueLoc4 = 0.08
		
	# Where seasonal velocity is gt 75th percentile
	ax4.plot(x[uWiGt75==True]/1e3,uWiGt75[uWiGt75==True]*(xMax4 - trueLoc1),lw=0,marker='o',color='c',alpha=0.5)
	ax4.plot(x[uSpGt75==True]/1e3,uSpGt75[uSpGt75==True]*(xMax4 - trueLoc2),lw=0,marker='o',color='y',alpha=0.5)
	ax4.plot(x[uSuGt75==True]/1e3,uSuGt75[uSuGt75==True]*(xMax4 - trueLoc3),lw=0,marker='o',color='g',alpha=0.5)
	ax4.plot(x[uFaGt75==True]/1e3,uFaGt75[uFaGt75==True]*(xMax4 - trueLoc4),lw=0,marker='o',color='r',alpha=0.5)		

	# Where is velocity lt 25th percentile
	ax4.plot(x[uWiLt25==True]/1e3,uWiLt25[uWiLt25==True]*(xMin4 + trueLoc1),lw=0,marker='o',color='c',alpha=0.5)
	ax4.plot(x[uSpLt25==True]/1e3,uSpLt25[uSpLt25==True]*(xMin4 + trueLoc2),lw=0,marker='o',color='y',alpha=0.5)
	ax4.plot(x[uSuLt25==True]/1e3,uSuLt25[uSuLt25==True]*(xMin4 + trueLoc3),lw=0,marker='o',color='g',alpha=0.5)
	ax4.plot(x[uFaLt25==True]/1e3,uFaLt25[uFaLt25==True]*(xMin4 + trueLoc4),lw=0,marker='o',color='r',alpha=0.5)

	saveData = np.vstack((x,obsNum[:,1],uAllStats[1,:],uAllStats[2,:],uAllStats[3,:],obsNum[:,2],uWiStats[2,:],obsNum[:,3],uSpStats[2,:],obsNum[:,4],uSuStats[2,:],obsNum[:,5],uFaStats[2,:]))
	saveDataT = saveData.transpose()
	np.savetxt(pathOut + transName +'_velocityStats.csv',saveDataT,delimiter=',')
		 
					
	## GET ELEVATION PROFILE
	#elevData = np.genfromtxt(elevationProfileFolder + shapefileName[:-4] + '_elevationProfile.csv',delimiter=',')
	if(0):
		elevData = np.genfromtxt('/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/dtm/elevationProfiles/' + transName + '_elevationProfile.csv',delimiter=',')
		ax3.plot(elevData[:,0]/1e3,elevData[:,1],color='gray')		

	ax3.plot(x/1e3,ndsi1,lw=2,alpha=0.7,c='b')
	ax3.plot(x/1e3,ndsi2,lw=2,alpha=0.7,c='b')	

	ax2.plot(obsNum[:,0]/1e3,obsNum[:,1],marker='.',color='k')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,2],marker='.',color='c')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,3],marker='.',color='y')	
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,4],marker='.',color='g')
	ax2.plot(obsNum[:,0]/1e3,obsNum[:,5],marker='.',color='r')		
				


	ax3.set_xlabel('Distance from headwall [km]',fontsize=16)
	ax1.set_ylabel('Velocity [m d$^{-1}$]',fontsize=16)
	ax2.set_title(transName,fontsize=18)

	ax1.plot(-1,-1,marker='.',color='c',linewidth=0,label='Winter')
	ax1.plot(-1,-1,marker='.',color='y',linewidth=0,label='Spring')
	ax1.plot(-1,-1,marker='.',color='g',linewidth=0,label='Summer')
	ax1.plot(-1,-1,marker='.',color='r',linewidth=0,label='Fall')


	#ax3.set_ylim([np.nanmin(elevData[:,1]),	np.nanmax(elevData[:,1])])
	ax2.set_ylabel('# Observations',fontsize=18)
	ax2.yaxis.set_label_position("right")
	ax3.set_ylabel('NDSI [-]',fontsize=18)
	#ax2.set_yticks(np.arange( round(elevData[:,1].min()/500)*500,round(elevData[:,1].max()/500)*500+500,500))
	ax3.set_yticks((np.arange(-0.25,1.25,0.25)))
	#ax1.set_ylim([0,np.nanmin((2.5,np.nanmax(uAllStats[4,:])))]) # auto-scaling limits
	ax1.set_ylim([0,2])
	ax1.set_xlim([0,x.max()/1e3])
	ax1lims = ax1.xaxis.get_view_interval()
	ax3.set_xlim([ax1lims[0],ax1lims[1]])
	ax2.set_xlim([ax1lims[0],ax1lims[1]])
	ax4.set_xlim([ax1lims[0],ax1lims[1]])
	ax4.set_ylim(xLims4)
	ax4.set_ylabel('Diff. from annual median',fontsize=18)
	ax4.yaxis.set_label_position("right")
	
	ax1.legend(loc='best',prop={'size':6},numpoints=1)
	plt.draw()
	plt.savefig(pathOut + transName + '_profilePlot_wSeasonalMedian_wElevation_wDataNum_wDiffFromMean.pdf')
	#plt.show()
	plt.close()		
		
	return saveData



fnIn = '/Users/anderson/kennCL_swathVelocitySampling_2016-10-05.json'
rastIn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/LC80640172016191LGN0_LC80640172016239LGN0_2016-215_48_24_40_10_vv_filt_corr0.15_delcorr0.05.tif'
summer = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/LC80640172016159LGN0_LC80640172016239LGN0_2016-199_80_24_40_10_vv_filt_corr0.15_delcorr0.05.tif'
spring = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/LC80640172016159LGN0_LC80640172016175LGN0_2016-167_16_24_40_10_vv_filt_corr0.15_delcorr0.05.tif'
winter = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/LC80640172014057LGN0_LC80640172014089LGN0_2014-073_32_24_40_10_vv_filt_corr0.15_delcorr0.05.tif'
#shpIn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/kennCL_samplePolygons.shp'
shpIn = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/swath/rootCL_samplePolygons.shp'
pathOut = 'test.pdf'
startSeasons=[304, 81, 121, 265]
axis='auto'

stats = zonal_stats(shpIn,summer,nodata_value=-1)
summerStats = readZonalStatsOutput(stats)

stats = zonal_stats(shpIn,winter,nodata_value=-1)
winterStats = readZonalStatsOutput(stats)

stats = zonal_stats(shpIn,spring,nodata_value=-1)
springStats = readZonalStatsOutput(stats)


ax1 = plt.subplot(2,1,1)
ax1.plot(np.array(winterStats['fids'])*500,winterStats['means'],lw=2,c='b',label='Winter: 2014 DOYs 057-089')
ax1.plot(np.array(winterStats['fids'])*500,springStats['means'],lw=2,c='c',label='Summer: 2016 DOYs 159-175')
ax1.plot(np.array(winterStats['fids'])*500,summerStats['means'],lw=2,c='r',label='Summer: 2016 DOYs 159-239')
plt.title('Root Glacier velocity variability',fontsize=20)
plt.legend(loc=1)
ax2 = plt.subplot(2,1,2)
ax2.plot(np.array(winterStats['fids'])*500,winterStats['counts'],lw=2,ls='--',c='b')
ax2.plot(np.array(winterStats['fids'])*500,summerStats['counts'],lw=2,ls='--',c='r')
ax2.plot(np.array(winterStats['fids'])*500,springStats['counts'],lw=2,ls='--',c='c')
plt.xlabel('Down-glacier distance [m]',fontsize=16)
ax1.set_ylabel('Velocity [m/d]',fontsize=16)
ax2.set_ylabel('Number of obs [-]',fontsize=16)
plt.savefig('rootGlacier_velocity_summerWinter.pdf')

plt.show()
plt.close()





