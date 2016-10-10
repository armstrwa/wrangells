## SCRIPT TO ITERATE RASTER REPROJECTION
#  Written by William Armstrong
#  13 November 2015
#
#  With help from: http://stackoverflow.com/questions/8715419/passing-a-variable-set-by-the-user-as-a-parameter-for-a-os-system-command
#  and http://gis.stackexchange.com/questions/92907/re-project-raster-image-from-mercator-to-equirectangular

#inDirectory = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/'
inDirectory = '/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/select/filtered/'
#outDirectory = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/annualCorrelations/'
#srsIn = 32607 # 32607 = utm w/ wgs84
srsOut = 102006 # 102006 = alaska albers conic conformal

daysBetweenToggle = 0 # turn to 1 to use days between images as an iteration filter
minDaysBetween = 0 # don't process if days between is less than this number
maxDaysBetween = 600 # don't process if days between is greater than this number


### PROCESSING ###

# Loading modules 
import os
import rasterio

allFiles = os.listdir(inDirectory)

for fileNow in allFiles:
	if os.path.isfile(fileNow) and fileNow.endswith('.tif') and fileNow[0:3] == 'LC8' and fileNow.split('_')[-1][0:4] != 'EPSG':
	
		rast = rasterio.open(inDirectory + fileNow)
		srsIn = int(rast.crs.to_dict()['init'][5:])

		
		fn1, fn2, centerDoy, daysBtwn = fileNow.split('_')[0:4] # parse correlation name
		wrsPath = fn1[3:6]
		wrsRow = fn1[6:9]
		
		# Days between threshold filter
		if daysBetweenToggle == 1:
			if int(daysBtwn) >= minDaysBetween and int(daysBtwn) <= maxDaysBetween:
				fnIn = inDirectory + fileNow
				fnOut = outDirectory + fileNow[0:-4] + '_EPSG' + str(srsOut) + '.tif'
				if not os.path.isfile(fnOut):
					os.system("gdalwarp -s_srs EPSG:" +  str(srsIn) + " -t_srs EPSG:" + str(srsOut) + " " + fnIn + " " + fnOut )
		else:
				fnIn = inDirectory + fileNow
				fnOut = outDirectory + fileNow[0:-4] + '_EPSG' + str(srsOut) + '.tif'
				if not os.path.isfile(fnOut):
					os.system("gdalwarp -s_srs EPSG:" +  str(srsIn) + " -t_srs EPSG:" + str(srsOut) + " " + fnIn + " " + fnOut )
					
					
					
					
					
					
					
								