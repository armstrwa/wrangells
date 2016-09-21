### SCRIPT TO ITERATE MARK'S L8 SAMPLER CODE
#   William Armstrong
#   17 Nov 2015

import os
import glob


shapeInDir='/Users/wiar9509/Documents/generalScripts/swath/shapefileInputs/'
myDataToggle = 1 # 1 for my data, 0 for Mark data

if myDataToggle == 1:
	rasterInDir='/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006'
	#rasterInDir='/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/annualCorrelations/'
	fnStart = 'LC8'
	outDir = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput'
	#outDir='/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/annualCorrelations/L8SamplerOutput'
	scriptPath = '/Users/wiar9509/Documents/generalScripts/corr/L8_sample_frames_for_shape_pts_v4_forAK.py'
	maxSpeed = 2.
elif myDataToggle == 0:
	rasterInDir = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/markL8correlations'
	fnStart = 'AK_ACC_dc0p2_akv0_'		
	outDir = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/markL8correlations/L8SamplerOutput/saturate_0.5'
	scriptPath = '/Users/wiar9509/Documents/generalScripts/fromMarkFahnestock/L8_sample_frames_for_shape_pts_v4_whaModified17nov2015.py'
	maxSpeed = 0.5
#rasterInDir='/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006'

#fnStart = 'LC8'

#outDir = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered/EPSG102006/L8SamplerOutput'


## Processing

filesIn = glob.glob(shapeInDir + '*_evenlySpace_100m.shp')
os.chdir(outDir) # so the output is written here

for file in filesIn:
	shapefileName = file.split('/')[-1]
	transName = shapefileName.split('_')[0]
	#commandText = 'python /Users/wiar9509/Documents/generalScripts/corr/L8_sample_frames_for_shape_pts_v4_forAK.py -key_for_distance_field dist -WRS2_shapefile_path /Users/wiar9509/Documents/generalScripts/fromMarkFahnestock/wrs2_descending/wrs2_descending.shp -indir /Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered -fn_start LC8 -plt_speed_max 2.0 -profile_shapefile_dir ' + shapeInDir + ' ' + shapefileName
	commandText='python ' + scriptPath + ' -indir ' + rasterInDir + ' -WRS2_shapefile_path /Users/wiar9509/Documents/generalScripts/fromMarkFahnestock/wrs2_descending/wrs2_descending.shp -fn_start ' + fnStart + ' -profile_shapefile_dir /Users/wiar9509/Documents/generalScripts/swath/shapefileInputs/ -plt_speed_max ' + str(maxSpeed) + ' ' + shapefileName
	#commandText = 'python /Users/wiar9509/Documents/generalScripts/fromMarkFahnestock/L8_sample_frames_for_shape_pts_v4.py -key_for_distance_field dist -WRS2_shapefile_path /Users/wiar9509/Documents/generalScripts/fromMarkFahnestock/wrs2_descending/wrs2_descending.shp -indir /Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/markL8correlations -fn_start AK_ACC_dc0p2_akv0_ -output_dir ' + outDir + ' -plt_speed_max 2.0 -profile_shapefile_dir ' + shapeInDir + ' ' + shapefileName
	os.system(commandText)
