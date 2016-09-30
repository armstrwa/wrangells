'''
Make RGB low-res preview image from separate bands of landsat imagery
'''

import os
import glob

def makeRGBfromBands(folderPath,outName,tr = 'auto', bs=[4,3,2]):
	'''
	Make color composite from three separate bands
	Inputs: folderPath = path containing image bands, tr = 100 100 = target resolution for 100 m pixels, bs = (4,3,2) = bands to use for R,G,B
	Outputs: tiff at outName
	'''
	# Remove output files if they exist 
	if os.path.isfile(outName):
		print "Output filename already exists; deleting"
		os.remove(outName)
	
	if os.path.isfile('tmp.tif'):
		os.remove('tmp.tif')
	
	# Get filenames for bands to be used for R,G,B
	rFn = glob.glob(folderPath + '*B' + str(bs[0]) + '.TIF' or folderPath + '*B' + str(bs[0]) + '.tif')[0]
	gFn = glob.glob(folderPath + '*B' + str(bs[1]) + '.TIF' or folderPath + '*B' + str(bs[1]) + '.tif')[0]
	bFn = glob.glob(folderPath + '*B' + str(bs[2]) + '.TIF' or folderPath + '*B' + str(bs[2]) + '.tif')[0]

	print "Red : " + rFn
	print "Green : " + gFn
	print "Blue : " + bFn
	
	if tr == 'auto': # assume output same resolution of not specified
		mergeCommandString = "gdal_merge.py -separate -co " + '"' + 'COMPRESS=LZW' + '"' + " -o " + "'" + outName  + "'" + " " + rFn + " " + gFn + " " + bFn
	else: # if target resolution specified
		mergeCommandString = "gdal_merge.py -separate -co " + '"' + 'COMPRESS=LZW' + '"' + " -ps " + str(tr) + " " + str(tr) + " -o " + "'" + 'tmp.tif'  + "'" + " " + rFn + " " + gFn + " " + bFn
	
	# run gdal_merge.py to make multiband
	os.system(mergeCommandString)
	
	# run gdal_translate to make RGB interpretation
	translateCommandString = "gdal_translate -of 'GTiff' -co " + '"' + 'PHOTOMETRIC=RGB' + '" -scale 0 100 -a_nodata 0 tmp.tif ' + outName
	os.system(translateCommandString)
	
	print "Created file: " + outName



folderPath = '/Volumes/oldEthan/l8_imagery/imagery/'
outPath = '/Volumes/oldEthan/l8_imagery/imagery/falseColorComposite/'
				
folderList = os.listdir(folderPath)

for folder in folderList:
	if os.path.isdir(folderPath + folder):
		subfolderList = os.listdir(folderPath + folder)
		for subfolder in subfolderList:
			if os.path.isdir(folderPath + folder + '/' + subfolder):
				subPath = folderPath + folder + '/' + subfolder + '/'

				bs = (7,5,3)
				outName = outPath + subfolder + '_rgb' + str(bs[0]) + str(bs[1]) + str(bs[2]) + '.tif'
				makeRGBfromBands(subPath,outName,250,bs)
				
				
				
