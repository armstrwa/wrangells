def intToBinary(intArr,bitNum=16):
	'''
	Convert an MxN integer array to an MxNxB bit array
	Inputs: intArr = integer array; bitNum = number of bits
	Outputs: outArr = bit array, where bits are now in 3rd dimension
	'''

	try:
		len(intArr) # this will balk if only one number
		
		# Get array shape
		rows,cols = np.shape(intArr)

		# Flatten array to next step doesn't balk
		flatArr = intArr.flatten()

		# Convert flattened integer array into bit array
		# answer from http://stackoverflow.com/questions/22227595/convert-integer-to-binary-array-with-suitable-padding
		bitArr = (((flatArr[:,None] & (1 << np.arange(bitNum)))) > 0).astype(int)
	
		# Reshape bitArr to be shape = [rows,cols,bitNum], so raster has original shape with bits in 3rd dimension
		outArr = np.reshape(bitArr,(rows,cols,bitNum))
				
	except TypeError: # if only one number
		outArr = "{0:b}".format(intArr)

	return outArr

def convertDoubleArrBitToIntArr(doubleBitArr):
	'''
	Takes landsat's double-bit info and converts into integer
	Input: double bit band (shape = row by col by 2 depth)
	Output: row by col array where: -1 = no, 0 = not determined, 1 = maybe, 2 = yes
	'''
	rows,cols,bits = doubleBitArr.shape
	unknown = np.logical_and(doubleBitArr[:,:,0] == 0, doubleBitArr[:,:,1] == 0)
	no = np.logical_and(doubleBitArr[:,:,0] == 0, doubleBitArr[:,:,1] == 1)
	maybe = np.logical_and(doubleBitArr[:,:,0] == 1, doubleBitArr[:,:,1] == 0)	
	yes = np.logical_and(doubleBitArr[:,:,0] == 1, doubleBitArr[:,:,1] == 1)
	
	intArrOut = np.zeros((rows,cols))
	intArrOut[no] = -1
	intArrOut[maybe] = 1
	intArrOut[yes] = 2
	
	return intArrOut
			
def binaryArrayToQualityDesignations(bitArr,designations='all'):
	'''
	Convert row by col by 16 bit depth array into landsat quality designations
	Input: bitArr = row by col by 16 bit depth array; designations='all','cloudsTerrainSnow' for which bits to return
	Output = outDict = dictionary where keys contain rasters of yes/maybe/no/unknown. For 2 bit items, -1 = no, 0 = unknown, 1 = maybe, 2 = yes. For 1 bit items (terrainOcclusion,droppedFrame,designatedFill), 1 = yes, 0 = no
	Modified from function named qualityAssessmentBand on https://github.com/poldrosky/alternar/blob/95b78fa0b8f4961caeb3b75ceb825349ae3f11ff/DNtoToAReflectanceL8/DNtoToAReflectanceL8.py
	'''
	valueMeanings = {'0':'No', '1':'Yes', '00':'Not Determined', '01':'No', '10':'Maybe', '11':'Yes'}

# 	if designations =='all':
# 		cloud = convertDoubleArrBitToIntArr(bitArr[:,:,0:2])
# 		cirrus = convertDoubleArrBitToIntArr(bitArr[:,:,2:4])
# 		snowIce = convertDoubleArrBitToIntArr(bitArr[:,:,4:6])
# 		vegetation = convertDoubleArrBitToIntArr(bitArr[:,:,6:8])
# 		cloudShadow = convertDoubleArrBitToIntArr(bitArr[:,:,8:10])
# 		water = convertDoubleArrBitToIntArr(bitArr[:,:,10:12])
# 		reserved = bitArr[:,:,12:13]
# 		terrainOcclusion = bitArr[:,:,13:14]
# 		droppedFrame = bitArr[:,:,14:15]
# 		designatedFill = bitArr[:,:,15:16]
# 	
# 		outDict = {'key':valueMeanings,'cloud':cloud, 'cirrus':cirrus,'snowIce':snowIce,'veg':vegetation,"cloudShadow":cloudShadow,"water":water,"reserved":reserved,"terrainOcclusion":terrainOcclusion,"droppedFrame":droppedFrame,"fill":designatedFill }			
# 	
# 	elif designations == 'cloudsTerrainSnow':
# 		cloud = convertDoubleArrBitToIntArr(bitArr[:,:,0:2])
# 		cirrus = convertDoubleArrBitToIntArr(bitArr[:,:,2:4])
# 		snowIce =convertDoubleArrBitToIntArr( bitArr[:,:,4:6])
# 		cloudShadow = convertDoubleArrBitToIntArr(bitArr[:,:,8:10])
# 		terrainOcclusion = convertDoubleArrBitToIntArr(bitArr[:,:,13:14])
# 		outDict = {'key':valueMeanings,'cloud':cloud, 'cirrus':cirrus,'snowIce':snowIce,"cloudShadow":cloudShadow,"terrainOcclusion":terrainOcclusion }			

	if designations =='all':
		cloud = convertDoubleArrBitToIntArr(bitArr[:,:,14:16])
		cirrus = convertDoubleArrBitToIntArr(bitArr[:,:,12:14])
		snowIce = convertDoubleArrBitToIntArr(bitArr[:,:,10:12])
		vegetation = convertDoubleArrBitToIntArr(bitArr[:,:,8:10])
		cloudShadow = convertDoubleArrBitToIntArr(bitArr[:,:,6:8])
		water = convertDoubleArrBitToIntArr(bitArr[:,:,4:6])
		reserved = bitArr[:,:,3:4]
		terrainOcclusion = bitArr[:,:,2:3]
		droppedFrame = bitArr[:,:,1:2]
		designatedFill = bitArr[:,:,0:1]
	
		outDict = {'key':valueMeanings,'cloud':cloud, 'cirrus':cirrus,'snowIce':snowIce,'veg':vegetation,"cloudShadow":cloudShadow,"water":water,"reserved":reserved,"terrainOcclusion":terrainOcclusion,"droppedFrame":droppedFrame,"fill":designatedFill }			
	
	elif designations == 'cloudsTerrainSnow':
		cloud = convertDoubleArrBitToIntArr(bitArr[:,:,0:2])
		cirrus = convertDoubleArrBitToIntArr(bitArr[:,:,2:4])
		snowIce =convertDoubleArrBitToIntArr( bitArr[:,:,4:6])
		cloudShadow = convertDoubleArrBitToIntArr(bitArr[:,:,8:10])
		terrainOcclusion = convertDoubleArrBitToIntArr(bitArr[:,:,13:14])
		outDict = {'key':valueMeanings,'cloud':cloud, 'cirrus':cirrus,'snowIce':snowIce,"cloudShadow":cloudShadow,"terrainOcclusion":terrainOcclusion }			


	return outDict
	
def makeQAplot(qaDict,plotToggle=True,saveToggle=True,outName='tmp.tif'):
	cloudMa = np.ma.masked_where(qaDict['cloud']<2,qaDict['cloud'])
	cirrusMa = np.ma.masked_where(qaDict['cirrus']<1,qaDict['cirrus'])
	snowMa = np.ma.masked_where(qaDict['snowIce']<1,qaDict['snowIce'])	
	terrainMa = np.ma.masked_where(qaDict['terrainOcclusion']<1,qaDict['terrainOcclusion'])		
	vegMa = np.ma.masked_where(qaDict['veg']<1,qaDict['veg'])
	waterMa =  np.ma.masked_where(qaDict['water']<1,qaDict['water'])

	#plt.figure()
	plt.imshow(cloudMa,cmap='autumn')
	plt.imshow(snowMa,cmap='cool_r')
	plt.imshow(vegMa,cmap='YlGn')
	plt.show()
	plt.close()

qaFn = '/Volumes/oldEthan/l8_imagery/imagery/62_17/LC80620172014219LGN00/LC80620172014219LGN00_BQA.TIF'
qa = gdal.Open(qaFn)
qaArr = qa.ReadAsArray(2500,4000,1000,1000)
bitArr = intToBinary(qaArr)
qaDict = binaryArrayToQualityDesignations(bitArr)

cloudMa = np.ma.masked_where(qaDict['cloud']<2,qaDict['cloud'])
cirrusMa = np.ma.masked_where(qaDict['cirrus']<1,qaDict['cirrus'])
snowMa = np.ma.masked_where(qaDict['snowIce']<1,qaDict['snowIce'])	
terrainMa = np.ma.masked_where(qaDict['terrainOcclusion']<1,qaDict['terrainOcclusion'])		
vegMa = np.ma.masked_where(qaDict['veg']<1,qaDict['veg'])
waterMa =  np.ma.masked_where(qaDict['water']<1,qaDict['water'])

#plt.figure()
plt.imshow(cloudMa,cmap='cool_r')
plt.imshow(snowMa,cmap='cool')
plt.imshow(vegMa,cmap='YlGn')
plt.imshow(waterMa,cmap='GnBu')
#plt.imshow(terrainMa,cmap='autumn_r')
plt.show()
plt.close()






