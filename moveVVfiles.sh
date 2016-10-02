! /bin/sh

# Script to move pycorr velocity files
# William Armstrong
# Original 09 July 2015 - Modified 22 Sep 2015 to accomodate folders with multiple correlation folders

writePath="/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/"

# Loop through all path/row folders in pycorr folder
for folders in /Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/p06*/; do
	echo "Cd to path-row folder: " $folders
	cd $folders
	# Loop through all correlation folders in path/row folder
	for subfolders in LC*/; do
		echo "Cd to image folder: " $subfolders
		cd $subfolders
		
		for corrFolders in LC*24_40_10/; do
			echo "Cd to corr param folder: " $corrFolders
			cd $corrFolders
			# Identify vv (total velocity, think this called euclidean norm too) file
			vvFile=$(ls *vv.tif)
			writeDestination=$writePath$vvFile
			# Copy the vv file to desired folder
			cp $vvFile $writeDestination
			cd .. # Go back to image folder
		done
		cd .. # Go back to path-row folder
	done
	cd .. # Go back to pycorr folder
done

