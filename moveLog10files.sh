! /bin/sh

# Script to move pycorr velocity files
# William Armstrong
# Original 09 July 2015 - modified 22 Sep 2015 to handle multiple correlation folders under image folders.

#writePath="/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/log10_files/"
writePath="/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/log10/"
#folderDirectoryWildcard=/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/p06*/
folderDirectoryWildcard=/Users/anderson/Desktop/ARMSTRONG/wrangells/corr/6*
# Loop through all path/row folders in pycorr folder
for folders in $folderDirectoryWildcard; do
	echo "Cd to path-row folder: " $folders
	cd $folders
	# Loop through all correlation folders in path/row folder
	for subfolders in LC*/; do
		echo "Cd to image folder: " $subfolders
		cd $subfolders
		for corrFolders in LC*24_40_10/; do # The 24_40_10 gets just those folders ending with those corr param flags
			echo "Cd to corr param folder: " $corrFolders
			cd $corrFolders
			# Identify log10 file
			log10File=$(ls *log10.tif)
			# Where to write to
			writeDestination=$writePath$log10File
			# Copy the vv file to desired folder
			cp $log10File $writeDestination 
			cd .. # Go back to Image folder
		done
		cd .. # Go back to path-row folder
	done
	cd .. # Go back to pycorr folder
done

