! /bin/sh

# Script to loop through all landsat imagery in directory and run all available correlations
# 15 September 2015
# William Armstrong

folderPath='/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/p062_r017/band8'

#echo $folderPath
cd $folderPath

image_num=$(ls -1 LC*.TIF | grep -v '/$' | wc -l) # number of l8 images in directory

for img1 in LC*.TIF
	do
	echo 'Master image: ' $img1
	for img2 in LC*.TIF
		do
		echo 'Slave image: ' $img2
		if [ "$img1" == "$img2" ]
		then
			echo 'Same image, not correlating'
		elif [ "$img1" != "$img2" ]
		then
			# Pull date information from filename. The many parens and + 0 makes the string an integer. 
			year1=$((${img1:9:4}+0))
			# Commenting out the line below and switching to current format fixed error:  value too great for base. This apparently because bash thought it was a hex number.
			#doy1=$((${img1:13:3}+0))
			doy1=${img1:13:3}
			year2=$((${img2:9:4}+0))
			doy2=${img2:13:3}
			#doy2=$((${img2:13:3}+0))
			
			
			# Gets the fractional year associated with the doy. This is to save on logic later. Answer from http://stackoverflow.com/questions/1088098/how-do-i-divide-in-the-linux-console
			yearFrac1=$(echo "scale=3; $doy1/365" | bc) 
			yearFrac2=$(echo "scale=3; $doy2/365" | bc)
			
			# This makes a fractional year date for the images
			date1=$(awk "BEGIN {printf \"%.3f\",${year1}+${yearFrac1}}")
			date2=$(awk "BEGIN {printf \"%.3f\",${year2}+${yearFrac2}}")
			
			# afterVal == 1 if image 2 is taken after image 1
			afterVal=`echo "$date2 > $date1" | bc` 
				
			if [ $afterVal -eq 1 ] # If image 2 taken after image 1
			then
				echo 'Correlating master image: ' $img1 ' and slave: ' $img2 
				
				# These lines generate the command for pycorr to correlate $img1 and $img2
				cmdStart="python wpycorr_v1.10.py "
				cmdSpace=" "
				cmdEnd1="-vf -sf -v -imgdir"
				cmdEnd2="-nsidc_out_names -use_hp -gfilt_sigma 3.0 -log10 -inc 10 -half_source_chip 24 -half_target_chip 40"
				cmdString=$cmdStart$img1$cmdSpace$img2$cmdSpace$cmdEnd1$cmdSpace$folderPath$cmdSpace$cmdEnd2
				echo $cmdString
				# Call pycorr for this image pair using above-defined string
				eval $cmdString
			fi
		fi
		done
done