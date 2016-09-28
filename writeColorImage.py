'''
Messing with Mark's scripts to get rgb image for Wrangells analysis
'''

import matplotlib.pyplot as plt
from osgeo import gdal
import numpy as np
import os

rasterFn = '/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/LC80630172014050LGN0_LC80630172014290LGN0_2014-170_240_24_40_10_hp_filt_3.0_vv.tif'
plotvmax = 2

rast = gdal.Open(rasterFn)
speed = rast.ReadAsArray()
dir_name, file_name = os.path.split(rasterFn)
file_name_base = file_name[:-4]
vv = np.ma.masked_where(speed<0.001,speed)

#vv=np.ma.masked_where(((del_corr_arr<dcam)&(corr_arr<cam))|(corr_arr<cam1)|(((offset_dist_ij_arr*img1.pix_x_m)/del_t_val)>plotvmax),((offset_dist_ij_arr*img1.pix_x_m)/del_t_val))
mapjet=plt.get_cmap('jet')
mmvv=vv.copy()
mmvv[(mmvv>plotvmax) & (~mmvv.mask)]=plotvmax
vv_zo=mmvv/plotvmax
vv_rgba=mapjet(vv_zo) # this produces a 4-band 0.0-1.0 image array using the jet colormap
#print vv_rgba.shape
(out_lines,out_pixels,out_bands)=vv_rgba.shape

# prepare a log10 version of the speed for output below
mmvv=vv.copy()
mmvv[(mmvv>plotvmax) & (~mmvv.mask)]=plotvmax
mmvv.mask[(mmvv==0) & (~mmvv.mask)]= True
lmmvv=np.log10(mmvv)

min_lmmvv=np.log10(0.01)
max_lmmvv=np.log10(plotvmax)
range_lmmvv=max_lmmvv - min_lmmvv
lvv_zo=(lmmvv - min_lmmvv)/(range_lmmvv)
lvv_rgba=mapjet(lvv_zo) # this produces a 4-band 0.0-1.0 image array using the jet colormap

format = "GTiff"
driver = gdal.GetDriverByName( format )
dst_filename = file_name_base + '_log10.tif'
dst_ds = driver.Create( dst_filename, out_pixels, out_lines, out_bands, gdal.GDT_Byte )
print 'out image %s %s %d, %d, %d'%(dst_filename, format, out_pixels, out_lines, out_bands)

# dst_ds.SetGeoTransform( [ com_min_x, inc * img1.pix_x_m, 0, com_max_y, 0, inc * img1.pix_y_m ] ) # note pix_y_m typically negative
dst_ds.SetGeoTransform( rast.GetGeoTransform() ) # note pix_y_m typically negative
dst_ds.SetProjection( rast.GetProjection() )
dst_ds.GetRasterBand(1).WriteArray( (lvv_rgba[:,:,0]*255).astype('ubyte') )
dst_ds.GetRasterBand(2).WriteArray( (lvv_rgba[:,:,1]*255).astype('ubyte') )
dst_ds.GetRasterBand(3).WriteArray( (lvv_rgba[:,:,2]*255).astype('ubyte') )
dst_ds.GetRasterBand(4).WriteArray( (lvv_rgba[:,:,3]*255).astype('ubyte') )
dst_ds = None # done, close the dataset
