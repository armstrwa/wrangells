import fiona
#import psycopg2
import pprint
import json
import collections  # need ordered dict for schema definition
import numpy as np
from matplotlib import pyplot as plt
import argparse


parser = argparse.ArgumentParser( \
    description= \
    	"""evenly+space_cl_from_shp_v1 reads input shapefile (line) and writes out evenly 
    		spaced point shapefile along the same line - with distance above the first 
    		point in dist field.  Output name is input name plus interval spec
    		
    		present version hard coded to EPSG:3338 - Albers Conformal Conic for Alaska for both
    		input and output shapefiles.  Other projections not dealt with (could be)""",
    formatter_class=argparse.RawDescriptionHelpFormatter)



parser.add_argument('-input_shapefile_dir', 
                    action='store', 
                    type=str, 
                    default='.',
                    help='directory containing profile to use [.]')
# parser.add_argument('-glacier_name', 
#                     action='store', 
#                     type=str, 
#                     default=None,
#                     help='name of glacier in database to pull cl from [None -> use specified shapefile instead]')
parser.add_argument('-key_for_distance_field', 
                    action='store', 
                    type=str,
                    default='dist',
                    help='field identifier in shapefile for point/profile distance [dist]')
parser.add_argument('-del_interval', 
                    action='store', 
                    type=float,
                    default=500.0,
                    help='sampling interval for output shapefile [500]')                                      
parser.add_argument('-output_dir', 
                    action='store', 
                    type=str, 
                    default='.',
                    help='output directory regularly sampled point shapefile [.]')
parser.add_argument('-v', 
                    action='store_true',  
                    default=False, 
                    help='verbose - extra diagnostic and image info [False if not raised]')
parser.add_argument('-input_shapefile', 
                    action='store', 
                    type=str,
                    default=None,
                    help='shapefile to be converted to regularly spaced points [None]')
parser.add_argument('-output_shapefile', 
                    action='store', 
                    type=str,
                    default=None,
                    help='ouput shapefile [None - make it up from glacier name or input shapefile name with interval appended]')
args = parser.parse_args()



print '#', args

print 'looking for centerline shapefile'
clsrid=3338   # forced to be AK...
with fiona.open(args.input_shapefile,'r') as inf:
	in_cl=inf.next()
	in_cl_crs=inf.crs    # but this isn't convertible to epsg...
cl=in_cl['geometry']
lengthorder = 0    # this assumes that input shapefile starts at terminus
glacier_name=args.input_shapefile.replace('.shp','')

# from osgeo import osr
# 
# def esriprj2standards(shapeprj_path):
#    prj_file = open(shapeprj_path, 'r')
#    prj_txt = prj_file.read()
#    srs = osr.SpatialReference()
#    srs.ImportFromESRI([prj_txt])
#    print 'Shape prj is: %s' % prj_txt
#    print 'WKT is: %s' % srs.ExportToWkt()
#    print 'Proj4 is: %s' % srs.ExportToProj4()
#    srs.AutoIdentifyEPSG()
#    print 'EPSG is: %s' % srs.GetAuthorityCode(None)


del_interval=args.del_interval   # meters spacing for regular sampling
if lengthorder==1:
	pts=cl['coordinates'][::-1]    # flip line coordinates end for end so that 0 is the terminus end
else:
	pts=cl['coordinates']
deldiffs=np.diff(pts,axis=0)   # delta x and delta y between adjacent points
deldist=[np.sqrt(x**2 + y**2) for x,y in deldiffs]
cl_dist=np.zeros(len(pts))
cl_dist[1:]=np.cumsum(deldist)  # set to be zero distance at start
mfdel=np.modf(cl_dist/del_interval)
ddint=np.diff(mfdel[1])    # integer part of ratio above - then find steps in that

num_samp_pnts=np.max(mfdel[1])+1    # maximum number of integer del_intervals that fit in cl_dist max
samp_pnts_x=np.zeros(num_samp_pnts)
samp_pnts_y=np.zeros(num_samp_pnts)
samp_pnts_x[0]=pts[0][0]
samp_pnts_y[0]=pts[0][1]
samp_pnts_dist=np.arange(0.0,del_interval * num_samp_pnts,del_interval)
samp_pnts_dist

for i in np.where(ddint>0)[0]:
	numpntstoadd=ddint[i]
	pntindex=mfdel[1][i]    # mfdel[1] is integer part of mod - so the point index for regularly spaced points
	delx=deldiffs[i][0]
	dely=deldiffs[i][1]
	startx=pts[i][0]
	starty=pts[i][1]
	for j in range(1,int(numpntstoadd)+1):
		frac= (deldist[i] - (del_interval*(np.float(numpntstoadd - j) + mfdel[0][i+1])))/deldist[i]    # fraction of total input pnt distance for the jth integer point
		samp_pnts_x[pntindex+j]=startx + (frac*delx)
		samp_pnts_y[pntindex+j]=starty + (frac*dely)
		if args.v:
			print i,j,np.sqrt((samp_pnts_x[pntindex+j]-samp_pnts_x[pntindex+j-1])**2 +\
				(samp_pnts_y[pntindex+j]-samp_pnts_y[pntindex+j-1])**2),\
				[pntindex+j],delx,dely,frac,\
				startx,samp_pnts_x[pntindex+j],frac*delx,starty,samp_pnts_y[pntindex+j],frac*dely
# 
if args.v:
	plt.figure()
	plt.plot(samp_pnts_x,samp_pnts_y,'+')
	plt.ion()
	plt.show()

sdriver=u'ESRI Shapefile'
scrs={u'datum': u'NAD83',
  u'lat_0': 50,
  u'lat_1': 55,
  u'lat_2': 65,
  u'lon_0': -154,
  u'no_defs': True,
  u'proj': u'aea',
  u'units': u'm',
  u'x_0': 0,
  u'y_0': 0}
# sschema={'geometry': 'Point',
#   'properties': OrderedDict([(u'fid', 'int:10'), (u'cng_(Meter', 'int:10')])}
sschema={'geometry': 'Point','properties': collections.OrderedDict([(u'fid', 'int:10'), (u'dist', 'float:16.8')])}

if (args.output_shapefile):
	glacier_name_out=args.output_shapefile
else:
	glacier_name_out=glacier_name.replace(' ','_')+'_%04dm.shp'%(del_interval)
with fiona.open(glacier_name_out,'w',driver=sdriver,crs=scrs,schema=sschema) as c:
    for idn,pnt in enumerate(zip(samp_pnts_x,samp_pnts_y)):
    	feature_to_write={'geometry': {'coordinates': pnt,
  					'type': 'Point'},
 					'id': str(idn),
 					'properties': collections.OrderedDict([(u'fid', idn), (u'dist', samp_pnts_dist[idn])]),
 					'type': 'Feature'}
        c.write(feature_to_write)
#         if args.v:
#         	print len(c)
print 'wrote shapefile: %s'%glacier_name_out

