import shapely.geometry
import fiona
from pyproj import Proj, transform
import glob
import os
import osr
import gdal
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker
import datetime as dt
import argparse
import json
import time
from operator import itemgetter


class GeoImg:
	"""geocoded image input and info
		a=GeoImg(in_file_name,indir='.')
			a.img will contain image
			a.parameter etc..."""
	def __init__(self, in_filename,in_dir='.',datestr=None,datefmt='%m/%d/%y'):
		self.filename = in_filename
		self.in_dir_path = in_dir  #in_dir can be relative...
		self.in_dir_abs_path=os.path.abspath(in_dir)  # get absolute path for later ref if needed
		self.gd=gdal.Open(self.in_dir_path + os.path.sep + self.filename)
		self.srs=osr.SpatialReference(wkt=self.gd.GetProjection())
		self.gt=self.gd.GetGeoTransform()
		self.proj=self.gd.GetProjection()
		self.intype=self.gd.GetDriver().ShortName
		self.min_x=self.gt[0]
		self.max_x=self.gt[0]+self.gd.RasterXSize*self.gt[1]
		self.min_y=self.gt[3]+self.gt[5]*self.gd.RasterYSize
		self.max_y=self.gt[3]
		self.pix_x_m=self.gt[1]
		self.pix_y_m=self.gt[5]
		self.num_pix_x=self.gd.RasterXSize
		self.num_pix_y=self.gd.RasterYSize
		self.XYtfm=np.array([self.min_x,self.max_y,self.pix_x_m,self.pix_y_m]).astype('float')
		if (datestr is not None):
			self.imagedatetime=dt.datetime.strptime(datestr,datefmt)
		elif ((self.filename.find('LC8') == 0) | (self.filename.find('LE7') == 0) | (self.filename.find('LT5') == 0) | (self.filename.find('LT4') == 0)):	# looks landsat like - try parsing the date from filename (contains day of year)
			self.sensor=self.filename[0:3]
			self.path=int(self.filename[3:6])
			self.row=int(self.filename[6:9])
			self.year=int(self.filename[9:13])
			self.doy=int(self.filename[13:16])
			self.imagedatetime=dt.date.fromordinal(dt.date(self.year-1,12,31).toordinal()+self.doy)
		else:
			self.imagedatetime=None  # need to throw error in this case...or get it from metadata
		self.img=self.gd.ReadAsArray().astype(np.float32)   # works for L8 and earlier - and openCV correlation routine needs float or byte so just use float...
# 		self.img=self.gd.ReadAsArray().astype(np.uint8)		# L7 and earlier - doesn't work with plt.imshow...
		self.img_ov2=self.img[0::2,0::2]
		self.img_ov10=self.img[0::10,0::10]
		self.srs=osr.SpatialReference(wkt=self.gd.GetProjection())
	def imageij2XY(self,ai,aj,outx=None,outy=None):
		it = np.nditer([ai,aj,outx,outy],
						flags = ['external_loop', 'buffered'],
						op_flags = [['readonly'],['readonly'],
									['writeonly', 'allocate', 'no_broadcast'],
									['writeonly', 'allocate', 'no_broadcast']])
		for ii,jj,ox,oy in it:
			ox[...]=(self.XYtfm[0]+((ii+0.5)*self.XYtfm[2]));
			oy[...]=(self.XYtfm[1]+((jj+0.5)*self.XYtfm[3]));
		return np.array(it.operands[2:4])
	def XY2imageij(self,ax,ay,outi=None,outj=None):
		it = np.nditer([ax,ay,outi,outj],
						flags = ['external_loop', 'buffered'],
						op_flags = [['readonly'],['readonly'],
									['writeonly', 'allocate', 'no_broadcast'],
									['writeonly', 'allocate', 'no_broadcast']])
		for xx,yy,oi,oj in it:
			oi[...]=((xx-self.XYtfm[0])/self.XYtfm[2])-0.5;  # if python arrays started at 1, + 0.5
			oj[...]=((yy-self.XYtfm[1])/self.XYtfm[3])-0.5;  # " " " " "
		return np.array(it.operands[2:4])


########################################################
# following from accepted answer here: http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years
# but I import datetime package as dt, not datetime.datetime as dt, so this is modfied to reflect that
########################################################

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt.datetime(year=year, month=1, day=1)
    startOfNextYear = dt.datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction
    



parser = argparse.ArgumentParser( \
    description="""L8_sample_frames_for_shape_pts_Greenland_v3 """,
    epilog='>>  <<',
    formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('-WRS2_shapefile_path', 
                    action='store', 
                    type=str, 
                    default='/Users/mark/Library/Mobile Documents/com~apple~CloudDocs/L8_list/wrs2_descending.shp',
                    help='reference Path,Row polygons from USGS [/Users/mark/Library/Mobile Documents/com~apple~CloudDocs/L8_list/wrs2_descending.shp]')
parser.add_argument('-profile_shapefile_dir', 
                    action='store', 
                    type=str, 
                    default='.',
                    help='directory containing profile to use [.]')
parser.add_argument('-key_for_distance_field', 
                    action='store', 
                    type=str,
                    default='dist',
                    help='field identifier in shapefile for point/profile distance [dist]')
parser.add_argument('-key_for_pointID', 
                    action='store', 
                    type=str,
                    default=None,
                    help='field identifier in shapefile for point id [None]')                    
parser.add_argument('-key_for_lineID', 
                    action='store', 
                    type=str,
                    default=None,
                    help='field identifier in shapefile for line id [None]')  
parser.add_argument('-Line_ID_value_to_plot', 
                    action='store', 
                    type=str,
                    default=None,
                    help='line id to use for plotting (if key_for_lineID specified, good idea to chose one line)[None]')                   
parser.add_argument('-indir', 
                    action='store', 
                    type=str,
                    default='/Volumes/tb2r8_data/L8/Greenland_out/grn_wpc0/repro_PS_2',
                    help='directory of reprojected velocity fields to scan[/Volumes/tb2r8_data/L8/Greenland_out/grn_wpc0/repro_PS_2]')
parser.add_argument('-fn_start', 
                    action='store', 
                    type=str, 
                    default='PS_OC_grn2_',
                    help='start of filenames for reprojected velocity fields [PS_OC_grn2_]')
parser.add_argument('-output_dir', 
                    action='store', 
                    type=str, 
                    default='.',
                    help='output directory for json and png files [.]')
parser.add_argument('-plt_year_min', 
                    action='store', 
                    type=float, 
                    default='2013.0',
                    help='min year on plots [2013.0]')
parser.add_argument('-plt_year_max', 
                    action='store', 
                    type=float, 
                    default='2016.0',
                    help='max year on plots [2016.0]')
parser.add_argument('-plt_speed_min', 
                    action='store', 
                    type=float, 
                    default='0.0',
                    help='min speed on plots [0.0]')
parser.add_argument('-plt_speed_max', 
                    action='store', 
                    type=float, 
                    default=None,
                    help='max speed on plots [None - set from data ceil(1.25 * 99.5 percentile value)]')
parser.add_argument('-plt_dist_min', 
                    action='store', 
                    type=float, 
                    default='0.0',
                    help='min dist on plots [0.0]')
parser.add_argument('-plt_dist_max', 
                    action='store', 
                    type=float, 
                    default=None,
                    help='max distance on plots [None - take from profile data]')
parser.add_argument('-plt_numpts_timeprofiles', 
                    action='store', 
                    type=int, 
                    default=15,
                    help='number of point timeseries to pick from profile points for middle plot [15]')
parser.add_argument('-v', 
                    action='store_true',  
                    default=False, 
                    help='verbose - extra diagnostic and image info [False if not raised]')
parser.add_argument('-nofigs', 
                    action='store_true',  
                    default=False, 
                    help='don\'t show figures (still output png) [False if not raised]')
parser.add_argument('profile_shapefile', 
                    action='store', 
                    type=str,
                    default=None,
                    help='shapefile with points to sample[None]')
args = parser.parse_args()


year_min=args.plt_year_min
year_max=args.plt_year_max
speed_min=args.plt_speed_min
speed_max=args.plt_speed_max
dist_min=args.plt_dist_min
dist_max=args.plt_dist_max
plt_numpts_timeprofiles=args.plt_numpts_timeprofiles   # use every fourth point along profile for middle plot



# bash script lines to fetch and process google list and process manually fetched usgs earthexplorer list
# datestr=`date +%Y_%m_%d`
# echo $datestr
# unzip scene_list.zip | grep L8 > L8_scene_list_$datestr
# unzip -c scene_list.zip | grep L8 > L8_scene_list_$datestr
# head -2 L8_scene_list_$datestr
# cp `ls -lrt LANDSAT_8_*.csv | tail -1 | awk '{print $9}'` LANDSAT_8_$datestr.csv

print '#', args




# from fiona.crs import from_epsg

# point = shapely.geometry.Point(-50.0,69.15)
# wrs_shapefile_dir='/Users/mark/Library/Mobile Documents/com~apple~CloudDocs/L8_list'
# wrs_shapefile='wrs2_descending.shp'
# wrs_shapefile='wrs2_asc_desc/wrs2_asc_desc.shp'
# 
# profile_shapefile_dir='/Users/mark/Library/Mobile Documents/com~apple~CloudDocs/L8_list'
# profile_shapefile='Rink_cl_1km_PS.shp'

# profile_shapefile_dir='/Volumes/tb2r8_data/L8/Greenland_out/profiles/Jak_profiles'
profile_shapefile_dir=args.profile_shapefile_dir
# profile_shapefile='helheim.shp'
# profile_shapefile='iker_north.shp'
# profile_shapefile='iker_south.shp'
# profile_shapefile='midgaard.shp'
# profile_shapefile='cl2_chain_1km.shp'
# profile_shapefile='Rink_cl_500m.shp'
profile_shapefile=args.profile_shapefile
# key_for_distance_field='Value'
key_for_distance_field=args.key_for_distance_field

indir=args.indir
# indir='/Volumes/tb2r8_data/L8/Greenland_JakVel/PS_v2'
# fn_start='PS_OC_grn2_'
fn_start=args.fn_start
# use fiona to read in points along centerline in local projection (called PS here)
# save the PS original points in sample_pts_PS, then reproject to 4326, push back into
# coordinates of input Point object (f) and make a list of these (out_pts_ll) and when done save
# them as a Shapely MultiPoint object
with fiona.open(profile_shapefile_dir+'/'+profile_shapefile,'r') as rcl_fc:
	original = Proj(rcl_fc.crs)
	destination = Proj(init='EPSG:4326')
	sample_pts_PS=[]
	sample_pts_lon_lat=[]
	sample_pts_frontdist=[]
	if (args.key_for_pointID):
		sample_pts_point_id=[]
	if (args.key_for_lineID):
		sample_pts_line_id=[]
	out_pts_ll=[]
	for f in rcl_fc:
		x,y =  f['geometry']['coordinates']
		sample_pts_PS.append((x,y))
# 		sample_pts_frontdist.append(f['properties']['cng_(Meter'])   #'cng_(Meter' a holdover from chainage QGIS plugin - will need to fix
		sample_pts_frontdist.append(float(f['properties'][key_for_distance_field]))
		if (args.key_for_pointID):
			sample_pts_point_id.append(f['properties'][args.key_for_pointID])
		if (args.key_for_lineID):
			sample_pts_line_id.append(f['properties'][args.key_for_lineID])
		lon,lat = transform(original, destination,x,y)
		sample_pts_lon_lat.append((lon,lat))
		# 		print (long,lat,x,y)
		# change only the coordinates of the feature
		f['geometry']['coordinates'] = (lon,lat)
		out_pts_ll.append(f)

points_ll = shapely.geometry.multipoint.MultiPoint([shapely.geometry.shape(f['geometry']) for f in out_pts_ll])

front_dist=[x['properties'][key_for_distance_field] for x in out_pts_ll]

# now run through the WRS2 outlines and find path,rows that intersect with the MultiPoint object
pr_list=[]
wrs2_fc = fiona.open(args.WRS2_shapefile_path)
for f in wrs2_fc:
    shape=shapely.geometry.asShape( f['geometry'] )
    if shape.intersects(points_ll):
        print f['properties']['PATH'],f['properties']['ROW']
        pr_list.append((f['properties']['PATH'],f['properties']['ROW']))


total_files=0
all_frames=[]
for path,row in pr_list:
	infiles=glob.glob(indir+'/'+fn_start+'S8_%03d_%03d_*vv.tif'%(path,row))
	print path,row,' number of frames: %d'%(len(infiles))
	frame_info={'path':path, 'row':row, 'infiles':infiles}
	all_frames.append(frame_info)
	total_files+=len(infiles)

# color_max=total_files
# color_index=range(color_max)
# scaled_z = np.array(color_index) / float(color_max)
# colors = plt.cm.jet(scaled_z)
# color_count=0
profiles=[]
for frame in all_frames:
	path=frame['path']
	row=frame['row']
	for inf in frame['infiles']:
		filepath,part,filename=inf.rpartition('/')
		fn_without_start=filename[len(fn_start):]
		series0,path_str,row_str,deldays_str,year1_str,doy1_str,year2_str,doy2_str=fn_without_start.split('_')[0:8]
		series_tag=fn_start+series0+'_'+path_str+'_'+row_str
		# should really comment out the next two lines - but check to see result from above is same?
		path=int(path_str)
		row=int(row_str)
		year1=int(year1_str)
		########################
		doy1=int(doy1_str)
		year2=int(year2_str)
		doy2=int(doy2_str)
		dt_1=dt.datetime.fromordinal(dt.date(year1-1,12,31).toordinal()+doy1)
		dt_2=dt.datetime.fromordinal(dt.date(year2-1,12,31).toordinal()+doy2)
		timedelta=dt_2-dt_1
		delts_days=timedelta.days
		mid_dt=dt_1+(timedelta/2)
		start_dec_year=toYearFraction(dt_1.date())
		stop_dec_year=toYearFraction(dt_2.date())
		mid_dec_year=toYearFraction(mid_dt.date())
		
		ina=GeoImg(filename,in_dir=filepath)
		ptsx,ptsy=zip(*sample_pts_PS)
		ptsi,ptsj=ina.XY2imageij(ptsx,ptsy)
		outvv=-1.0*np.ones_like(ptsx)
		valid_pts_index=np.where((np.rint(np.array(ptsi))<np.array(ina.num_pix_x)) & (np.array(ptsi)>=0) & (np.rint(np.array(ptsj))<ina.num_pix_y) & (np.array(ptsj)>=0))
		outvv[valid_pts_index]=ina.img[np.rint(ptsj[valid_pts_index]).astype(int),np.rint(ptsi[valid_pts_index]).astype(int)]
		if (args.v & (len(valid_pts_index[0])<len(ptsx))):
			print 'path %d row %d %d %s dropped %d points'%(path,row,len(valid_pts_index[0]),filename,len(ptsx)-len(valid_pts_index[0]))
		profile_dict={ 'delts_days':delts_days,\
						'series':series_tag, 'start_dec_year':start_dec_year, 'stop_dec_year':stop_dec_year, 'mid_dec_year':mid_dec_year,\
						'speed':outvv.tolist(), 'path':path, 'row':row, 'infile':filename, 'indir':filepath }
		profiles.append(profile_dict)
		ina=None


#convert from profiles on individual days to time series at each point
delts_days=[]
indir=[]
infile=[]
start_dec_year=[]
mid_dec_year=[]
stop_dec_year=[]
path=[]
row=[]
series=[]
pt_speed=[[] for a in range(len(sample_pts_PS))]  # create a list of blank lists to append the speeds of each point into...
for prof in profiles:
	delts_days.append(prof['delts_days'])
	indir.append(prof['indir'])
	infile.append(prof['infile'])
	start_dec_year.append(prof['start_dec_year'])
	mid_dec_year.append(prof['mid_dec_year'])
	stop_dec_year.append(prof['stop_dec_year'])
	path.append(prof['path'])
	row.append(prof['row'])
	series.append(prof['series'])
	for indx,sp in enumerate(prof['speed']):
		pt_speed[indx].append(sp)

pts={ 'delts_days':delts_days, 'indir':indir, 'infile':infile, 'start_dec_year':start_dec_year, 'mid_dec_year':mid_dec_year,\
		'stop_dec_year':stop_dec_year, 'path':path, 'row':row, 'series':series, 'pt_speed':pt_speed }

current_time=dt.datetime.now()

if (args.key_for_pointID):
	if (args.key_for_lineID):
		profile_data={ 'sample_pts_lon_lat':sample_pts_lon_lat, 'sample_pts_PS':sample_pts_PS, 'sample_pts_frontdist':sample_pts_frontdist,\
					'sampled_on_date_iso':current_time.isoformat(), 'profile_shapefile':profile_shapefile, 'profile_shapefile_dir':profile_shapefile_dir,\
					'profile_shapefile_projection_Proj':original.srs, 'profiles':profiles, 'pts':pts, \
					'sample_pts_point_id':sample_pts_point_id, 'sample_pts_line_id':sample_pts_line_id }
	else:
		profile_data={ 'sample_pts_lon_lat':sample_pts_lon_lat, 'sample_pts_PS':sample_pts_PS, 'sample_pts_frontdist':sample_pts_frontdist,\
					'sampled_on_date_iso':current_time.isoformat(), 'profile_shapefile':profile_shapefile, 'profile_shapefile_dir':profile_shapefile_dir,\
					'profile_shapefile_projection_Proj':original.srs, 'profiles':profiles, 'pts':pts, \
					'sample_pts_point_id':sample_pts_point_id }
else:	
	profile_data={ 'sample_pts_lon_lat':sample_pts_lon_lat, 'sample_pts_PS':sample_pts_PS, 'sample_pts_frontdist':sample_pts_frontdist,\
					'sampled_on_date_iso':current_time.isoformat(), 'profile_shapefile':profile_shapefile, 'profile_shapefile_dir':profile_shapefile_dir,\
					'profile_shapefile_projection_Proj':original.srs, 'profiles':profiles, 'pts':pts }

#  note that numpy arrays are not serializable in json (just lists) so speed arrays sampled for each profile have been converted to lists
#  above using outvv.tolist()...  just be aware that is was an array, but will come back as a list.
json_outfile='%s_profiles_sampled_%s.json'%(profile_shapefile.split('.')[0],time.strftime("%Y-%m-%d"))
with open(args.output_dir + '/' + json_outfile, 'w') as f: 
    f.write(json.dumps(profile_data))
# here is how you might get this back:
#  with open('Rink_cl_1km_PS_profiles_sampled_04-09-2015.json') as data_file:
#    	data = json.load(data_file)



if not(speed_max):
	all_speeds=[]
	for spds in profile_data['pts']['pt_speed']:
		all_speeds.extend(spds)
	speed_max=np.ceil(1.25 * np.percentile(all_speeds,99.5))  # set from data

if not(dist_max):   # max distance for plot not specified, find from input profile
	dist_max=np.max(profile_data['sample_pts_frontdist'])
	
	
fig4=plt.figure(figsize=[26,16])
# year_min=2013.0
# year_max=2016.0
# speed_min=0.0
# speed_max=7.0
# dist_min=0.0
# dist_max=50000.0
# skippnts_timeprofiles=8   # use every fourth point along profile for middle plot
if ((dist_max!=None) | (dist_min!=None)):
	if ((dist_min == None) & (dist_max!=None)):
		pt_indx_to_plot=range(len(np.where(np.array(profile_data['sample_pts_frontdist'])<dist_max)[0]))
	elif ((dist_min != None) & (dist_max!=None)):
		pt_indx_to_plot=range(len(np.where(np.array(profile_data['sample_pts_frontdist'])<dist_min)[0]),len(np.where(np.array(profile_data['sample_pts_frontdist'])<dist_max)[0]))
	elif ((dist_min != None) & (dist_max==None)):
		pt_indx_to_plot=range(len(np.where(np.array(profile_data['sample_pts_frontdist'])<dist_min)[0]),len(profile_data['sample_pts_frontdist']))
else:
	pt_indx_to_plot=range(len(profile_data['sample_pts_frontdist']))
	
if(args.Line_ID_value_to_plot):  # don't plot all points, only points from one line
	pt_indx_to_plot=np.where(np.array(sample_pts_line_id)==args.Line_ID_value_to_plot)[0]

color_max=year_max
color_min=year_min
# scaled_z = (np.array(middates) - color_min) / float(color_max - color_min)
# colors = plt.cm.jet(scaled_z)
cm=plt.cm.get_cmap('jet')
# fig0=plt.figure(figsize=[20,16])
ax0=fig4.add_subplot('131',title='profiles (color is profile date)',ylabel='speed (m/d)',xlabel='postition along profile (m)')
ax0.set_autoscaley_on(False)
ax0.set_ylim([speed_min, speed_max])
ax0.set_autoscalex_on(False)
ax0.set_xlim([dist_min, dist_max])

color_count=0
for sp_profile in profiles:
# 	sc=plt.scatter(front_dist,sp_profile['speed'],marker='+',c=np.array(mid_dec_year[color_count]).repeat(len(front_dist)),s=40.0, vmin=color_min, vmax=color_max, cmap=cm)
	sc=plt.scatter(np.array(front_dist)[pt_indx_to_plot],np.array(sp_profile['speed'])[pt_indx_to_plot],marker='+',c=np.array(mid_dec_year[color_count]).repeat(len(pt_indx_to_plot)),s=40.0, vmin=color_min, vmax=color_max, cmap=cm)
	color_count+=1
cb=plt.colorbar(sc,format='%6.1f')
for item in ([ax0.title, ax0.xaxis.label, ax0.yaxis.label] +
             ax0.get_xticklabels() + ax0.get_yticklabels() + cb.ax.get_yticklabels()):
    item.set_fontsize(17)

if not(args.nofigs):
	plt.ion()
	plt.show()
# plt.savefig(profile_shapefile.split('.')[0] + '_lp.png')

# num_pts=len(profile_data['sample_pts_PS'])
# values = range(len(pt_indx_to_plot))
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=np.min(pt_indx_to_plot), vmax=np.max(pt_indx_to_plot))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
print scalarMap.get_clim()
# cstr='rbcmrygbcmrygbcmrygbcmrygbcrbcmrygbcmrygbcmrygbcmrygbcrbcmrygbcmrygbcmrygbcmrygbcrbcmrygbcmrygbcmrygbcmrygbcrbcmrygbcmrygbcmrygbcmrygbcrbcmrygbcmrygbcmrygbcmrygbcrbcmrygbcmrygbcmrygbcmrygbcrbcmrygbcmrygbcmrygbcmrygbc'
# cm=plt.cm.get_cmap('jet')
# fig1=plt.figure(figsize=[20,16])
ax1=fig4.add_subplot('132',title='point speeds over time (color is position along profile)',ylabel='speed (m/d)',xlabel='year')
ax1.set_autoscaley_on(False)
ax1.set_ylim([speed_min, speed_max])
ax1.set_autoscalex_on(False)
ax1.set_xlim([year_min, year_max])
ax1.grid(b=None,axis='x',which='major',color='k',linestyle='-', alpha=0.75)
# color_min=0
# color_max=num_pts
# for indx in range(0,num_pts,skippnts_timeprofiles):
point_indicies_subsampled=np.array(pt_indx_to_plot)[range(0,len(pt_indx_to_plot),int(np.floor(len(pt_indx_to_plot)/plt_numpts_timeprofiles)))]
for indx in point_indicies_subsampled:
	point = profile_data['pts']['pt_speed'][indx]
# 	colorVal = scalarMap.to_rgba(values[indx])
	colorVal = scalarMap.to_rgba(indx)
	colorText = (
        '%6.0f m'%(sample_pts_frontdist[indx])
        )
	for idx2,pt in enumerate(point):
		if idx2==0:
			ax1.plot([pts['start_dec_year'][idx2], pts['stop_dec_year'][idx2]],[pt,pt],linewidth=2.0,color=colorVal,label=colorText)
		else:
			ax1.plot([pts['start_dec_year'][idx2], pts['stop_dec_year'][idx2]],[pt,pt],linewidth=2.0,color=colorVal)
# cb=plt.colorbar(sc,format='%6.1f')
x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax1.xaxis.set_major_formatter(x_formatter)
handles,labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='upper left')
for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(17)
# ax1.title(profile_shapefile.split('.')[0])
if not(args.nofigs):
	plt.show()
# plt.savefig(profile_shapefile.split('.')[0] + '.png')

# creates the plot with filled polygons - no gaps show up when plot resized
#
# maxspd=speed_max
# num_pts=len(profile_data['sample_pts_PS'])
# values = range(maxspd)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=speed_max)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
print scalarMap.get_clim()
ax2=fig4.add_subplot('133',title='speed evolution (color is speed (m/d))',xlabel='position along profile (m)',ylabel='year')
ax2.set_autoscalex_on(False)
ax2.set_xlim([dist_min, dist_max])
ax2.set_autoscaley_on(False)
ax2.set_ylim([year_min, year_max])
ax2.invert_yaxis()
ax2.grid(b=None,axis='y',which='major',color='k',linestyle='-', alpha=0.75)
color_min=speed_min
color_max=speed_max
plot_points_front_dist=np.array(profile_data['sample_pts_frontdist'])[pt_indx_to_plot]
plot_points_speed=np.array(profile_data['pts']['pt_speed'])[pt_indx_to_plot]
plot_points_x=np.array(profile_data['sample_pts_frontdist'])[pt_indx_to_plot]
plotted_points=[]
for indx in range(len(pt_indx_to_plot)):
	point_speeds = plot_points_speed[indx]
	pntx=plot_points_x[indx]
	if indx==0:
		halfdistpos=(plot_points_front_dist[indx+1] - plot_points_front_dist[indx])/2.0
		halfdistneg=-halfdistpos
	elif indx==(len(pt_indx_to_plot)-1):
		halfdistneg=(plot_points_front_dist[indx-1] - plot_points_front_dist[indx])/2.0
		halfdistpos=halfdistneg
	else:
		halfdistpos=(plot_points_front_dist[indx+1] - plot_points_front_dist[indx])/2.0
		halfdistneg=(plot_points_front_dist[indx-1] - plot_points_front_dist[indx])/2.0
	for idx2,epoch_speed in enumerate(point_speeds):
		cind=epoch_speed
		if cind>speed_max:
			cind=speed_max
		colorVal = scalarMap.to_rgba(cind)
		if epoch_speed!=-1.0:
			xvec=[pntx+halfdistneg,pntx+halfdistpos,pntx+halfdistpos,pntx+halfdistneg,pntx+halfdistneg]
			yvec=[pts['start_dec_year'][idx2],pts['start_dec_year'][idx2],pts['stop_dec_year'][idx2],pts['stop_dec_year'][idx2],pts['start_dec_year'][idx2]]
			plotted_points.append((pts['delts_days'][idx2],xvec,yvec,colorVal))
# for indx in pt_indx_to_plot:
# 	point_speeds = profile_data['pts']['pt_speed'][indx]
# 	pntx=profile_data['sample_pts_frontdist'][indx]
# 	if indx==0:
# 		halfdistpos=(profile_data['sample_pts_frontdist'][indx+1] - profile_data['sample_pts_frontdist'][indx])/2.0
# 		halfdistneg=-halfdistpos
# 	elif indx==(num_pts-1):
# 		halfdistneg=(profile_data['sample_pts_frontdist'][indx-1] - profile_data['sample_pts_frontdist'][indx])/2.0
# 		halfdistpos=halfdistneg
# 	else:
# 		halfdistpos=(profile_data['sample_pts_frontdist'][indx+1] - profile_data['sample_pts_frontdist'][indx])/2.0
# 		halfdistneg=(profile_data['sample_pts_frontdist'][indx-1] - profile_data['sample_pts_frontdist'][indx])/2.0
# 	for idx2,epoch_speed in enumerate(point_speeds):
# 		cind=epoch_speed
# 		if cind>speed_max:
# 			cind=speed_max
# 		colorVal = scalarMap.to_rgba(cind)
# 		if epoch_speed!=-1.0:
# 			xvec=[pntx+halfdistneg,pntx+halfdistpos,pntx+halfdistpos,pntx+halfdistneg,pntx+halfdistneg]
# 			yvec=[pts['start_dec_year'][idx2],pts['start_dec_year'][idx2],pts['stop_dec_year'][idx2],pts['stop_dec_year'][idx2],pts['start_dec_year'][idx2]]
# 			plotted_points.append((pts['delts_days'][idx2],xvec,yvec,colorVal))
sorted_points=sorted(plotted_points,key=itemgetter(0),reverse=True)
for point in sorted_points:
	ax2.fill(point[1],point[2],color=point[3])
scalarMap._A=[]   # need to set this array even though not used in order to get colormap function to work
cb2=plt.colorbar(scalarMap)
y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax2.yaxis.set_major_formatter(y_formatter)
for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
             ax2.get_xticklabels() + ax2.get_yticklabels() + cb2.ax.get_yticklabels()):
    item.set_fontsize(17)

fig4.subplots_adjust(left=0.06,right=0.97,top=0.93, bottom=0.07, wspace=0.25, hspace=0.15)

figtitlestr=profile_shapefile.split('.')[0]
if(args.Line_ID_value_to_plot):
	figtitlestr=figtitlestr + ' - Profile ' + args.Line_ID_value_to_plot
figtitle=fig4.suptitle(figtitlestr)
figtitle.set_fontsize(20)

graphic_outfile='%s_profiles_sampled_%s.png'%(profile_shapefile.split('.')[0],time.strftime("%Y-%m-%d"))
plt.savefig(graphic_outfile)
