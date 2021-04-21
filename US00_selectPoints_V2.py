import sys, glob, os, re, math, time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
import pylab as pl
import scipy.optimize as optim
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable

pl.ioff()
pl.rc('mathtext',default='regular') # use normal font for math expressions in rendered text
MONTH_SHORT_NAMES = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
# Dictionary mapping observation names to model output names
OBS_MOD = {'FAPAR':'FAPAR', 'FAPAR_V1':'DFAPARC', 'FAPAR_V2':'FAPAR', 'FAPAR300':'FAPAR',
           'LAI':'LAI', 'LAI_V1':'LAI', 'LAI_V2':'LAI', 'LAI300':'LAI',
           'LST':'TS', 'LST-6':'TS', 'LST-12':'TS',
           'SA':'TALB','raw_SWI_nc':'WG2','raw_SWI':'WG2',
           'SA':'TALB',
           'SWI':'WG2', 'SWI_nc':'WG2', 'SSM':'WG2', 'SSM_COMBINED_v4':'WG2', 'SSM_COMBINED_v4.5':'WG2', 'SSM_ACTIVE_v4':'WG2', 'SSM_PASSIVE_v4':'WG2',
           'SIF':'GPP', 'SIF-L3':'GPP', 'SIF-F':'GPP', 'FLCGPP':'GPPC',
           'SIF-L3_V27':'GPPC','SIF-L3_V27_d':'GPPC','SIF-L3_V27_p':'GPPC',
           'E_ACTUAL':'EVAPC', 'E_t':'LETRC', 'E_b':'LEGC', 'E_i':'LERC', 'FLLE':'EVAPC', 'FLGPP':'GPPC',
           'SM_QD':'WG2', 'SM_QD_2540':'WG2', 'SM_QD_1030':'WG2','VODC':'LAI','VODX':'LAI','VODKu':'LAI',
           'VOD_QD_1030':'LAI', 'VOD_QD_2540':'LAI', 'VOD_L_ASC':'LAI','SIGMA_40_QD':'LAI','IMS':'PSNG','NCALDAS':'EVAPC','FLDAS':'EVAPC','NCALDAS_LAI':'EVAPC','SCFV-MODIS':'PSNG',
           'ET_04':'EVAPC','ET_05':'EVAPC'}
           #'VOD_QD_1030':'GPPC', 'VOD_QD_2540':'GPPC', 'VOD_L_ASC':'GPPC','SIGMA_40_QD':'LAI'}
OBS_MOD_sameUnit = {'FAPAR':True, 'FAPAR_V1':True, 'FAPAR_V2':True, 'FAPAR300':True,
           'LAI':True, 'LAI_V1':True, 'LAI_V2':True, 'LAI300':True,
           'LST':True, 'LST-6':True, 'LST-12':True,
           'SA':True,'VODC': False,'VODX': False,'VODKu': False,
           'VOD_QD_1030': False, 'VOD_QD_2540':False, 'VOD_L_ASC':False,
           'SWI':True, 'SWI_nc':True,'SSM':True, 'SSM_COMBINED_v4':True, 'SSM_COMBINED_v4.5':True, 'SSM_ACTIVE_v4':True, 'SSM_PASSIVE_v4':True,
           'SM_QD':True,'SM_QD_2540':True,'SM_QD_1030':True,'SIGMA_40_QD':False,
           'SIF':False, 'SIF-L3':False, 'SIF-F':False, 'FLCGPP':True,
           'SIF-L3_V27':False,'SIF-L3_V27_d':False,'SIF-L3_V27_p':False,
           'E_ACTUAL':True, 'E_t':True, 'E_b':True, 'E_i':True, 'FLLE':True, 'FLGPP':True,
           'ET_04':True, 'ET_05':True}


numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts

def check_post_zones(post_zones):
    '''
    Function used to check and format the post_zones dict.
    '''
    tmp_zones = dict()
    for zname,zone in post_zones.items():
        if zname=='def':
            tmp_zones['model zone'] = None
        elif zname=='point':
            for i,xy in zip(range(len(zone)),zone):
                pnt = 'pnt_{:02d}'.format(i+1)
                tmp_zones[pnt] = [[xy[0],xy[1]],[xy[0],xy[1]]]
        elif isinstance(zone,str): # If a netcdf file is given
            if not os.path.isfile(zone):
                raise Exception('\nNetcdf file given for mask does not exist: {0}'.format(zone))
            with Dataset(zone) as nc:
                if zname=='all':
                    for name in nc.variables['Name']: tmp_zones[name] = zone
                elif zname not in nc.variables['Name']:
                    raise Exception('Mask name {0} not in {1}'.format(zname,zone))
                else: tmp_zones[zname] = zone
        else:
            tmp_zones[zname] = zone
    return tmp_zones

def prepareZone(zname,zone,map_set,graphics_dir):
    '''
    Prepares variables for the selected zone
    '''
    if zname=='model zone': # By default, zone = model grid
        map_set.zone = None
        graphics_dir_zone = graphics_dir
        IPNT = None
    else:
        graphics_dir_zone = graphics_dir+zname+'/'
        if isinstance(zone,str): # If a netcdf file is given
            nc = Dataset(zone)
            zone_id = nc.variables['Id'][pl.argwhere(nc.variables['Name'][:]==zname)[0][0]]
            IPNT = pl.where(nc.variables['Mask'][::-1,:]==zone_id,True,False).ravel()
            map_set.zone = [[min(map_set.X.ravel()[IPNT]),max(map_set.Y.ravel()[IPNT])],
                            [max(map_set.X.ravel()[IPNT]),min(map_set.Y.ravel()[IPNT])]]
            nc.close()
        else:                    # If zone is a rectangular area [[lon_min,lat_max],[lon_max,lat_min]]
            map_set.zone = zone
            IPNT = pl.logical_and(map_set.mod_lon.ravel()>=zone[0][0],map_set.mod_lon.ravel()<=zone[1][0])
            IPNT = pl.logical_and(IPNT,\
                   pl.logical_and(map_set.mod_lat.ravel()<=zone[0][1],map_set.mod_lat.ravel()>=zone[1][1]))
        if map_set.mod_grid_type!='uniform': IPNT = IPNT[map_set.remap]
    # Directories for the selected zone
    graphics_dir_zone = graphics_dir_zone.replace(' ','_')
    makeDirIfNeeded(graphics_dir_zone)
    return map_set,graphics_dir_zone,IPNT

from ldasModule import *
from ldasMapSet import mapSet, uniformTicks, Basemap, Figure


def postObs(obs, data=None, minim=None, maxim=None):
    importLocals(currentframe())
    print('===================================================')
    print('Post-processing for observation: {0}'.format(obs)+(', model variable: {1}'.format(obs, OBS_MOD[obs]) if obs in OBS_MOD.keys() else ''))
    # Read all data
    if data is None:
        data = readAllData(obs, analysis_period, assim_hour, patch_frac, out_all_dir, \
                            post_from_pickle, mod_ana_dir, ['Model','Analysis','Obs'], \
                            IPNT=IPNT, map_set=map_set, filter_nan=not('SIF' in obs or 'E_ACTUAL' in obs or 'FLGPP' in obs))
    if 0 in data.shape:
        print('No dataset or no time step or no point. Continuing with next obs...')
        return
    items = data.coords['name'].values
    col_labels = [modexp_names[n] if n in modexp_names.keys() else n for n in items]
    if 'Analysis' in data.coords['name'].values and 'Model' in data.coords['name'].values:
        col_labels.append(col_labels[pl.find(items=='Analysis')[0]]+'-'+col_labels[pl.find(items=='Model')[0]])
    for time_window,t in zip(list_time_windows,range(len(list_time_windows))):
        if any([var_spatial_maps_out[t],var_time_series_out[t]]):
            print('Time window applied for {0}: {1}'.format(obs, time_window))
            # Apply time window
            data_tw = applyTimeWindow(data, time_window)
            # Additional diagnostis: analysis-model
            items_tw = data_tw.coords['name'].values
            if 'Analysis' in items_tw and 'Model' in items_tw:
                diff = data_tw.sel(name='Analysis')-data_tw.sel(name='Model')
                diff.values[diff.values==0] = np.nan
                diff = diff.assign_coords(name='Analysis-Model')
                diff = diff.expand_dims('name')
                data_tw = xr.concat([data_tw, diff], dim='name')
            #
            # Draw maps of variables: model, obs, analysis, diff(analysis-model)
            if minim is None and obs in CBAR_MINMAX: minim = CBAR_MINMAX[obs][0]
            if maxim is None and obs in CBAR_MINMAX: maxim = CBAR_MINMAX[obs][1]
            if var_spatial_maps_out[t]: draw_map_vars(map_set,obs,data_tw,graphics_dir=graphics_dir_zone,time_window=time_window,fig_dpi=fig_dpi,minim=minim,maxim=maxim,col_labels=col_labels)
            # Draw time series of variables (spatial average over the whole domain)
            if var_time_series_out[t]: draw_series_vars(data_tw.mean('point'), obs, graphics_dir=graphics_dir_zone, time_window=time_window, fig_dpi=fig_dpi, col_labels=col_labels)
            #
        # Compute statistics ('R','bias','std','rmse') and show
        if 'Obs' in items and any([stats_spatial_maps_out[t],stats_time_series_out[t],stats_table_out[t]]):
            # Draw maps of statistic scores: model vs obs, analysis vs obs, diff(analysis-model)
            if stats_spatial_maps_out[t]:
                if 'SIF' in obs or 'E_ACTUAL' in obs or 'FLGPP' in obs:
                    scores_per_point = compute_scores_per_point(applyTimeWindow(data,'month'), time_window)
                else:
                    scores_per_point = compute_scores_per_point(data, time_window)
                draw_map_stats(map_set, scores_per_point, graphics_dir=graphics_dir_zone, obs_var=obs, time_window=time_window, fig_dpi=fig_dpi, col_labels=col_labels)
            # Draw timeseries of statistic scores (spatial average over the whole domain)
            if 'SIF' in obs: scores_for_all_points = compute_scores_for_all_points(applyTimeWindow(data,'month'), time_window)
            ### Here, I compute both scores_per_point, and scores for all points at the same time - Take out if not needed, as this slows the process
            else           : scores_for_all_points = compute_scores_for_all_points(data, time_window); scores_per_point = compute_scores_per_point(data, time_window)
            if stats_time_series_out[t]: draw_series_stats(scores_for_all_points, obs, graphics_dir=graphics_dir_zone, time_window=time_window, fig_dpi=fig_dpi, col_labels=col_labels)
            if stats_table_out[t]: drawScoresTable(scores_for_all_points, obs=obs, graphics_dir=graphics_dir_zone, time_window=time_window, fig_dpi=fig_dpi)
            #if stats_table_out[t]: drawScoresTable(scores_per_point.mean(axis=3), obs=obs, graphics_dir=graphics_dir_zone, time_window=time_window, fig_dpi=fig_dpi)
            #pdb.set_trace()

            #scores_for_all_points.to_netcdf('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/ldasPostScores/V3/{0}_{1}_{2}_scores_for_all_points.nc'.format(zname,obs,fname))
            #scores_per_point.to_netcdf('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/ldasPostScores/V3/{0}_{1}_{2}_scores_per_point.nc'.format(zname,obs,fname))

            #return scores_for_all_points
            #return scores_per_point



options_file = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/options_US00.py'

#initPost(options_file)
if options_file is not None:
    local_vars = dict()
    #try: execfile(options_file,local_vars)
    try: exec(open(options_file).read(),local_vars)
    except: raise Exception('\nError in the options file!')
else:
    local_vars = currentframe(1).f_globals
# Update local variables and check options
#for key,value in local_vars.iteritems(): exec(key+"=value")
for key,value in local_vars.items(): exec(key+"=value")
checkOptions(locals(),options_file)
#
print('\nInitializing LDAS post-processing ({0} - {1}).\n'.format(analysis_period[0],analysis_period[1]))
#
# Make initializations
# If LST is in the obs list, then LST is analyzed at 6am and 12pm
if 'LST' in obs_names:
    id_LST = obs_names.index('LST')
    del to_assim[id_LST]
    del rescale_calibs[id_LST]
    del rescale_calib_periods[id_LST]
    obs_names = filter(lambda x: x != 'LST',obs_names) + ['LST-6', 'LST-12']
    to_assim  = to_assim + [False, False]
    rescale_calibs = rescale_calibs + [None,None]
    rescale_calib_periods = rescale_calib_periods+[[],[]]
#
# Main parameters definition
out_all_dir = {'Obs': out_obs_dir, 'Model': out_mod_dir, 'Analysis': out_ana_dir}
post_from_pickle = {'Obs': post_obs_from_pickle, 'rawObs': post_obs_from_pickle, \
                    'Model': post_mod_from_pickle, 'Analysis': post_ana_from_pickle}
obs_assim   = [obs_names[i] for i in range(len(obs_names)) if to_assim[i]==True ]
obs_NOassim = [obs_names[i] for i in range(len(obs_names)) if to_assim[i]==False]
obs_assim_rescaled = [obs_names[i] for i in range(len(obs_names)) if to_assim[i] and rescale_calibs[i]]
patch_frac = getModelPatchFractions(mod_pgd_path)
mod_grid,trip_grid = parseOptionsNam(options_path)
mod_ana_dir = {'Model': openloop_dir, 'Analysis': analysis_dir}
if modexp_names=='def': modexp_names = {'Model':'Model', 'Analysis':'Analysis'}
else: modexp_names = {'Model':modexp_names[0], 'Analysis':modexp_names[1]}
patch_post = pl.array(patch_out)-1
#
# Create instance of mapSet
map_set = mapSet(mod_grid,pgd=mod_pgd_path)
if len(trip_vars)>0:
    map_set_trip = mapSet(trip_grid)
#local_
# Check for netcdf files given as mask
post_zones = check_post_zones(post_zones)
#
# Initialize loop over zones
#zname,zone = post_zones.items().next()
zname,zone = next(iter(post_zones.items()))
map_set,graphics_dir_zone,IPNT = prepareZone(zname,zone,map_set,graphics_dir)
#
# Update caller local variables
#local_vars = currentframe(1).f_globals
local_vars = currentframe().f_globals
for key,value in locals().items(): local_vars[key] = value

def get_zone_grid(mod_grid, zone):
    zone_grid = dict()
    zone_grid['type'] = mod_grid['type']
    zone_grid['res'] = mod_grid['res']
    zone_grid['zone'] = zone

    return zone_grid

### Chop up Zones
if zone is not None:
    map_set = mapSet(mod_grid, pgd=mod_pgd_path)
    map_set,graphics_dir_zone,IPNT = prepareZone(zname,zone,map_set,graphics_dir)

    zone_grid = get_zone_grid(mod_grid, zone)
    map_set = mapSet(zone_grid, pgd=mod_pgd_path)
    map_set,graphics_dir_zone,_ = prepareZone(zname,zone,map_set,graphics_dir)
else:
    map_set,graphics_dir_zone,IPNT = prepareZone(zname,zone,map_set,graphics_dir)



### Load Data
vod = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/LAIfromVODCAX_V8_2003-2018.PData')
vodx_int = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/LAIfromVODCAX_V7_2003-2018_interpolated.PData')
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/sfx-trip/LAI_V2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),pd.read_pickle(file[18])]
obs = pd.concat(df)

date = pd.date_range(start="2003-01-01 09:00:00",end='2018-12-31 09:00:00',freq='D')
obs = obs.reindex(date,fill_value=np.nan)
lai_cgls = lai_cgls.interpolate(inplace=False,limit_area='inside')

###################################################################
### Data Processing 
### Only takes data from the same time as observations 
###################################################################

vod[np.invert(~np.isnan(obs))]=np.nan
obs[np.invert(~np.isnan(vod))]=np.nan


vodx_r = vod.corrwith(obs,axis=0)
vodx_int_r = vodx_int.corrwith(obs,axis=0)
###################################################################
### Graphing Function(s)
###################################################################

vodm = vod.mean(axis=1)
obsm = obs.mean(axis=1)
vodx_intm = vodx_int.mean(axis=1)

vod2 = vodm['2012-01-01':'2012-12-31']
obs2 = obsm['2012-01-01':'2012-12-31']
vodx_int2 = vodx_intm['2012-01-01':'2012-12-31']

plt.title('LAI and Matched LAI from VODX : 2003-2018')
#plt.plot(obsm,color='black',label='CGLS LAI Obs',marker='.',markersize=12)
plt.plot(obsm,color='black',label='CGLS LAI Obs')
plt.plot(vodm,color='green',label='LAI from VODX')
plt.plot(vodx_intm,color='red',label='LAI from VODX - Interpolated Obs')

plt.legend()
plt.show()

plt.title('LAI and Matched LAI from VODX : 2012')
#plt.plot(obs2,color='black',label='CGLS LAI Obs',marker='.',markersize=12)
plt.plot(obs2,color='black',label='CGLS LAI Obs')
plt.plot(vod2,color='green',label='LAI from VODX')
plt.plot(vodx_int2,color='red',label='LAI from VODX - Interpolated Obs')

plt.legend()
plt.show()

###################################################################
### Mapping Function(s)
###################################################################

tmp = obs.copy()*np.nan
t = tmp.mean(axis=0)

for i in range(len(t)):
    print(i)
    t.loc[i] = i

v1 = t.values.reshape((280,140))

#fig, axes = plt.subplots(1)

map = Basemap(llcrnrlon=-130, llcrnrlat=20, urcrnrlon=-60, urcrnrlat=55,lat_0=40,lon_0=-98,resolution='l')
#map.shadedrelief()
#cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.5)
map.drawcountries(1.5)
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
cs1 = map.imshow(v1,interpolation='none')
#cs1 = map.imshow(v1,cmap='RdYlGn')

#cbar1 = map.colorbar(v1,location='bottom',pad="10%")
#cbar1.set_label("LAI [$m^2$/$m^2$]")
#tick_locator = ticker.MaxNLocator(nbins=4)
#cbar1.locator = tick_locator
#cbar1.update_ticks()

plt.show()


