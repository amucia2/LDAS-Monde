# Figures for paper LDAS-Monde EnKF for HESS
# author: B. Bonan
# Last modified: March 2021

## patchFrac_figures_Tony.py
import sys, glob, os, re, math, time
os.chdir('/cnrm/vegeo/muciaa/ldas_chain_python/LDAS_v1.3.0')
#os.chdir('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr')
from ldasPost import *
from ldasMapSet import *
os.chdir('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing')

options_file = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/options_US00.py'

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

# Load options files
initPost(options_file)

##################################################################################
##
## Create figures for LAI 

# Read all data related to ekf_lai for LAI
obs = 'LAI_V2'
mod_ana_dir['Analysis'] = '/cnrm/vegeo/muciaa/NO_SAVE/US00/ekf_lai/'
out_all_dir['Analysis'] = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/ekf_lai/'
### Originally data_SEKF
data_LAI = readAllData(obs,analysis_period, assim_hour, patch_frac, out_all_dir, \
                        post_from_pickle, mod_ana_dir, ['Model','Analysis','Obs'], \
                        IPNT=IPNT, map_set=map_set, filter_nan=not('SIF' in obs or 'E_ACTUAL' in obs or 'FLGPP' in obs))

# Read all data related to ekf_vod for LAI
obs = 'LAI_V2'
mod_ana_dir['Analysis'] = '/cnrm/vegeo/muciaa/NO_SAVE/US00/ekf_vod/'
out_all_dir['Analysis'] = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/ekf_vod/'
### Originially data_EnKF
data_VOD = readAllData(obs,analysis_period, assim_hour, patch_frac, out_all_dir, \
                        post_from_pickle, mod_ana_dir, ['Model','Analysis','Obs'], \
                        IPNT=IPNT, map_set=map_set, filter_nan=not('SIF' in obs or 'E_ACTUAL' in obs or 'FLGPP' in obs))

###################################################################################
##
## Displaying time series of LAI for patches (with patch fraction > 0.5)

# Reading patch fractions
PREPdir = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/sfx-trip/pgd_prep/'
file_prep = 'PREP.nc'
PREP_set = Dataset(PREPdir+file_prep,'r')
Auxi = PREP_set.variables['PATCH'][:]
patch_frac = Auxi.data
patch_frac[patch_frac == Auxi.fill_value] = np.nan
PREP_set.close()

data_reord_SEKF = data_LAI.reindex(minor=pl.arange(map_set.X.size),copy=False)
data_reord_EnKF = data_VOD.reindex(minor=pl.arange(map_set.X.size),copy=False)

#data_reord_SEKF = data_SEKF.reindex('point',pl.arange(map_set.X.size))
#data_reord_EnKF = data_EnKF.reindex(pl.arange(map_set.X.size))

#data_reord_SEKF = data_SEKF
#data_reord_EnKF = data_EnKF

patch_num = [3,4,6,9]
patch_name = ['Deciduous','Coniferous','c3_crops','Grasslands']
indx = 0



for ip in patch_num:
    
    patch_test = patch_frac[ip,:,:].reshape((280*140))
    
    fig = pl.figure(figsize=[9/3.*4.5,7/6.*4.5])

    #t = data_SEKF.major_axis.to_datetime()
    #t = data_SEKF['time']
    # = pd.date_range(['2000-01-01', '2018-12-31'])
    t = pd.date_range('2000-01-01', freq='10D',periods=684)
    xlim = [t[0]-pd.DateOffset(months=1), t[-1]+pd.DateOffset(months=1)]
    xtickpos = pd.date_range(t[0]-pd.DateOffset(months=1),t[-1]+pd.DateOffset(months=1),freq='A')+pd.DateOffset(days=1)
    xticklbl = [xtickpos[i].strftime('%Y') for i in range(len(xtickpos))]
    
    ax = fig.add_subplot(1,1,1)
    ph = []
    p_tmp = ax.plot(t, data_reord_SEKF.loc['Model'].loc[:,patch_test > 0.5].mean(axis=1), 'b', label='Model')
    ph.append(p_tmp)
    p_tmp = ax.plot(t, data_reord_SEKF.loc['Analysis'].loc[:,patch_test > 0.5].mean(axis=1), '--', color='purple', linewidth=2.0, label='EKF LAI')
    ph.append(p_tmp)
    p_tmp = ax.plot(t, data_reord_EnKF.loc['Analysis'].loc[:,patch_test > 0.5].mean(axis=1), 'r', label='EKF VOD')
    ph.append(p_tmp)
    p_tmp = ax.plot(t, data_reord_SEKF.loc['Obs'].loc[:,patch_test > 0.5].mean(axis=1), 'g.:', label='Obs')
    ph.append(p_tmp)

    ax.set_xlim(xlim)
    ax.set_xticks(xtickpos)
    ax.set_xticklabels(xticklbl)
    ax.set_ylabel('LAI [$m^2.m^{-2}$]')

    labs = [ph[i][0].get_label() for i in range(0,4)]
    fig.legend(ax.lines,labs,ncol=4,loc=9)

    fig_name_tmp = 'LAI_daily_'+patch_name[indx]+'_time_series.png'
    fig.savefig(fig_name_tmp, dpi=DEF_FIG_DPI, bbox_inches='tight')
    pl.close()
    
    indx = indx+1

sys.exit()

patch_num = 3

patch_test = patch_frac[ip,:,:].reshape((280*140))

fig = pl.figure(figsize=[9/3.*4.5,7/6.*4.5])

t = data_SEKF['time']
t = pd.date_range('2000-01-01', freq='D',periods=684)
xlim = [t[0]-pd.DateOffset(months=1), t[-1]+pd.DateOffset(months=1)]
xtickpos = pd.date_range(t[0]-pd.DateOffset(months=1),t[-1]+pd.DateOffset(months=1),freq='A')+pd.DateOffset(days=1)
xticklbl = [xtickpos[i].strftime('%Y') for i in range(len(xtickpos))]

ax = fig.add_subplot(1,1,1)
ph = []
p_tmp = ax.plot(t, data_reord_SEKF.loc['Model'].loc[:,patch_test > 0.5].mean(axis=1), 'b', label='Model')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_SEKF.loc['Analysis'].loc[:,patch_test > 0.5].mean(axis=1), '--', color='purple', linewidth=2.0, label='EKF LAI')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_EnKF.loc['Analysis'].loc[:,patch_test > 0.5].mean(axis=1), 'r', label='EKF VOD')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_SEKF.loc['Obs'].loc[:,patch_test > 0.5].mean(axis=1), 'g.:', label='Obs')
ph.append(p_tmp)

ax.set_xlim(xlim)
ax.set_xticks(xtickpos)
ax.set_xticklabels(xticklbl)
ax.set_ylabel('LAI [$m^2.m^{-2}$]')

labs = [ph[i][0].get_label() for i in range(0,4)]
fig.legend(ax.lines,labs,ncol=4,loc=9)

fig_name_tmp = 'LAI_daily_'+patch_name[indx]+'_time_series.png'
fig.savefig(fig_name_tmp, dpi=DEF_FIG_DPI, bbox_inches='tight')
pl.close()

indx = indx+1
###################################################################################
##
## Displaying RSMD for pixels with patch fraction > 0.5

## 280x140

# Compute scores for LAI
scores_SEKF = compute_scores_per_point(data_SEKF, 'season')
scores_EnKF = compute_scores_per_point(data_EnKF, 'season')

data_reord_SEKF = scores_SEKF.reindex(minor=pl.arange(map_set.X.size),copy=False)
data_reord_EnKF = scores_EnKF.reindex(minor=pl.arange(map_set.X.size),copy=False)

# Reading patch fractions
PREPdir = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/sfx-trip/pgd_prep/'
file_prep = 'PREP.nc'
PREP_set = Dataset(PREPdir+file_prep,'r')
Auxi = PREP_set.variables['PATCH'][:]
patch_frac = Auxi.data
patch_frac[patch_frac == Auxi.fill_value] = np.nan
PREP_set.close()

t = range(0, scores_SEKF.major_axis.size)
xlim = [-.25, scores_SEKF.major_axis.size-.75]
xtickpos = range(0,11,2)
xticklbl = [m[0:3] for m in scores_SEKF.major_axis[range(0,11,2)]]

grid_spec, fig, title_rel_h, cbar_aspect = map_set.figureLayout(title=None,n_rows=2,n_cols=3)

ax = pl.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)

# RMSD or R for LAI for all the grid

#sco = 'RMSD'
sco = 'Correlation'

ph = []
p_tmp = ax.plot(t, data_reord_SEKF[sco]['Model'].mean(axis=1).values, 'b', label='Model')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_SEKF[sco]['Analysis'].mean(axis=1).values, '--', color='purple', linewidth=2.0, label='EKF LAI')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_EnKF[sco]['Analysis'].mean(axis=1).values, 'r', label='EKF VOD')
ph.append(p_tmp)

ax.set_title('R for LAI')
ax.set_xlim(xlim)
ax.set_ylim([0,1])
ax.set_ylabel('[$m^2.m^{-2}$]')
ax.set_xticks(xtickpos)
ax.set_xticklabels(xticklbl)
ax.grid()

pl.tight_layout()

ax = pl.subplot2grid((2,6), (0,2), colspan=2)

##patch_test = patch_frac[4,:,:].reshape((202*296))
patch_test = patch_frac[4,:,:].reshape((280*140))
ph = []
p_tmp = ax.plot(t, data_reord_SEKF[sco]['Model'].loc[:,patch_test > 0.5].mean(axis=1).values, 'b', label='Model')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_SEKF[sco]['Analysis'].loc[:,patch_test > 0.5].mean(axis=1).values, '--', color='purple', linewidth=2.0, label='EKF LAI')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_EnKF[sco]['Analysis'].loc[:,patch_test > 0.5].mean(axis=1).values, 'r', label='EKF VOD')
ph.append(p_tmp)

ax.set_title('R (Coniferous)')
ax.set_xlim(xlim)
ax.set_ylim([0,1])
ax.set_xticks(xtickpos)
ax.set_xticklabels(xticklbl)
ax.grid()

pl.tight_layout()

ax = pl.subplot2grid((2,6), (0,4), colspan=2)

patch_test = patch_frac[6,:,:].reshape((280*140))
ph = []
p_tmp = ax.plot(t, data_reord_SEKF[sco]['Model'].loc[:,patch_test > 0.5].mean(axis=1).values, 'b', label='Model')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_SEKF[sco]['Analysis'].loc[:,patch_test > 0.5].mean(axis=1).values, '--', color='purple', linewidth=2.0, label='EKF LAI')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_EnKF[sco]['Analysis'].loc[:,patch_test > 0.5].mean(axis=1).values, 'r', label='EKF VOD')
ph.append(p_tmp)

ax.set_title('R (c3 crops)')
ax.set_xlim(xlim)
ax.set_ylim([0,1])
ax.set_xticks(xtickpos)
ax.set_xticklabels(xticklbl)
ax.grid()

pl.tight_layout()

ax = pl.subplot2grid((2,6), (1,0), colspan=2)

patch_test = patch_frac[9,:,:].reshape((280*140))
ph = []
p_tmp = ax.plot(t, data_reord_SEKF[sco]['Model'].loc[:,patch_test > 0.5].mean(axis=1).values, 'b', label='Model')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_SEKF[sco]['Analysis'].loc[:,patch_test > 0.5].mean(axis=1).values, '--', color='purple', linewidth=2.0, label='EKF LAI')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_EnKF[sco]['Analysis'].loc[:,patch_test > 0.5].mean(axis=1).values, 'r', label='EKF VOD')
ph.append(p_tmp)

ax.set_title('R (Grasslands)')
ax.set_xlim(xlim)
ax.set_ylim([0,1])
ax.set_ylabel('[$m^2.m^{-2}$]')
ax.set_xticks(xtickpos)
ax.set_xticklabels(xticklbl)
ax.grid()

pl.tight_layout()

ax = pl.subplot2grid((2,6), (1,4), colspan=2)

patch_test = patch_frac[3,:,:].reshape((280*140))

ph = []
p_tmp = ax.plot(t, data_reord_SEKF[sco]['Model'].loc[:,patch_test > 0.5].mean(axis=1).values, 'b', label='Model')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_SEKF[sco]['Analysis'].loc[:,patch_test > 0.5].mean(axis=1).values, '--', color='purple', linewidth=2.0, label='EKF LAI')
ph.append(p_tmp)
p_tmp = ax.plot(t, data_reord_EnKF[sco]['Analysis'].loc[:,patch_test > 0.5].mean(axis=1).values, 'r', label='EKF VOD')
ph.append(p_tmp)

ax.set_title('R (Deciduous)')
ax.set_xlim(xlim)
ax.set_ylim([0,1])
ax.set_ylabel('[$m^2.m^{-2}$]')
ax.set_xticks(xtickpos)
ax.set_xticklabels(xticklbl)
ax.grid()

pl.tight_layout()

labs = [ph[i][0].get_label() for i in range(0,3)]
fig.legend(ax.lines, labs, ncol=1, loc=4, bbox_to_anchor=(.55, 0.17), borderaxespad=0.)

pl.tight_layout()

fig_name_tmp = 'CGLS_LAI_R.png'
fig.savefig(fig_name_tmp, dpi=DEF_FIG_DPI, bbox_inches='tight')
pl.close()

