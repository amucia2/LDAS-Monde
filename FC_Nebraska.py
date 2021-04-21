import sys, glob, os, re, math, time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
import pylab as pl
import scipy.optimize as optim
import scipy as sp
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import zscore
from scipy import stats


FILL = 999.
###################################################################
### Useful Functions
###################################################################
def readFromFile(file_name,var_names,row_slice=None,col_slice=None,row_invert=False):
    '''
    Read values from netcdf or hdf5 file
    Arguments:
        file_name (str): full path to the file to read
        var_names (list): names of variables to read
        row_slice/col_slice (slices): indices to slice the variables to the output region
        row_invert (bool): it True, row order is inverted
    Returns:
        obs_values (list): values read from file
    '''
    # determine format
    if file_name[-2:] not in ('h5','nc'):
        raise Exception('filename extension must be h5 of nc')
    file_format = file_name[-2:]
    if isinstance(var_names,basestring): var_names = [var_names]
    obs_values = []
    if file_format == 'nc':
        fid = Dataset(file_name)
        for i in range(len(var_names)):
            v = fid.variables[var_names[i]]
            tmp = v[:].squeeze()
            if row_invert: tmp = tmp[::-1,:]
            if row_slice is not None and col_slice is not None: tmp = tmp[row_slice,col_slice]
            tmp = tmp.astype(float)
            if not type(tmp) == pl.ma.core.MaskedArray: tmp = pl.ma.masked_array(tmp)
            for a in v.ncattrs():
                if a.lower() in ['missing_value','_fillvalue']: tmp.mask = pl.logical_or(tmp.mask,(tmp==float(v.getncattr(a))))
            tmp = tmp.filled(FILL)
            obs_values.append(tmp)
        fid.close()
    return obs_values


def pearsonr_ci(x,y,alpha=0.05):
    ''' calculate Pearson correlation along with the confidence interval using scipy and numpy
    Parameters
    ----------
    x, y : iterable object such as a list or np.array
      Input for correlation calculation
    alpha : float
      Significance level. 0.05 by default
    Returns
    -------
    r : float
      Pearson's correlation coefficient
    pval : float
      The corresponding p value
    lo, hi : float
      The lower and upper bound of confidence intervals
    '''
    r, p = stats.pearsonr(x,y)
    r_z = np.arctanh(r)
    se = 1/np.sqrt(x.size-3)
    z = stats.norm.ppf(1-alpha/2)
    lo_z, hi_z = r_z-z*se, r_z+z*se
    lo, hi = np.tanh((lo_z, hi_z))
    return r, p, lo, hi

def save(path, ext='png', close=True, verbose=True):
    """Save a figure from pyplot.
    Parameters
    ----------
    path : string
    The path (and filename, without the extension) to save the
    figure to.
    ext : string (default='png')
    The file extension. This must be supported by the active
    matplotlib backend (see matplotlib.backends module).  Most
    backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.
    close : boolean (default=True)
    Whether to close the figure after saving.  If you want to save
    the figure multiple times (e.g., to multiple formats), you
    should NOT close it in between saves or you will have to
    re-plot it.
    verbose : boolean (default=True)
    whether to print information about when and where the image
    has been saved.
    """
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
       directory = '.'
    #If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)
    # The final path to save to
    savepath = os.path.join(directory, filename)
    if verbose:
        print("Saving figure to '%s'..." % savepath),
    # Actually save the figure
    plt.savefig(savepath)
    # Close it
    if close:
        plt.close()
    if verbose:
        print("Done")

def ZscoreM(x):
    for ii in [1,2,3,4,5,6,7,8,9,10,11,12]:
        for name in x.columns:
            season = (x[name].index.month == ii)
            x[name].ix[season] = (x[name].ix[season]-x[name].ix[season].mean())/x[name].ix[season].std()
    return x

def ZscoreY(x):
    for ii in [1,2,3,4,5,6,7,8,9,10,11,12]:
        season = (x.index.month == ii)
        x.ix[season] = (x.ix[season]-x.ix[season].mean())/x.ix[season].std()
    return x

def getModelPatchFractions(mod_pgd_path):
    '''
    Read the patch fractions of model grid points from file.
    Arguments:
        mod_pgd_path (str): directory of PREP.nc
    Returns:
        mod_patch_frac (numpy.ndarray): patch fractions of each model grid point; shape = n_points,n_patches
    '''
    prep_file = mod_pgd_path+'PREP.nc'
    if not os.path.isfile(prep_file):
        raise Exception(prep_file+' does not exist.')
    frac = readFromFile(prep_file,['PATCH'])[0]
    return frac.reshape((frac.shape[0],-1))


###################################################################
### Import Data  
###################################################################


### Import individual yearly files
lai_obs = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/observations/operational/LAI_V2_2018-01-01_2018-12-31.PData')
swi_obs = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/observations/operational/SWI_nc_2018-01-01_2018-12-31.PData')


lai_ol = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/ol/LAI_2018-01-01_2018-12-31.PData')
swi_ol = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/ol/WG2_2018-01-01_2018-12-31.PData')
lai_ekf = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/ekf/LAI_2018-01-01_2018-12-31.PData')
swi_ekf = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/ekf/WG2_2018-01-01_2018-12-31.PData')
### -- -- -- -- -- -- -- 
lai_fc2 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_2/LAI_2018-01-01_2018-12-31.PData')
lai_fc3 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_3/LAI_2018-01-01_2018-12-31.PData')
lai_fc4 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_4/LAI_2018-01-01_2018-12-31.PData')
lai_fc5 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_5/LAI_2018-01-01_2018-12-31.PData')
lai_fc6 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_6/LAI_2018-01-01_2018-12-31.PData')
lai_fc7 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_7/LAI_2018-01-01_2018-12-31.PData')
lai_fc8 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_8/LAI_2018-01-01_2018-12-31.PData')
lai_fc9 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_9/LAI_2018-01-01_2018-12-31.PData')

swi_fc2 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_2/WG2_2018-01-01_2018-12-31.PData')
swi_fc3 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_3/WG2_2018-01-01_2018-12-31.PData')
swi_fc4 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_4/WG2_2018-01-01_2018-12-31.PData')
swi_fc5 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_5/WG2_2018-01-01_2018-12-31.PData')
swi_fc6 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_6/WG2_2018-01-01_2018-12-31.PData')
swi_fc7 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_7/WG2_2018-01-01_2018-12-31.PData')
swi_fc8 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_8/WG2_2018-01-01_2018-12-31.PData')
swi_fc9 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/results/fc_9/WG2_2018-01-01_2018-12-31.PData')




### Import entire folders of yearly data
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    return parts
'''
### It is required to have as many 'read_pickle's as years in the folder 
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_SGF/results/ol/LAI_2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
lai_fc1_us = pd.concat(df)


'''

###################################################################
### Data Processing 
###################################################################

patch_frac = getModelPatchFractions('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/sfx-trip/pgd_prep/')

patch_frac[patch_frac==999.0]=0.

c4_t = patch_frac[7]
c4 = c4_t.reshape((40,100))

w = np.where(c4_t>0.5)
'''
lai_obs = lai_obs[w[0][:]].mean(axis=1)
lai_ol = lai_ol[w[0][:]].mean(axis=1)
lai_ekf = lai_ekf[w[0][:]].mean(axis=1)
lai_fc2 = lai_fc2[w[0][:]].mean(axis=1)
lai_fc3 = lai_fc3[w[0][:]].mean(axis=1)
lai_fc4 = lai_fc4[w[0][:]].mean(axis=1)
lai_fc5 = lai_fc5[w[0][:]].mean(axis=1)
lai_fc6 = lai_fc6[w[0][:]].mean(axis=1)
lai_fc7 = lai_fc7[w[0][:]].mean(axis=1)
lai_fc8 = lai_fc8[w[0][:]].mean(axis=1)
lai_fc9 = lai_fc9[w[0][:]].mean(axis=1)

swi_obs = swi_obs[w[0][:]].mean(axis=1)
swi_ol = swi_ol[w[0][:]].mean(axis=1)
swi_ekf = swi_ekf[w[0][:]].mean(axis=1)
swi_fc2 = swi_fc2[w[0][:]].mean(axis=1)
swi_fc3 = swi_fc3[w[0][:]].mean(axis=1)
swi_fc4 = swi_fc4[w[0][:]].mean(axis=1)
swi_fc5 = swi_fc5[w[0][:]].mean(axis=1)
swi_fc6 = swi_fc6[w[0][:]].mean(axis=1)
swi_fc7 = swi_fc7[w[0][:]].mean(axis=1)
swi_fc8 = swi_fc8[w[0][:]].mean(axis=1)
swi_fc9 = swi_fc9[w[0][:]].mean(axis=1)
'''

lai_obs = lai_obs.mean(axis=1)
lai_ol = lai_ol.mean(axis=1)
lai_ekf = lai_ekf.mean(axis=1)
lai_fc2 = lai_fc2.mean(axis=1)
lai_fc3 = lai_fc3.mean(axis=1)
lai_fc4 = lai_fc4.mean(axis=1)
lai_fc5 = lai_fc5.mean(axis=1)
lai_fc6 = lai_fc6.mean(axis=1)
lai_fc7 = lai_fc7.mean(axis=1)
lai_fc8 = lai_fc8.mean(axis=1)
lai_fc9 = lai_fc9.mean(axis=1)

swi_obs = swi_obs.mean(axis=1)
swi_ol = swi_ol.mean(axis=1)
swi_ekf = swi_ekf.mean(axis=1)
swi_fc2 = swi_fc2.mean(axis=1)
swi_fc3 = swi_fc3.mean(axis=1)
swi_fc4 = swi_fc4.mean(axis=1)
swi_fc5 = swi_fc5.mean(axis=1)
swi_fc6 = swi_fc6.mean(axis=1)
swi_fc7 = swi_fc7.mean(axis=1)
swi_fc8 = swi_fc8.mean(axis=1)
swi_fc9 = swi_fc9.mean(axis=1)



lai_obs = lai_obs.dropna()
lai_tmp = lai_ol.copy()*np.nan
lai_tmp.ix[lai_obs.index] = lai_obs

lai_fc2 = lai_fc2[~np.isnan(lai_tmp)].dropna()
lai_fc3 = lai_fc3[~np.isnan(lai_tmp)].dropna()
lai_fc4 = lai_fc4[~np.isnan(lai_tmp)].dropna()
lai_fc5 = lai_fc5[~np.isnan(lai_tmp)].dropna()
lai_fc6 = lai_fc6[~np.isnan(lai_tmp)].dropna()
lai_fc7 = lai_fc7[~np.isnan(lai_tmp)].dropna()
lai_fc8 = lai_fc8[~np.isnan(lai_tmp)].dropna()
lai_fc9 = lai_fc9[~np.isnan(lai_tmp)].dropna()
lai_ol = lai_ol[~np.isnan(lai_tmp)].dropna()
lai_ekf = lai_ekf[~np.isnan(lai_tmp)].dropna()

swi_obs = swi_obs.dropna()
swi_tmp = swi_ol.copy()*np.nan
swi_tmp.ix[swi_obs.index] = swi_obs

swi_fc2 = swi_fc2[~np.isnan(swi_tmp)].dropna()
swi_fc3 = swi_fc3[~np.isnan(swi_tmp)].dropna()
swi_fc4 = swi_fc4[~np.isnan(swi_tmp)].dropna()
swi_fc5 = swi_fc5[~np.isnan(swi_tmp)].dropna()
swi_fc6 = swi_fc6[~np.isnan(swi_tmp)].dropna()
swi_fc7 = swi_fc7[~np.isnan(swi_tmp)].dropna()
swi_fc8 = swi_fc8[~np.isnan(swi_tmp)].dropna()
swi_fc9 = swi_fc9[~np.isnan(swi_tmp)].dropna()
swi_ol = swi_ol[~np.isnan(swi_tmp)].dropna()
swi_ekf = swi_ekf[~np.isnan(swi_tmp)].dropna()


#swi_obs = swi_obs.dropna()



#sys.exit()
###################################################################
### Computer Z Scores
###################################################################


###################################################################
### Computer Correlations
###################################################################

laicor_ol_obs = np.corrcoef(lai_ol.values,lai_obs.values)
laicor_ekf_obs = np.corrcoef(lai_ekf.values,lai_obs.values)
laicor_fc2_obs = np.corrcoef(lai_fc2.values,lai_obs.values)
laicor_fc3_obs = np.corrcoef(lai_fc3.values,lai_obs.values)
laicor_fc4_obs = np.corrcoef(lai_fc4.values,lai_obs.values)
laicor_fc5_obs = np.corrcoef(lai_fc5.values,lai_obs.values)
laicor_fc6_obs = np.corrcoef(lai_fc6.values,lai_obs.values)
laicor_fc7_obs = np.corrcoef(lai_fc7.values,lai_obs.values)
laicor_fc8_obs = np.corrcoef(lai_fc8.values,lai_obs.values)
laicor_fc9_obs = np.corrcoef(lai_fc9.values,lai_obs.values)

swicor_ol_obs = np.corrcoef(swi_ol.values,swi_obs.values)
swicor_ekf_obs = np.corrcoef(swi_ekf.values,swi_obs.values)
swicor_fc2_obs = np.corrcoef(swi_fc2.values,swi_obs.values)
swicor_fc3_obs = np.corrcoef(swi_fc3.values,swi_obs.values)
swicor_fc4_obs = np.corrcoef(swi_fc4.values,swi_obs.values)
swicor_fc5_obs = np.corrcoef(swi_fc5.values,swi_obs.values)
swicor_fc6_obs = np.corrcoef(swi_fc6.values,swi_obs.values)
swicor_fc7_obs = np.corrcoef(swi_fc7.values,swi_obs.values)
swicor_fc8_obs = np.corrcoef(swi_fc8.values,swi_obs.values)
swicor_fc9_obs = np.corrcoef(swi_fc9.values,swi_obs.values)

print("LAI Correlation OL vs Observations: " + str(round(laicor_ol_obs[0,1],4)))
print("LAI Correlation EKF vs Observations: " + str(round(laicor_ekf_obs[0,1],4)))
print("LAI Correlation FC2 vs Observations: " + str(round(laicor_fc2_obs[0,1],4)))
print("LAI Correlation FC3 vs Observations: " + str(round(laicor_fc3_obs[0,1],4)))
print("LAI Correlation FC4 vs Observations: " + str(round(laicor_fc4_obs[0,1],4)))
print("LAI Correlation FC5 vs Observations: " + str(round(laicor_fc5_obs[0,1],4)))
print("LAI Correlation FC6 vs Observations: " + str(round(laicor_fc6_obs[0,1],4)))
print("LAI Correlation FC7 vs Observations: " + str(round(laicor_fc7_obs[0,1],4)))
print("LAI Correlation FC8 vs Observations: " + str(round(laicor_fc8_obs[0,1],4)))
print("LAI Correlation FC9 vs Observations: " + str(round(laicor_fc9_obs[0,1],4)))

print("-----------------------------------------------")

print("SWI Correlation OL vs Observations: " + str(round(swicor_ol_obs[0,1],4)))
print("SWI Correlation EKF vs Observations: " + str(round(swicor_ekf_obs[0,1],4)))
print("SWI Correlation FC2 vs Observations: " + str(round(swicor_fc2_obs[0,1],4)))
print("SWI Correlation FC3 vs Observations: " + str(round(swicor_fc3_obs[0,1],4)))
print("SWI Correlation FC4 vs Observations: " + str(round(swicor_fc4_obs[0,1],4)))
print("SWI Correlation FC5 vs Observations: " + str(round(swicor_fc5_obs[0,1],4)))
print("SWI Correlation FC6 vs Observations: " + str(round(swicor_fc6_obs[0,1],4)))
print("SWI Correlation FC7 vs Observations: " + str(round(swicor_fc7_obs[0,1],4)))
print("SWI Correlation FC8 vs Observations: " + str(round(swicor_fc8_obs[0,1],4)))
print("SWI Correlation FC9 vs Observations: " + str(round(swicor_fc9_obs[0,1],4)))

#sys.exit()
###################################################################
### Graphing Functions 
###################################################################
### plt.subplots(num_rows,num_cols)
#fig, axs = plt.subplots(2,1)

raw_input("Press Enter to continue ...")
'''
#axs[0].set_title('LAI OL, EKF, and Forecasts (patches with >0.5 c4) over Nebraska')
axs[0].set_title('LAI OL, EKF, and Forecasts over Nebraska')
axs[0].plot(lai_ol,label='lai_ol',linewidth=2.5,color='blue')
axs[0].plot(lai_ekf,label='lai_ekf',linewidth=2.5,color='red')
axs[0].plot(lai_fc2,label='lai_fc2',linestyle='--')
#axs[0].plot(lai_fc3,label='lai_fc3',linestyle='--')
axs[0].plot(lai_fc4,label='lai_fc4',linestyle='--')
#axs[0].plot(lai_fc5,label='lai_fc5',linestyle='--')
axs[0].plot(lai_fc6,label='lai_fc6',linestyle='--')
#axs[0].plot(lai_fc7,label='lai_fc7',linestyle='--')
axs[0].plot(lai_fc8,label='lai_fc8',linestyle='--')
#axs[0].plot(lai_fc9,label='lai_fc9',linestyle='--')
axs[0].plot(lai_obs,label='lai_Obs',marker='*',linewidth=0,color='green',markersize=10)
axs[0].legend(loc='upper left')
#axs[0].set_ylim([-2.5,2.5])
axs[0].axhline(color='black',linewidth=.5)
axs[0].set_ylabel('LAI [$m^2$/$m^2$]')


#axs[1].set_title('SWI OL, EKF, and Forecasts (patches with >0.5 c4) over Nebraska')
axs[1].set_title('SWI OL, EKF, and Forecasts over Nebraska')
axs[1].plot(swi_ol,label='swi_ol',linewidth=2.5,color='blue')
axs[1].plot(swi_ekf,label='swi_ekf',linewidth=2.5,color='red')
axs[1].plot(swi_fc2,label='swi_fc2',linestyle='--')
#axs[1].plot(swi_fc3,label='swi_fc3',linestyle='--')
axs[1].plot(swi_fc4,label='swi_fc4',linestyle='--')
#axs[1].plot(swi_fc5,label='swi_fc5',linestyle='--')
axs[1].plot(swi_fc6,label='swi_fc6',linestyle='--')
#axs[1].plot(swi_fc7,label='swi_fc7',linestyle='--')
axs[1].plot(swi_fc8,label='swi_fc8',linestyle='--')
#axs[1].plot(swi_fc9,label='swi_fc9',linestyle='--')
axs[1].plot(swi_obs,label='swi_Obs',marker='*',linewidth=0,color='green',markersize=10)
axs[1].legend(loc='upper left')
axs[1].axhline(color='black',linewidth=.5)
axs[1].set_ylabel('SWI [%]')
axs[1].set_ylim([0.1,.45])

'''

plt.title('LAI over Nebraska')
plt.plot(lai_ol,label='Open-Loop',linewidth=2.5,color='blue')
plt.plot(lai_ekf,label='Analysis',linewidth=2.5,color='red')
#plt.plot(lai_fc2,label='2 Day FC',linestyle='--')
#plt.plot(lai_fc3,label='lai_fc3',linestyle='--')
#plt.plot(lai_fc4,label='4 Day FC',linestyle='--')
#plt.plot(lai_fc5,label='lai_fc5',linestyle='--')
#plt.plot(lai_fc6,label='6 Day FC',linestyle='--')
#plt.plot(lai_fc7,label='lai_fc7',linestyle='--')
#plt.plot(lai_fc8,label='8 Day FC',linestyle='--')
plt.plot(lai_fc9,label='lai_fc9',linestyle='--',color='magenta')
plt.plot(lai_obs,label='Observed LAI',marker='*',linewidth=0,color='green',markersize=10)
plt.legend(loc='upper left')
#plt.set_ylim([-2.5,2.5])
plt.axhline(color='black',linewidth=.5)
plt.ylabel('LAI [$m^2$/$m^2$]')



plt.axhline(color='black',linewidth=.5)
plt.legend(loc='upper left')

plt.rcParams.update({'font.size': 12})
#fig.tight_layout()
#save('vegeoFigures/LAI_NE10_2018_FC', ext='ps', close=True, verbose=True)
plt.show()

