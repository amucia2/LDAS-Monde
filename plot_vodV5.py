import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpl_toolkits.basemap import Basemap
import pandas as pd
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr
import pylab as pl
import sys, glob, os, re
import math
import time, tqdm
from tqdm import tqdm

def filter_nan(s,o):
    """
    this functions removed the data  from simulated and observed data
    whereever the observed data contains nan
    this is used by all other functions, otherwise they will produce nan as
    output
    """
    data = np.array([s.flatten(),o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]
    #data = data[~np.isnan(data)]
    return data[:,0],data[:,1]

def rmse(s,o):
    """
    Root Mean Squared Error
    input:
    s: simulated
    o: observed
    output: 
    rmses: root mean squared error
    """
    s,o = filter_nan(s,o)
    return np.sqrt(np.mean((s-o)**2))

def correlation(s,o):      
    """
    correlation coefficient
    input:
    s: simulated
    o: observed
    output:
    correlation: correlation coefficient
    """
    s,o = filter_nan(s,o)
    if s.size == 0:
        corr = np.NaN
    else:
        #corr = np.corrcoef(o, s)[0,1]  
        corr = pearsonr(o, s)
    return corr

def bias(s, o):
    """
        Bias
        input:
        s: simulated
        o: observed
        output:
        bias: bias
        """
    s,o = filter_nan(s,o)
    return np.mean(s-o)

def ubrmse(s,o):
    """
    Unbiased root Mean Squared Error
    input:
    s: simulated
    o: observed
    output: 
    ubrmses: unbiased root mean squared error
    """
    s,o = filter_nan(s,o)
    n = len(o)
    if (n!=0):
        o_moy = np.mean(o)
        s_moy = np.mean(s)
        somme = 0.
        for i in range(n):
            somme = somme + ((s[i]-s_moy)-(o[i]-o_moy))**2
        return np.sqrt(somme/n)
    else:
        return np.nan

FILL = 999.

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
    if isinstance(var_names,str): var_names = [var_names]
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

def getModelPatchFractions(mod_pgd_path):
    '''
    Read the patch fractions of model grid points from file.
    Arguments:
        mod_pgd_path (str): directory of PREP.nc
    Returns:
        mod_patch_frac (numpy.ndarray): patch fractions of each model grid point; shape = n_points,n_patches
    '''
    mod_pgd_path = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/sfx-trip/pgd_prep/'
    prep_file = mod_pgd_path+'PREP.nc'
    if not os.path.isfile(prep_file):
        raise Exception(prep_file+' does not exist.')
    frac = readFromFile(prep_file,['PATCH'])[0]
    return frac.reshape((frac.shape[0],-1))

patch_frac = getModelPatchFractions('')
patch_frac[patch_frac==999.0]=0.

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

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

import pickle

#tt1 = pd.read_pickle('../CONUS_025_VODX_2003_2006.PData') 
#tt2 = pd.read_pickle('../CONUS_025_VODX_2007_2011.PData')
#tt3 = pd.read_pickle('../CONUS_025_VODX_2012_2016.PData')
#tt_lai = pd.concat([tt1,tt2,tt3],axis=1)
#tt_vod = tt_lai
#sys.exit()

#'''
#tt_lai = pd.read_pickle('CONUS_025_VODX_2003_2016.PData')
#tt_vod = pd.read_pickle('CONUS_025_VODX_2003_2016.PData')

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts

### Importing all the info seperately, that is not in a Panel, but will require to rewrite massive amounts of code that only uses Panel selection
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/Panels/LAI_V2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),pd.read_pickle(file[18])]
lai_cgls_ = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/Panels/LAI_20*.PData'),key=numericalSort)
df = [pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),pd.read_pickle(file[18])]
lai_isba = pd.concat(df)

date = pd.date_range(start="2003-01-01 09:00:00",end='2018-12-31 09:00:00',freq='D')
lai_cgls = lai_cgls_.reindex(date,fill_value=np.nan)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/sfx-trip/VODC*.PData'),key=numericalSort)
df = [pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),pd.read_pickle(file[18])]
vodx = pd.concat(df)

#mvodx = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/LAIfromVODCAX_V4.PData')
### Currently cannot import Panels and pickles are all f**ked up too - I'm frustrated, if you can't tell...
#tt_lai = pd.read_pickle('CONUS_2000_2018_LAI_VOD_Panel.PData')
#tt_vod = pd.read_pickle('CONUS_2000_2018_LAI_VOD_Panel.PData')
#tt_lai = pickle.load('CONUS_2000_2018_LAI_VOD_Panel.PData',protocol=2)


'''
Model : LAI ISBA
Analysis : LAI CGLS
Obs : VOD
'''

#test = lai_cgls.mean(axis=1)[~np.isnan(lai_cgls.mean(axis=1))]

#for i in range((tt_vod.shape[2])):
#    vodx[i] = pd.rolling_mean(vodx[i],30)

'''

fig, ax1 =  plt.subplots()
ax1.plot(test,label='LAI obs.',color='r',linewidth=0,marker='*')
ax1.plot(lai_isba.mean(axis=1),label='LAI isba',color='red')
plt.rcParams.update({'font.size': 15})
ax1.legend(loc='best')
ax1.set_ylabel('LAI',color='red')
ax2 = ax1.twinx()
ax2.plot(vodx.mean(axis=1),label='Cband VOD',color='blue')
ax2.set_ylabel('VOD',color='blue')
plt.show()
'''



'''

corr_tmp_VOD_LAIisba = vodx.corrwith(lai_isba,axis=0)
corr_tmp_VOD_LAIobs = vodx.corrwith(lai_cgls,axis=0)
corr_tmp_LAIisba_LAIobs = lai_isba.corrwith(lai_cgls,axis=0)

corr_tmp_VOD_LAIisba.values[pl.where(pl.isinf(corr_tmp_VOD_LAIisba.values))] = pl.nan
corr_tmp_VOD_LAIobs.values[pl.where(pl.isinf(corr_tmp_VOD_LAIobs.values))] = pl.nan
corr_tmp_LAIisba_LAIobs.values[pl.where(pl.isinf(corr_tmp_LAIisba_LAIobs.values))] = pl.nan

corr_tmp_MVOD_LAIobs = mvodx.corrwith(lai_cgls,axis=0)
corr_tmp_MVOD_LAIobs.values[pl.where(pl.isinf(corr_tmp_MVOD_LAIobs.values))] = pl.nan

v = lai_cgls.mean().values.reshape((140,280))
v2 = vodx.mean().values.reshape((140,280))
v3 = mvodx.mean().values.reshape((140,280))
#v4 = corr_tmp_VOD_LAIobs.values.reshape((140,280))
v4 = corr_tmp_MVOD_LAIobs.values.reshape((140,280))
'''
'''
plt.subplot(221) ; plt.imshow(v,origin='lower',vmin=0,vmax=6,cmap="RdYlGn") ; plt.title('GEOV2 LAI') ; cbar = plt.colorbar(orientation='horizontal',fraction=0.08, pad=0.1)
plt.subplot(222) ; plt.imshow(v2,origin='lower',vmin=0,cmap="RdYlGn"); plt.title('VODCA VOD-X') ; cbar = plt.colorbar(orientation='horizontal',fraction=0.08, pad=0.1)
plt.subplot(223) ; plt.imshow(v3,origin='lower',vmin=0,vmax=6,cmap="RdYlGn"); plt.title('Matched VOD-X') ; cbar = plt.colorbar(orientation='horizontal',fraction=0.08, pad=0.1)
plt.subplot(224) ; plt.imshow(v4,origin='lower',vmin=0,vmax=1,cmap="RdYlGn"); plt.title('R : Matched VOD-X vs LAI GEOV2') ; cbar = plt.colorbar(orientation='horizontal',fraction=0.08, pad=0.1)
plt.rcParams.update({'font.size': 10})
plt.show()
'''
'''
v = corr_tmp_VOD_LAIisba.values.reshape((140,280))
v2 = corr_tmp_VOD_LAIobs.values.reshape((140,280))
v3 = corr_tmp_LAIisba_LAIobs.values.reshape((140,280))
plt.subplot(221) ; plt.imshow(v,origin='lower',vmin=0,vmax=1.,cmap='bwr') ; plt.title('R : VODCA-X vs. LAI ISBA') ; cbar = plt.colorbar(orientation='horizontal')
plt.subplot(222) ; plt.imshow(v2,origin='lower',vmin=0,vmax=1.,cmap='bwr'); plt.title('R : VODCA-X vs. LAI GEOV2') ; cbar = plt.colorbar(orientation='horizontal')
plt.subplot(223) ; plt.imshow(v3,origin='lower',vmin=0,vmax=1,cmap='bwr'); plt.title('R : LAI ISBA vs. LAI GEOV2') ; cbar = plt.colorbar(orientation='horizontal')
plt.rcParams.update({'font.size': 10})
plt.show()

'''
'''
print('R : VOD vs LAIisba (mean / median): ', corr_tmp_VOD_LAIisba.mean() ,' : ',corr_tmp_VOD_LAIisba.median())
print('R : VOD vs LAIobs (mean / median): ', corr_tmp_VOD_LAIobs.mean() ,' : ',corr_tmp_VOD_LAIobs.median())
print('R : LAIisba vs LAIobs (mean / median): ', corr_tmp_LAIisba_LAIobs.mean() ,' : ',corr_tmp_LAIisba_LAIobs.median())
'''
toto = vodx.copy()*np.nan

#sys.exit()
#'''
#for i in tqdm(range((vodx.shape[2]))):
for i in tqdm(range((vodx.shape[1]))):
#    season=((vodx.index.month>=0))
#    season=((vodx.index.month>=4) & (vodx.index.month<=9))
#    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & ((patch_frac[0,i]+patch_frac[1,i]+patch_frac[2,i]) < 0.80)):
#        try:
#            aa = np.polyfit(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])],\
#                    vodx[i].ix[season][~np.isnan(vodx[i].ix[season])],2)
#            toto['Obs'][i].ix[season] = aa[0]*vodx[i].ix[season]**2+aa[1]*vodx[i].ix[season]+aa[2]
#            toto['Model'][i] = correlation(toto['Obs'][i].ix[season].values,vodx[i].ix[season].values)[0]
#            print i, correlation(toto['Obs'][i].ix[season].values,vodx[i].ix[season].values)[0]
#        except ValueError:
#            print('ValueError')
#
#    season2=((vodx.index.month<=3) | (vodx.index.month>=10))
#    if ((len(vodx[i].ix[season2][~np.isnan(vodx[i].ix[season2])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.10)):
#        try:
#            aa2 = np.polyfit(vodx[i].ix[season2][~np.isnan(vodx[i].ix[season2])],\
#                    vodx[i].ix[season2][~np.isnan(vodx[i].ix[season2])],2)
#            toto['Obs'][i].ix[season2] = aa2[0]*vodx[i].ix[season2]**2+aa2[1]*vodx[i].ix[season2]+aa2[2]
#            toto['Model'][i] = correlation(toto['Obs'][i].ix[season2].values,vodx[i].ix[season2].values)[0]
#            print i, correlation(toto['Obs'][i].ix[season2].values,vodx[i].ix[season2].values)[0]
#        except ValueError:
#            print('ValueError')

#    if ( (len(toto['Obs'][i][~np.isnan(toto['Obs'][i])])> 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.10)):
#        toto['Model'][i] = correlation(toto['Obs'][i].values,vodx[i].values)[0]
#        toto['Analysis'][i] = ubrmse(toto['Obs'][i].values,vodx[i].values)
#        print i, toto['Model'][i].ix[0],toto['Analysis'][i].ix[0]


    season=((vodx.index.month==12) | (vodx.index.month==1) | (vodx.index.month==2))
    season2=(vodx.index.month==1)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            #aa = np.polyfit(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])],\
            #        vodx[i].ix[season][~np.isnan(vodx[i].ix[season])],2)
            #toto['Obs'][i].ix[season2] = aa[0]*vodx[i].ix[season2]**2+aa[1]*vodx[i].ix[season2]+aa[2]
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

    season=((vodx.index.month==1) | (vodx.index.month==2) | (vodx.index.month==3))
    season2=(vodx.index.month==2)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            #aa = np.polyfit(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])],\
            #        vodx[i].ix[season][~np.isnan(vodx[i].ix[season])],2)
            #toto['Obs'][i].ix[season2] = aa[0]*vodx[i].ix[season2]**2+aa[1]*vodx[i].ix[season2]+aa[2]
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

    season=((vodx.index.month==2) | (vodx.index.month==3) | (vodx.index.month==4))
    season2=(vodx.index.month==3)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

    season=((vodx.index.month==3) | (vodx.index.month==4) | (vodx.index.month==5))
    season2=(vodx.index.month==4)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

    season=((vodx.index.month==4) | (vodx.index.month==5) | (vodx.index.month==6))
    season2=(vodx.index.month==5)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

    season=((vodx.index.month==5) | (vodx.index.month==6) | (vodx.index.month==7))
    season2=(vodx.index.month==6)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

    season=((vodx.index.month==6) | (vodx.index.month==7) | (vodx.index.month==8))
    season2=(vodx.index.month==7)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

    season=((vodx.index.month==7) | (vodx.index.month==8) | (vodx.index.month==9))
    season2=(vodx.index.month==8)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

    season=((vodx.index.month==8) | (vodx.index.month==9) | (vodx.index.month==10))
    season2=(vodx.index.month==9)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

    season=((vodx.index.month==9) | (vodx.index.month==10) | (vodx.index.month==11))
    season2=(vodx.index.month==10)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

    season=((vodx.index.month==10) | (vodx.index.month==11) | (vodx.index.month==12))
    season2=(vodx.index.month==11)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

    season=((vodx.index.month==11) | (vodx.index.month==12) | (vodx.index.month==1))
    season2=(vodx.index.month==12)
    if ((len(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]) > 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        try:
            b = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).std()/\
                    (vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).std()
            a = (lai_cgls[i].ix[season][~np.isnan(lai_cgls[i].ix[season])]).mean() - \
                    b*(vodx[i].ix[season][~np.isnan(vodx[i].ix[season])]).mean()
            toto[i].ix[season2] = b*(vodx[i].ix[season2]) + a
        except ValueError:
            print('ValueError')

toto[toto>10]=np.nan
toto[toto<0]=np.nan

'''
for i in range((vodx.shape[2])):
    if ( (len(toto['Obs'][i][~np.isnan(toto['Obs'][i])])> 0) & (patch_frac[0,i]+patch_frac[1,i] < 0.80)):
        #toto['Obs'][i] = pd.rolling_mean(toto['Obs'][i],30)
        toto['Model'][i] = ubrmse(toto['Obs'][i].values,vodx[i].values)
        #toto['Model'][i] = correlation(toto['Obs'][i].values,vodx[i].values)[0]
        toto['Analysis'][i] = ubrmse(toto['Obs'][i].values,lai_isba[i].values)
        #toto['Analysis'][i] = correlation(toto['Obs'][i].values,lai_isba[i].values)[0]
        print(i, toto['Model'][i].ix[0],toto['Analysis'][i].ix[0])
'''
'''
v = toto['Obs'].mean(axis=0).values.reshape((140,280))
plt.subplot(221) ; plt.imshow(v,origin='lower') ; plt.title('LAI from VOD')
cbar = plt.colorbar(orientation='horizontal')
v = lai_cgls.mean(axis=0).values.reshape((140,280))
plt.subplot(222) ; plt.imshow(v,origin='lower') ; plt.title('LAI GEOV2')
cbar = plt.colorbar(orientation='horizontal')
v = toto['Model'].ix[0].values.reshape((140,280))
plt.subplot(223) ; plt.imshow(v,origin='lower') ; plt.title('ubRMSD : LAI GEOV2 vs LAI from VOD')
cbar = plt.colorbar(orientation='horizontal')
v = toto['Analysis'].ix[0].values.reshape((140,280))
plt.subplot(224) ; plt.imshow(v,origin='lower') ; plt.title('ubRMSD : LAI ISBA vs LAI from VOD')
#plt.subplot(224) ; plt.imshow(v,origin='lower') ; plt.title('R : LAI ISBA vs LAI from VOD')
cbar = plt.colorbar(orientation='horizontal')
plt.rcParams.update({'font.size': 10})
plt.show()
#save('VOD_FIT_CDF', ext="ps", close=True, verbose=True)

v = toto['Obs'].mean(axis=0).values.reshape((140,280))
plt.subplot(221) ; plt.imshow(v,origin='lower') ; plt.title('LAI from VOD')
cbar = plt.colorbar(orientation='horizontal')
v = lai_cgls.mean(axis=0).values.reshape((140,280))
plt.subplot(222) ; plt.imshow(v,origin='lower') ; plt.title('LAI GEOV2')
cbar = plt.colorbar(orientation='horizontal')
v = toto['Obs'].mean(axis=0).values.reshape((140,280))-lai_cgls.mean(axis=0).values.reshape((140,280))
plt.subplot(223) ; plt.imshow(v,origin='lower')
cbar = plt.colorbar(orientation='horizontal')
#v = toto['Model'].ix[0].reshape((140,280))
#plt.subplot(223) ; plt.imshow(v,origin='lower') ; plt.title('ubRMSD : LAI GEOV2 vs LAI from VOD')
#cbar = plt.colorbar(orientation='horizontal')
#v = toto['Analysis'].ix[0].reshape((140,280))
#plt.subplot(224) ; plt.imshow(v,origin='lower') ; plt.title('ubRMSD : LAI ISBA vs LAI from VOD')
##plt.subplot(224) ; plt.imshow(v,origin='lower') ; plt.title('R : LAI ISBA vs LAI from VOD')
#cbar = plt.colorbar(orientation='horizontal')
plt.rcParams.update({'font.size': 10})
plt.show()
'''

'''
#plt.plot(pd.rolling_mean(toto['Obs'],window=30,min_periods=5,center=True),\
#        label='LAI from VOD',marker='*',linewidth=0) ; plt.plot(test, label='LAI_obs',marker='*',linewidth=0) 
#plt.plot(lai_isba.mean(axis=1),label='LAI_ISBA') ; plt.legend(loc='best') ; plt.show()

#for ii in [1,2,3,4,5,6,7,8,9,10,11,12]:
#    season = (vodx.index.month == ii)
#
#    corr_tmp_VOD_LAIisba = vodx.ix[season,:].corrwith(lai_isba.ix[season,:],axis=0)
#    corr_tmp_VOD_LAIobs = vodx.ix[season,:].corrwith(vodx.ix[season,:],axis=0)
#    corr_tmp_LAIisba_LAIobs = lai_isba.ix[season,:].corrwith(vodx.ix[season,:],axis=0)
#
#    corr_tmp_VOD_LAIisba.values[pl.where(pl.isinf(corr_tmp_VOD_LAIisba.values))] = pl.nan
#    corr_tmp_VOD_LAIobs.values[pl.where(pl.isinf(corr_tmp_VOD_LAIobs.values))] = pl.nan
#    corr_tmp_LAIisba_LAIobs.values[pl.where(pl.isinf(corr_tmp_LAIisba_LAIobs.values))] = pl.nan
#
#    v = corr_tmp_VOD_LAIisba.values.reshape((140,280))
#    v2 = corr_tmp_VOD_LAIobs.values.reshape((139,280))
#    v3 = corr_tmp_LAIisba_LAIobs.values.reshape((140,280))
#    plt.subplot(221) ; plt.imshow(v,origin='lower') ; plt.title('R : VOD vs. LAI ISBA') ; cbar = plt.colorbar(orientation='horizontal')
#    plt.subplot(222) ; plt.imshow(v2,origin='lower'); plt.title('R : VOD vs. LAI OBS') ; cbar = plt.colorbar(orientation='horizontal')
#    plt.subplot(223) ; plt.imshow(v3,origin='lower'); plt.title('R : LAI ISBA vs. LAI OBS') ; cbar = plt.colorbar(orientation='horizontal')
#    plt.rcParams.update({'font.size': 12})
#    plt.show()
#    print 'R : VOD vs LAIisba (mean / median): ', corr_tmp_VOD_LAIisba.mean() ,' : ',corr_tmp_VOD_LAIisba.median()
#    print 'R : VOD vs LAIobs (mean / median): ', corr_tmp_VOD_LAIobs.mean() ,' : ',corr_tmp_VOD_LAIobs.median()
#    print 'R : LAIisba vs LAIobs (mean / median): ', corr_tmp_LAIisba_LAIobs.mean() ,' : ',corr_tmp_LAIisba_LAIobs.median()

'''

#A = pd.rolling_mean(toto,window=30,min_periods=5,center=True)
A = toto.rolling(window=30,min_periods=5,center=True).mean()
#A.to_pickle('LAIfromVODCAX_V5_rollingMean.PData')
#A = toto['Obs']
#toto['Obs'].to_pickle('LAIfromVODCAX_V4.PData') 
#A = toto['Obs']
#sys.exit()
#A = pd.read_pickle('LAIfromVODCAX_V4.PData') 


for d in tqdm(range(len(A))): 
    if d < 5844:
        #print(A.index[d])
        try:
            #print(min(A.ix[d][~np.isnan(A.ix[d])]), max(A.ix[d][~np.isnan(A.ix[d])]))
            # Clean nan values
            A.ix[d][np.invert(~np.isnan(A.ix[d]))] = 999
            # Specify date for file name
            YY = str(A.index[d])[2:4]
            MM = str(A.index[d])[5:7]
            DD = str(A.index[d])[8:10]
            # create file
            f= open('VODC_2_LAI_cdf/CANARI/CANARI_NATURE_'+YY+MM+DD+'H09.DAT','a+')
            # write in file
            for i in range(len(A.ix[d])):
                #print(A.ix[d][i])
                f.write(str(A.ix[d][i])+'\n')
                #time.sleep(0.01)
            # close file
            f.close()
        except ValueError:
            print('No values')
            # Clean nan values
            A.ix[d][np.invert(~np.isnan(A.ix[d]))] = 999
            # Specify date for file name
            YY = str(A.index[d])[2:4]
            MM = str(A.index[d])[5:7]
            DD = str(A.index[d])[8:10]
            # create file
            f= open('VODC_2_LAI_cdf/CANARI/CANARI_NATURE_'+YY+MM+DD+'H09.DAT','a+')
            # write in file
            for i in range(len(A.ix[d])):
                f.write(str(A.ix[d][i])+'\n')
            # close file
            f.close()

#A.to_pickle('LAIfromVODCAX_V5_rollingMean.PData')

import joblib
joblib.dump(A,'LAIfromVODCAC_V7_2003-2018.PData')
'''        

a= []

for d in range(len(A)): 
    #if d > 2191:
    print(A.index[d])
    try:
        print(min(A.ix[d][~np.isnan(A.ix[d])]), max(A.ix[d][~np.isnan(A.ix[d])]))
        # Clean nan values
        print(len(A.ix[d]))
        A.ix[d][np.invert(~np.isnan(A.ix[d]))] = 999
        # Specify date for file name
        YY = str(A.index[d])[2:4]
        MM = str(A.index[d])[5:7]
        DD = str(A.index[d])[8:10]
        # create file
        time.sleep(0.01)
        f= open('VOD_2_LAI_cdf/CANARI/CANARI_NATURE_'+YY+MM+DD+'H09.DAT','w+')
        print(len(A.ix[d]))
        # write in file
        for i in range(len(A.ix[d])):
            #print(A.ix[d][i])
            f.write(str(A.ix[d][i])+'\n')
        # close file
        time.sleep(0.01)
        f.close()
    except ValueError:
        print('No values')
        # Clean nan values
        A.ix[d][np.invert(~np.isnan(A.ix[d]))] = 999
        # Specify date for file name
        YY = str(A.index[d])[2:4]
        MM = str(A.index[d])[5:7]
        DD = str(A.index[d])[8:10]
        # create file
        f= open('VOD_2_LAI_cdf/CANARI/CANARI_NATURE_'+YY+MM+DD+'H09.DAT','w+')
        # write in file
        for i in range(len(A.ix[d])):
            f.write(str(A.ix[d][i])+'\n')
        # close file
        f.close()

'''
'''
for d in range(len(A)):
    if len(A.ix[d]) != 39200:
        print("Error")




for d in range(len(A)):
    #if d > 2191:
    print(A.index[d])
    try:
        print(min(A.ix[d][~np.isnan(A.ix[d])]), max(A.ix[d][~np.isnan(A.ix[d])]))
        # Clean nan values
        print(len(A.ix[d]))
        A.ix[d][np.invert(~np.isnan(A.ix[d]))] = 999
        # Specify date for file name
        YY = str(A.index[d])[2:4]
        MM = str(A.index[d])[5:7]
        DD = str(A.index[d])[8:10]
        # create file
        #time.sleep(0.01)
        f= open('VOD_2_LAI_cdf/CANARI/CANARI_NATURE_'+YY+MM+DD+'H09.DAT','w+')
        print(len(A.ix[d]))
        # write in file
        w = A.ix[d]
        f.write(str(w)+'\n')
        # close file
        #time.sleep(0.01)
        f.close()
    except ValueError:
        print('No values')
        # Clean nan values
        A.ix[d][np.invert(~np.isnan(A.ix[d]))] = 999
        # Specify date for file name
        YY = str(A.index[d])[2:4]
        MM = str(A.index[d])[5:7]
        DD = str(A.index[d])[8:10]
        # create file
        f= open('VOD_2_LAI_cdf/CANARI/CANARI_NATURE_'+YY+MM+DD+'H09.DAT','w+')
        # write in file
        for i in range(len(A.ix[d])):
            f.write(str(A.ix[d][i])+'\n')
        # close file
        f.close()
'''
