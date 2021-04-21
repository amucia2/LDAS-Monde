# author: Anthony Mucia
# Last modified: March 2021

## patchFrac_figures_Tony.py
import sys, glob, os, re, math, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
import pylab as pl
import scipy.optimize as optim
from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr, linregress
from numpy.polynomial.polynomial import polyfit


def prep_reg(x,y,xlims):
    x = x.mean(axis=1)
    y = y.mean(axis=1)
    x[np.invert(~np.isnan(y))]=np.nan
    y[np.invert(~np.isnan(x))]=np.nan
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    slope, intercept, r_value, p_value, std_err = linregress(x,y)
    ext = np.arange(xlims[0],xlims[1],1)
    line = [slope*xi + intercept for xi in ext]
    return x, y, ext, line

def nancorr(x,y,patch,perc):
    x = x.loc[:,patch >perc].mean(axis=1)
    y = y.loc[:,patch >perc].mean(axis=1)
    y[np.invert(~np.isnan(x))]=np.nan
    x = x.dropna(); y = y.dropna()
    slope, intercept, r_value, p_value, std_err = linregress(x,y)
    return slope, intercept, r_value, p_value, std_err

def filter_nan(s,o):
    s = np.transpose(s)
    o = np.transpose(o)
    s[np.invert(~np.isnan(o))] = np.nan
    o[np.invert(~np.isnan(s))] = np.nan
    s = s[~pd.isnull(s)]
    o = o[~pd.isnull(o)]
    return s,o

def correlation(x,y,patch,perc):
    x = x.loc[:,patch >perc].mean(axis=1)
    y = y.loc[:,patch >perc].mean(axis=1)
    x, y = filter_nan(x,y)
    if x.size == 0:
        corr = np.NaN
    else:
        corr = pearsonr(x, y) 
    return corr



### ECOSG Prep
PREPdir = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/sfx-trip/pgd_prep/'
file_prep = 'PREP.nc'
PREP_set = Dataset(PREPdir+file_prep,'r')
Auxi = PREP_set.variables['PATCH'][:]
patch_frac = Auxi.data
patch_frac[patch_frac == Auxi.fill_value] = np.nan
PREP_set.close()

### ECOII Prep
'''
PREPdir = '/cnrm/vegeo/bonanb/LDAS/LDAS_curr/US_Tony/pgd_prep_ECOII/'
file_prep = 'PREP_2003010100.nc'
PREP_set = Dataset(PREPdir+file_prep,'r')
Auxi = PREP_set.variables['PATCH'][:]
patch_frac = Auxi.data
patch_frac[patch_frac == Auxi.fill_value] = np.nan
PREP_set.close()
'''

### Import Data
### Raw Obs
#vodx = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2/CONUS/VODX.Pdata')
#vodc = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2/CONUS/VODC.Pdata')
lai_cgls = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/CONUS/LAI_CGLS.Pdata')

### Matched Obs
#'''
vodx = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2/CONUS/LAIfromVODCAX_V8_2003-2018.PData')
vodc = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2/CONUS/LAIfromVODCAC_V8_2003-2018.PData')
#'''

### Organize data by season
idx = pd.DatetimeIndex(vodx.index)
lai_cgls = lai_cgls.reindex(idx, fill_value=np.nan)

lai_cgls_win = lai_cgls[lai_cgls.index.month.isin([12,1,2])]
vodx_win = vodx[vodx.index.month.isin([12,1,2])]
vodc_win = vodc[vodc.index.month.isin([12,1,2])]

lai_cgls_spr = lai_cgls[lai_cgls.index.month.isin([3,4,5])]
vodx_spr = vodx[vodx.index.month.isin([3,4,5])]
vodc_spr = vodc[vodc.index.month.isin([3,4,5])]

lai_cgls_sum = lai_cgls[lai_cgls.index.month.isin([6,7,8])]
vodx_sum = vodx[vodx.index.month.isin([6,7,8])]
vodc_sum = vodc[vodc.index.month.isin([6,7,8])]

lai_cgls_aut = lai_cgls[lai_cgls.index.month.isin([9,10,11])]
vodx_aut = vodx[vodx.index.month.isin([9,10,11])]
vodc_aut = vodc[vodc.index.month.isin([9,10,11])]

patch_num = [3,4,6,9]
patch_name = ['Deciduous','Coniferous','c3_crops','Grasslands']
indx = 0 

### Adding all non-vegetation patches together (need to do so carefully because of the nans)
patch_0 = np.nan_to_num(patch_frac[0,:,:].reshape((280*140)))
patch_1 = np.nan_to_num(patch_frac[1,:,:].reshape((280*140)))
patch_2 = np.nan_to_num(patch_frac[2,:,:].reshape((280*140)))
patch_nonveg = patch_0+patch_1+patch_2


patch_3 = patch_frac[3,:,:].reshape((280*140)) ### Deciduous
patch_4 = patch_frac[4,:,:].reshape((280*140)) ### Coniferous
patch_6 = patch_frac[6,:,:].reshape((280*140)) ### C3 Crops
patch_7 = patch_frac[7,:,:].reshape((280*140)) ### C4 Crops
patch_9 = patch_frac[9,:,:].reshape((280*140)) ### C3 Herb
patch_10 = patch_frac[10,:,:].reshape((280*140)) ### Irrigated Crops

###################################################################
###################################################################
### Scatter Plots, with seasons and 4 (2x2) Patch types
###################################################################
###################################################################

imgdir = '/cnrm/vegeo/muciaa/Images/US00/VOD_Comparisons/Patch/V3/'
perc = 0.5
eco = "ECOII"
ecol = "ECOCLIMAP II"
match = "matched"

###################################################################
### VODX
###################################################################

#plt.rcParams['figure.figsize'] = (11,11)
#plt.rcParams['figure.figsize'] = (11,7)
plt.rcParams['figure.figsize'] = (7,11)
fig, axs = plt.subplots(3,2)
fig.suptitle('LAI vs Matched VODX\n{0}'.format(ecol),fontsize=14)
#fig.suptitle('LAI vs VODX\n{0}'.format(ecol),fontsize=14)

axs[0,0].set_title('Deciduous Forests')
axs[0,0].set_xlim(0,6.5)
axs[0,0].set_ylim(0,6.5)
xlims3=(0,6.5)
x3, y3, ext3, line3 = prep_reg(lai_cgls.loc[:,patch_3 >perc],vodx.loc[:,patch_3 >perc],xlims3)
axs[0,0].plot(ext3,line3,color='purple',linestyle='--')
axs[0,0].scatter(lai_cgls.loc[:,patch_3 >perc],vodx.loc[:,patch_3 >perc],c='black',s=2)
axs[0,0].scatter(lai_cgls_win.loc[:,patch_3 >perc].mean(axis=1),vodx_win.loc[:,patch_3 >perc].mean(axis=1),c='blue',s=4)
axs[0,0].scatter(lai_cgls_spr.loc[:,patch_3 >perc].mean(axis=1),vodx_spr.loc[:,patch_3 >perc].mean(axis=1),c='green',s=4)
axs[0,0].scatter(lai_cgls_sum.loc[:,patch_3 >perc].mean(axis=1),vodx_sum.loc[:,patch_3 >perc].mean(axis=1),c='red',s=4)
axs[0,0].scatter(lai_cgls_aut.loc[:,patch_3 >perc].mean(axis=1),vodx_aut.loc[:,patch_3 >perc].mean(axis=1),c='yellow',s=4)
axs[0,0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[0,0].set_ylabel('Matched VODX [m^2/m^-2]',color='red')

axs[0,1].set_title('Coniferous Forests')
axs[0,1].set_xlim(0,6.5)
axs[0,1].set_ylim(0,6.5)
xlims4=(0,6.5)
x4, y4, ext4, line4 = prep_reg(lai_cgls.loc[:,patch_4 >perc],vodx.loc[:,patch_4 >perc],xlims4)
axs[0,1].plot(ext4,line4,color='purple',linestyle='--')
axs[0,1].scatter(lai_cgls.loc[:,patch_4 >perc],vodx.loc[:,patch_4 >perc],c='black',s=2)
axs[0,1].scatter(lai_cgls_win.loc[:,patch_4 >perc].mean(axis=1),vodx_win.loc[:,patch_4 >perc].mean(axis=1),c='blue',s=4)
axs[0,1].scatter(lai_cgls_spr.loc[:,patch_4 >perc].mean(axis=1),vodx_spr.loc[:,patch_4 >perc].mean(axis=1),c='green',s=4)
axs[0,1].scatter(lai_cgls_sum.loc[:,patch_4 >perc].mean(axis=1),vodx_sum.loc[:,patch_4 >perc].mean(axis=1),c='red',s=4)
axs[0,1].scatter(lai_cgls_aut.loc[:,patch_4 >perc].mean(axis=1),vodx_aut.loc[:,patch_4 >perc].mean(axis=1),c='yellow',s=4)
axs[0,1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[0,1].set_ylabel('Matched VODX [m^2/m^-2]',color='red')

axs[1,0].set_title('C3 Crops')
axs[1,0].set_xlim(0,6.5)
axs[1,0].set_ylim(0,6.5)
xlims6=(0,6.5)
x6, y6, ext6, line6 = prep_reg(lai_cgls.loc[:,patch_6 >perc],vodx.loc[:,patch_6 >perc],xlims6)
axs[1,0].plot(ext6,line6,color='purple',linestyle='--')
axs[1,0].scatter(lai_cgls.loc[:,patch_6 >perc],vodx.loc[:,patch_6 >perc],c='black',s=2)
axs[1,0].scatter(lai_cgls_win.loc[:,patch_6 >perc].mean(axis=1),vodx_win.loc[:,patch_6 >perc].mean(axis=1),c='blue',s=4)
axs[1,0].scatter(lai_cgls_spr.loc[:,patch_6 >perc].mean(axis=1),vodx_spr.loc[:,patch_6 >perc].mean(axis=1),c='green',s=4)
axs[1,0].scatter(lai_cgls_sum.loc[:,patch_6 >perc].mean(axis=1),vodx_sum.loc[:,patch_6 >perc].mean(axis=1),c='red',s=4)
axs[1,0].scatter(lai_cgls_aut.loc[:,patch_6 >perc].mean(axis=1),vodx_aut.loc[:,patch_6 >perc].mean(axis=1),c='yellow',s=4)
axs[1,0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[1,0].set_ylabel('Matched VODX [m^2/m^-2]',color='red')

axs[1,1].set_title('C4 Crops')
axs[1,1].set_xlim(0,6.5)
axs[1,1].set_ylim(0,6.5)
xlims7=(0,6.5)
x7, y7, ext7, line7 = prep_reg(lai_cgls.loc[:,patch_7 >perc],vodx.loc[:,patch_7 >perc],xlims7)
axs[1,1].plot(ext7,line7,color='purple',linestyle='--')
axs[1,1].scatter(lai_cgls.loc[:,patch_7 >perc],vodx.loc[:,patch_7 >perc],c='black',s=2)
axs[1,1].scatter(lai_cgls_win.loc[:,patch_7 >perc].mean(axis=1),vodx_win.loc[:,patch_7 >perc].mean(axis=1),c='blue',s=4)
axs[1,1].scatter(lai_cgls_spr.loc[:,patch_7 >perc].mean(axis=1),vodx_spr.loc[:,patch_7 >perc].mean(axis=1),c='green',s=4)
axs[1,1].scatter(lai_cgls_sum.loc[:,patch_7 >perc].mean(axis=1),vodx_sum.loc[:,patch_7 >perc].mean(axis=1),c='red',s=4)
axs[1,1].scatter(lai_cgls_aut.loc[:,patch_7 >perc].mean(axis=1),vodx_aut.loc[:,patch_7 >perc].mean(axis=1),c='yellow',s=4)
axs[1,1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[1,1].set_ylabel('Matched VODX [m^2/m^-2]',color='red')

axs[2,0].set_title('C3 Herbaceous')
axs[2,0].set_xlim(0,6.5)
axs[2,0].set_ylim(0,6.5)
xlims9=(0,6.5)
x9, y9, ext9, line9 = prep_reg(lai_cgls.loc[:,patch_9 >perc],vodx.loc[:,patch_9 >perc],xlims9)
axs[2,0].plot(ext9,line9,color='purple',linestyle='--')
axs[2,0].scatter(lai_cgls.loc[:,patch_9 >perc],vodx.loc[:,patch_9 >perc],c='black',s=2)
axs[2,0].scatter(lai_cgls_win.loc[:,patch_9 >perc].mean(axis=1),vodx_win.loc[:,patch_9 >perc].mean(axis=1),c='blue',s=4)
axs[2,0].scatter(lai_cgls_spr.loc[:,patch_9 >perc].mean(axis=1),vodx_spr.loc[:,patch_9 >perc].mean(axis=1),c='green',s=4)
axs[2,0].scatter(lai_cgls_sum.loc[:,patch_9 >perc].mean(axis=1),vodx_sum.loc[:,patch_9 >perc].mean(axis=1),c='red',s=4)
axs[2,0].scatter(lai_cgls_aut.loc[:,patch_9 >perc].mean(axis=1),vodx_aut.loc[:,patch_9 >perc].mean(axis=1),c='yellow',s=4)
axs[2,0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[2,0].set_ylabel('Matched VODX [m^2/m^-2]',color='red')

axs[2,1].set_title('Irrigated Crops')
axs[2,1].set_xlim(0,6.5)
axs[2,1].set_ylim(0,6.5)
xlims10=(0,6.5)
x10, y10, ext10, line10 = prep_reg(lai_cgls.loc[:,patch_10 >perc],vodx.loc[:,patch_10 >perc],xlims10)
axs[2,1].plot(ext10,line10,color='purple',linestyle='--')
axs[2,1].scatter(lai_cgls.loc[:,patch_10 >perc],vodx.loc[:,patch_10 >perc],c='black',s=2)
axs[2,1].scatter(lai_cgls_win.loc[:,patch_10 >perc].mean(axis=1),vodx_win.loc[:,patch_10 >perc].mean(axis=1),c='blue',s=4)
axs[2,1].scatter(lai_cgls_spr.loc[:,patch_10 >perc].mean(axis=1),vodx_spr.loc[:,patch_10 >perc].mean(axis=1),c='green',s=4)
axs[2,1].scatter(lai_cgls_sum.loc[:,patch_10 >perc].mean(axis=1),vodx_sum.loc[:,patch_10 >perc].mean(axis=1),c='red',s=4)
axs[2,1].scatter(lai_cgls_aut.loc[:,patch_10 >perc].mean(axis=1),vodx_aut.loc[:,patch_10 >perc].mean(axis=1),c='yellow',s=4)
axs[2,1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[2,1].set_ylabel('Matched VODX [m^2/m^-2]',color='red')

plt.tight_layout(rect=[0,0,1,0.95])

#plt.show()
plt.savefig(imgdir+'LAI_vs_{0}VODX_Patch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)
plt.close()
###################################################################
### VODC
###################################################################
#plt.rcParams['figure.figsize'] = (11,7)
plt.rcParams['figure.figsize'] = (7,11)
fig, axs = plt.subplots(3,2)
fig.suptitle('LAI vs Matched VODC\n{0}'.format(ecol),fontsize=14)
#fig.suptitle('LAI vs VODX\n{0}'.format(ecol),fontsize=14)

#fig.suptitle('LAI vs vodc\nDominant Vegetation > {0}%\n{1}'.format(round(perc*100),ecol),fontsize=14)

axs[0,0].set_title('Deciduous Forests')
axs[0,0].set_xlim(0,6.5)
axs[0,0].set_ylim(0,6.5)
xlims3=(0,6.5)
x3, y3, ext3, line3 = prep_reg(lai_cgls.loc[:,patch_3 >perc],vodc.loc[:,patch_3 >perc],xlims3)
axs[0,0].plot(ext3,line3,color='purple',linestyle='--')
axs[0,0].scatter(lai_cgls.loc[:,patch_3 >perc],vodc.loc[:,patch_3 >perc],c='black',s=2)
axs[0,0].scatter(lai_cgls_win.loc[:,patch_3 >perc].mean(axis=1),vodc_win.loc[:,patch_3 >perc].mean(axis=1),c='blue',s=4)
axs[0,0].scatter(lai_cgls_spr.loc[:,patch_3 >perc].mean(axis=1),vodc_spr.loc[:,patch_3 >perc].mean(axis=1),c='green',s=4)
axs[0,0].scatter(lai_cgls_sum.loc[:,patch_3 >perc].mean(axis=1),vodc_sum.loc[:,patch_3 >perc].mean(axis=1),c='red',s=4)
axs[0,0].scatter(lai_cgls_aut.loc[:,patch_3 >perc].mean(axis=1),vodc_aut.loc[:,patch_3 >perc].mean(axis=1),c='yellow',s=4)
axs[0,0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[0,0].set_ylabel('Matched VODC [m^2/m^-2]',color='red')

axs[0,1].set_title('Coniferous Forests')
axs[0,1].set_xlim(0,6.5)
axs[0,1].set_ylim(0,6.5)
xlims4=(0,6.5)
x4, y4, ext4, line4 = prep_reg(lai_cgls.loc[:,patch_4 >perc],vodc.loc[:,patch_4 >perc],xlims4)
axs[0,1].plot(ext4,line4,color='purple',linestyle='--')
axs[0,1].scatter(lai_cgls.loc[:,patch_4 >perc],vodc.loc[:,patch_4 >perc],c='black',s=2)
axs[0,1].scatter(lai_cgls_win.loc[:,patch_4 >perc].mean(axis=1),vodc_win.loc[:,patch_4 >perc].mean(axis=1),c='blue',s=4)
axs[0,1].scatter(lai_cgls_spr.loc[:,patch_4 >perc].mean(axis=1),vodc_spr.loc[:,patch_4 >perc].mean(axis=1),c='green',s=4)
axs[0,1].scatter(lai_cgls_sum.loc[:,patch_4 >perc].mean(axis=1),vodc_sum.loc[:,patch_4 >perc].mean(axis=1),c='red',s=4)
axs[0,1].scatter(lai_cgls_aut.loc[:,patch_4 >perc].mean(axis=1),vodc_aut.loc[:,patch_4 >perc].mean(axis=1),c='yellow',s=4)
axs[0,1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[0,1].set_ylabel('Matched VODC [m^2/m^-2]',color='red')

axs[1,0].set_title('C3 Crops')
axs[1,0].set_xlim(0,6.5)
axs[1,0].set_ylim(0,6.5)
xlims6=(0,6.5)
x6, y6, ext6, line6 = prep_reg(lai_cgls.loc[:,patch_6 >perc],vodc.loc[:,patch_6 >perc],xlims6)
axs[1,0].plot(ext6,line6,color='purple',linestyle='--')
axs[1,0].scatter(lai_cgls.loc[:,patch_6 >perc],vodc.loc[:,patch_6 >perc],c='black',s=2)
axs[1,0].scatter(lai_cgls_win.loc[:,patch_6 >perc].mean(axis=1),vodc_win.loc[:,patch_6 >perc].mean(axis=1),c='blue',s=4)
axs[1,0].scatter(lai_cgls_spr.loc[:,patch_6 >perc].mean(axis=1),vodc_spr.loc[:,patch_6 >perc].mean(axis=1),c='green',s=4)
axs[1,0].scatter(lai_cgls_sum.loc[:,patch_6 >perc].mean(axis=1),vodc_sum.loc[:,patch_6 >perc].mean(axis=1),c='red',s=4)
axs[1,0].scatter(lai_cgls_aut.loc[:,patch_6 >perc].mean(axis=1),vodc_aut.loc[:,patch_6 >perc].mean(axis=1),c='yellow',s=4)
axs[1,0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[1,0].set_ylabel('Matched VODC [m^2/m^-2]C',color='red')

axs[1,1].set_title('C4 Crops')
axs[1,1].set_xlim(0,6.5)
axs[1,1].set_ylim(0,6.5)
xlims7=(0,6.5)
x7, y7, ext7, line7 = prep_reg(lai_cgls.loc[:,patch_7 >perc],vodc.loc[:,patch_7 >perc],xlims7)
axs[1,1].plot(ext7,line7,color='purple',linestyle='--')
axs[1,1].scatter(lai_cgls.loc[:,patch_7 >perc],vodc.loc[:,patch_7 >perc],c='black',s=2)
axs[1,1].scatter(lai_cgls_win.loc[:,patch_7 >perc].mean(axis=1),vodc_win.loc[:,patch_7 >perc].mean(axis=1),c='blue',s=4)
axs[1,1].scatter(lai_cgls_spr.loc[:,patch_7 >perc].mean(axis=1),vodc_spr.loc[:,patch_7 >perc].mean(axis=1),c='green',s=4)
axs[1,1].scatter(lai_cgls_sum.loc[:,patch_7 >perc].mean(axis=1),vodc_sum.loc[:,patch_7 >perc].mean(axis=1),c='red',s=4)
axs[1,1].scatter(lai_cgls_aut.loc[:,patch_7 >perc].mean(axis=1),vodc_aut.loc[:,patch_7 >perc].mean(axis=1),c='yellow',s=4)
axs[1,1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[1,1].set_ylabel('Matched VODC [m^2/m^-2]',color='red')

axs[2,0].set_title('C3 Herbaceous')
axs[2,0].set_xlim(0,6.5)
axs[2,0].set_ylim(0,6.5)
xlims9=(0,6.5)
x9, y9, ext9, line9 = prep_reg(lai_cgls.loc[:,patch_9 >perc],vodc.loc[:,patch_9 >perc],xlims9)
axs[2,0].plot(ext9,line9,color='purple',linestyle='--')
axs[2,0].scatter(lai_cgls.loc[:,patch_9 >perc],vodc.loc[:,patch_9 >perc],c='black',s=2)
axs[2,0].scatter(lai_cgls_win.loc[:,patch_9 >perc].mean(axis=1),vodc_win.loc[:,patch_9 >perc].mean(axis=1),c='blue',s=4)
axs[2,0].scatter(lai_cgls_spr.loc[:,patch_9 >perc].mean(axis=1),vodc_spr.loc[:,patch_9 >perc].mean(axis=1),c='green',s=4)
axs[2,0].scatter(lai_cgls_sum.loc[:,patch_9 >perc].mean(axis=1),vodc_sum.loc[:,patch_9 >perc].mean(axis=1),c='red',s=4)
axs[2,0].scatter(lai_cgls_aut.loc[:,patch_9 >perc].mean(axis=1),vodc_aut.loc[:,patch_9 >perc].mean(axis=1),c='yellow',s=4)
axs[2,0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[2,0].set_ylabel('Matched VODC [m^2/m^-2]',color='red')

axs[2,1].set_title('Irrigated Crops')
axs[2,1].set_xlim(0,6.5)
axs[2,1].set_ylim(0,6.5)
xlims10=(0,6.5)
x10, y10, ext10, line10 = prep_reg(lai_cgls.loc[:,patch_10 >perc],vodc.loc[:,patch_10 >perc],xlims10)
axs[2,1].plot(ext10,line10,color='purple',linestyle='--')
axs[2,1].scatter(lai_cgls.loc[:,patch_10 >perc],vodc.loc[:,patch_10 >perc],c='black',s=2)
axs[2,1].scatter(lai_cgls_win.loc[:,patch_10 >perc].mean(axis=1),vodc_win.loc[:,patch_10 >perc].mean(axis=1),c='blue',s=4)
axs[2,1].scatter(lai_cgls_spr.loc[:,patch_10 >perc].mean(axis=1),vodc_spr.loc[:,patch_10 >perc].mean(axis=1),c='green',s=4)
axs[2,1].scatter(lai_cgls_sum.loc[:,patch_10 >perc].mean(axis=1),vodc_sum.loc[:,patch_10 >perc].mean(axis=1),c='red',s=4)
axs[2,1].scatter(lai_cgls_aut.loc[:,patch_10 >perc].mean(axis=1),vodc_aut.loc[:,patch_10 >perc].mean(axis=1),c='yellow',s=4)
axs[2,1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[2,1].set_ylabel('Matched VODC [m^2/m^-2]',color='red')

plt.tight_layout(rect=[0,0,1,0.95])
plt.savefig(imgdir+'LAI_vs_{0}VODC_Patch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)
plt.close()
###################################################################
### Scatter Plots, with seasons, only Vegetation >10%
###################################################################

perc=0.1
plt.rcParams['figure.figsize'] = (11,11)
fig, axs = plt.subplots(1,1)
axs.set_title('LAI vs Matched VODX : Vegetation > {0}%\n{1}'.format(round(perc*100),ecol),fontsize=16)
axs.set_xlim(0,6.5)
axs.set_ylim(0,1.5)
xlims=(0,6.5)
xv, yv, extv, linev = prep_reg(lai_cgls.loc[:,patch_nonveg <perc],vodx.loc[:,patch_nonveg <perc],xlims)
axs.plot(extv,linev,color='purple',linestyle='--')
axs.scatter(lai_cgls.loc[:,patch_nonveg <perc],vodx.loc[:,patch_nonveg <perc],c='black',s=2)
axs.scatter(lai_cgls_win.loc[:,patch_nonveg <perc].mean(axis=1),vodx_win.loc[:,patch_nonveg <perc].mean(axis=1),c='blue',s=4)
axs.scatter(lai_cgls_spr.loc[:,patch_nonveg <perc].mean(axis=1),vodx_spr.loc[:,patch_nonveg <perc].mean(axis=1),c='green',s=4)
axs.scatter(lai_cgls_sum.loc[:,patch_nonveg <perc].mean(axis=1),vodx_sum.loc[:,patch_nonveg <perc].mean(axis=1),c='red',s=4)
axs.scatter(lai_cgls_aut.loc[:,patch_nonveg <perc].mean(axis=1),vodx_aut.loc[:,patch_nonveg <perc].mean(axis=1),c='yellow',s=4)
axs.set_xlabel('LAI [m^2/m^-2]',color='green')
axs.set_ylabel('VODX',color='red')
#plt.show()
plt.savefig(imgdir+'LAI_vs_{0}VODX_VEGPatch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)
plt.close()

perc=0.1
plt.rcParams['figure.figsize'] = (11,11)
fig, axs = plt.subplots(1,1)
axs.set_xlim(0,6.5)
axs.set_ylim(0,1.5)
xlims=(0,6.5)
x, y, ext, line = prep_reg(lai_cgls.loc[:,patch_nonveg <perc],vodc.loc[:,patch_nonveg <perc],xlims)
axs.plot(ext,line,color='purple',linestyle='--')
axs.set_title('LAI vs Matched VODC : Vegetation > {0}%\n{1}'.format(round(perc*100),ecol),fontsize=16)
axs.scatter(lai_cgls.loc[:,patch_nonveg <perc],vodc.loc[:,patch_nonveg <perc],c='black',s=2)
axs.scatter(lai_cgls_win.loc[:,patch_nonveg <perc].mean(axis=1),vodc_win.loc[:,patch_nonveg <perc].mean(axis=1),c='blue',s=4)
axs.scatter(lai_cgls_spr.loc[:,patch_nonveg <perc].mean(axis=1),vodc_spr.loc[:,patch_nonveg <perc].mean(axis=1),c='green',s=4)
axs.scatter(lai_cgls_sum.loc[:,patch_nonveg <perc].mean(axis=1),vodc_sum.loc[:,patch_nonveg <perc].mean(axis=1),c='red',s=4)
axs.scatter(lai_cgls_aut.loc[:,patch_nonveg <perc].mean(axis=1),vodc_aut.loc[:,patch_nonveg <perc].mean(axis=1),c='yellow',s=4)
axs.set_xlabel('LAI [m^2/m^-2]',color='green')
axs.set_ylabel('VODC',color='red')

#plt.show()
plt.savefig(imgdir+'LAI_vs_{0}VODC_VEGPatch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)
plt.close()

###################################################################
### Scatter Plots, with seasons, LOW VEGETATION (<10%)
###################################################################

perc=0.1
plt.rcParams['figure.figsize'] = (11,11)
fig, axs = plt.subplots(1,1)
axs.set_title('LAI vs VODX : Vegetation < {0}%\n{1}'.format(round(perc*100),ecol),fontsize=16)
axs.set_xlim(0,6.5)
axs.set_ylim(0,1.5)
xlims=(0,6.5)
xv, yv, extv, linev = prep_reg(lai_cgls.loc[:,patch_nonveg >perc],vodx.loc[:,patch_nonveg >perc],xlims)
axs.plot(extv,linev,color='purple',linestyle='--')
axs.scatter(lai_cgls.loc[:,patch_nonveg >perc],vodx.loc[:,patch_nonveg >perc],c='black',s=2)
axs.scatter(lai_cgls_win.loc[:,patch_nonveg >perc].mean(axis=1),vodx_win.loc[:,patch_nonveg >perc].mean(axis=1),c='blue',s=4)
axs.scatter(lai_cgls_spr.loc[:,patch_nonveg >perc].mean(axis=1),vodx_spr.loc[:,patch_nonveg >perc].mean(axis=1),c='green',s=4)
axs.scatter(lai_cgls_sum.loc[:,patch_nonveg >perc].mean(axis=1),vodx_sum.loc[:,patch_nonveg >perc].mean(axis=1),c='red',s=4)
axs.scatter(lai_cgls_aut.loc[:,patch_nonveg >perc].mean(axis=1),vodx_aut.loc[:,patch_nonveg >perc].mean(axis=1),c='yellow',s=4)
axs.set_xlabel('LAI [m^2/m^-2]',color='green')
axs.set_ylabel('VODX',color='red')
#plt.show()
plt.savefig(imgdir+'LAI_vs_{0}VODX_LowVEGPatch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)
plt.close()

perc=0.1
plt.rcParams['figure.figsize'] = (11,11)
fig, axs = plt.subplots(1,1)
axs.set_xlim(0,6.5)
axs.set_ylim(0,1.5)
xlims=(0,6.5)
x, y, ext, line = prep_reg(lai_cgls.loc[:,patch_nonveg >perc],vodc.loc[:,patch_nonveg >perc],xlims)
axs.plot(ext,line,color='purple',linestyle='--')
axs.set_title('LAI vs VODC : Vegetation < {0}%\n{1}'.format(round(perc*100),ecol),fontsize=16)
axs.scatter(lai_cgls.loc[:,patch_nonveg >perc],vodc.loc[:,patch_nonveg >perc],c='black',s=2)
axs.scatter(lai_cgls_win.loc[:,patch_nonveg >perc].mean(axis=1),vodc_win.loc[:,patch_nonveg >perc].mean(axis=1),c='blue',s=4)
axs.scatter(lai_cgls_spr.loc[:,patch_nonveg >perc].mean(axis=1),vodc_spr.loc[:,patch_nonveg >perc].mean(axis=1),c='green',s=4)
axs.scatter(lai_cgls_sum.loc[:,patch_nonveg >perc].mean(axis=1),vodc_sum.loc[:,patch_nonveg >perc].mean(axis=1),c='red',s=4)
axs.scatter(lai_cgls_aut.loc[:,patch_nonveg >perc].mean(axis=1),vodc_aut.loc[:,patch_nonveg >perc].mean(axis=1),c='yellow',s=4)
axs.set_xlabel('LAI [m^2/m^-2]',color='green')
axs.set_ylabel('VODC',color='red')

#plt.show()
plt.savefig(imgdir+'LAI_vs_{0}VODC_LowVEGPatch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)
plt.close()

###################################################################
### Scatter Plots, with seasons, INDIVIDUAL Plots for each Patch
###################################################################

for ip in patch_num:
    patch_test = patch_frac[ip,:,:].reshape((280*140))
    fig, axs = plt.subplots(1,1)
    fig.suptitle('CONUS : {0}'.format(patch_name[indx]),fontsize=16)
    axs.set_title('LAI CGLS vs VODX')
    axs.scatter(lai_cgls.loc[:,patch_test >0.5],vodx.loc[:,patch_test >0.5],c='black',)
    axs.scatter(lai_cgls_win.loc[:,patch_test >0.5].mean(axis=1),vodx_win.loc[:,patch_test >0.5].mean(axis=1),c='blue')
    axs.scatter(lai_cgls_spr.loc[:,patch_test >0.5].mean(axis=1),vodx_spr.loc[:,patch_test >0.5].mean(axis=1),c='green')
    axs.scatter(lai_cgls_sum.loc[:,patch_test >0.5].mean(axis=1),vodx_sum.loc[:,patch_test >0.5].mean(axis=1),c='red')
    axs.scatter(lai_cgls_aut.loc[:,patch_test >0.5].mean(axis=1),vodx_aut.loc[:,patch_test >0.5].mean(axis=1),c='yellow')
    axs.set_xlabel('LAI [m^2/m^-2]',color='green')
    axs.set_ylabel('VODX',color='red')
    #plt.show()
    plt.savefig(imgdir+'VODX_Scatter_{0}.png'.format(patch_name[indx]),format='png',dpi=300)
    indx = indx+1

###################################################################
### STATISTICS CALCULTAIONS
###################################################################
perc = 0.5

### VODX
r_decid_x = correlation(lai_cgls,vodx,patch_3,perc)
r_decid_win_x = correlation(lai_cgls_win,vodx_win,patch_3,perc)
r_decid_spr_x = correlation(lai_cgls_spr,vodx_spr,patch_3,perc)
r_decid_sum_x = correlation(lai_cgls_sum,vodx_sum,patch_3,perc)
r_decid_aut_x = correlation(lai_cgls_aut,vodx_aut,patch_3,perc)

r_conif_x = correlation(lai_cgls,vodx,patch_4,perc)
r_conif_win_x = correlation(lai_cgls_win,vodx_win,patch_4,perc)
r_conif_spr_x = correlation(lai_cgls_spr,vodx_spr,patch_4,perc)
r_conif_sum_x = correlation(lai_cgls_sum,vodx_sum,patch_4,perc)
r_conif_aut_x = correlation(lai_cgls_aut,vodx_aut,patch_4,perc)

r_c3_x = correlation(lai_cgls,vodx,patch_6,perc)
r_c3_win_x = correlation(lai_cgls_win,vodx_win,patch_6,perc)
r_c3_spr_x = correlation(lai_cgls_spr,vodx_spr,patch_6,perc)
r_c3_sum_x = correlation(lai_cgls_sum,vodx_sum,patch_6,perc)
r_c3_aut_x = correlation(lai_cgls_aut,vodx_aut,patch_6,perc)

r_c4_x = correlation(lai_cgls,vodx,patch_7,perc)
r_c4_win_x = correlation(lai_cgls_win,vodx_win,patch_7,perc)
r_c4_spr_x = correlation(lai_cgls_spr,vodx_spr,patch_7,perc)
r_c4_sum_x = correlation(lai_cgls_sum,vodx_sum,patch_7,perc)
r_c4_aut_x = correlation(lai_cgls_aut,vodx_aut,patch_7,perc)

r_grass_x = correlation(lai_cgls,vodx,patch_9,perc)
r_grass_win_x = correlation(lai_cgls_win,vodx_win,patch_9,perc)
r_grass_spr_x = correlation(lai_cgls_spr,vodx_spr,patch_9,perc)
r_grass_sum_x = correlation(lai_cgls_sum,vodx_sum,patch_9,perc)
r_grass_aut_x = correlation(lai_cgls_aut,vodx_aut,patch_9,perc)

r_irr_x = correlation(lai_cgls,vodx,patch_10,perc)
r_irr_win_x = correlation(lai_cgls_win,vodx_win,patch_10,perc)
r_irr_spr_x = correlation(lai_cgls_spr,vodx_spr,patch_10,perc)
r_irr_sum_x = correlation(lai_cgls_sum,vodx_sum,patch_10,perc)
r_irr_aut_x = correlation(lai_cgls_aut,vodx_aut,patch_10,perc)

### VODC

r_decid_c = correlation(lai_cgls,vodc,patch_3,perc)
r_decid_win_c = correlation(lai_cgls_win,vodc_win,patch_3,perc)
r_decid_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_3,perc)
r_decid_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_3,perc)
r_decid_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_3,perc)

r_conif_c = correlation(lai_cgls,vodc,patch_4,perc)
r_conif_win_c = correlation(lai_cgls_win,vodc_win,patch_4,perc)
r_conif_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_4,perc)
r_conif_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_4,perc)
r_conif_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_4,perc)

r_c3_c = correlation(lai_cgls,vodc,patch_6,perc)
r_c3_win_c = correlation(lai_cgls_win,vodc_win,patch_6,perc)
r_c3_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_6,perc)
r_c3_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_6,perc)
r_c3_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_6,perc)

r_c4_c = correlation(lai_cgls,vodc,patch_7,perc)
r_c4_win_c = correlation(lai_cgls_win,vodc_win,patch_7,perc)
r_c4_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_7,perc)
r_c4_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_7,perc)
r_c4_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_7,perc)

r_grass_c = correlation(lai_cgls,vodc,patch_9,perc)
r_grass_win_c = correlation(lai_cgls_win,vodc_win,patch_9,perc)
r_grass_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_9,perc)
r_grass_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_9,perc)
r_grass_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_9,perc)

r_irr_c = correlation(lai_cgls,vodc,patch_10,perc)
r_irr_win_c = correlation(lai_cgls_win,vodc_win,patch_10,perc)
r_irr_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_10,perc)
r_irr_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_10,perc)
r_irr_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_10,perc)

### Vegetated Patches only, both VODX and VODC

perc = 0.1
r_veg_x = correlation(lai_cgls,vodx,patch_nonveg,perc)
r_veg_win_x = correlation(lai_cgls_win,vodx_win,patch_nonveg,perc)
r_veg_spr_x = correlation(lai_cgls_spr,vodx_spr,patch_nonveg,perc)
r_veg_sum_x = correlation(lai_cgls_sum,vodx_sum,patch_nonveg,perc)
r_veg_aut_x = correlation(lai_cgls_aut,vodx_aut,patch_nonveg,perc)

r_veg_c = correlation(lai_cgls,vodc,patch_nonveg,perc)
r_veg_win_c = correlation(lai_cgls_win,vodc_win,patch_nonveg,perc)
r_veg_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_nonveg,perc)
r_veg_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_nonveg,perc)
r_veg_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_nonveg,perc)

### Output correlations to files - so I don't need to calculate them again
### Or just print out the stats - these don't take too long


headers = ['Veg Type','All Seasons','Winter','Spring','Summer','Autumn']
veg_x = ['All Vegetation',r_veg_x[0],r_veg_win_x[0],r_veg_spr_x[0],r_veg_sum_x[0],r_veg_aut_x[0]] 
decid_x = ['Deciduous',r_decid_x[0],r_decid_win_x[0],r_decid_spr_x[0],r_decid_sum_x[0],r_decid_aut_x[0]] 
conif_x = ['Coniferous',r_conif_x[0],r_conif_win_x[0],r_conif_spr_x[0],r_conif_sum_x[0],r_conif_aut_x[0]] 
c3_x = ['C3 Crops',r_c3_x[0],r_c3_win_x[0],r_c3_spr_x[0],r_c3_sum_x[0],r_c3_aut_x[0]]
c4_x = ['C4 Crops',r_c4_x[0],r_c4_win_x[0],r_c4_spr_x[0],r_c4_sum_x[0],r_c4_aut_x[0]]
grass_x = ['C3 Herbaceous',r_grass_x[0],r_grass_win_x[0],r_grass_spr_x[0],r_grass_sum_x[0],r_grass_aut_x[0]]
irr_x = ['Irrigated Crops',r_irr_x[0],r_irr_win_x[0],r_irr_spr_x[0],r_irr_sum_x[0],r_irr_aut_x[0]]

datax = [veg_x,decid_x,conif_x,c3_x,c4_x,grass_x,irr_x]
cortable_x = pd.DataFrame(data=datax,columns=headers)

headers = ['Veg Type','All Seasons','Winter','Spring','Summer','Autumn']
veg_c = ['All Vegetation',r_veg_c[0],r_veg_win_c[0],r_veg_spr_c[0],r_veg_sum_c[0],r_veg_aut_c[0]] 
decid_c = ['Deciduous',r_decid_c[0],r_decid_win_c[0],r_decid_spr_c[0],r_decid_sum_c[0],r_decid_aut_c[0]] 
conif_c = ['Coniferous',r_conif_c[0],r_conif_win_c[0],r_conif_spr_c[0],r_conif_sum_c[0],r_conif_aut_c[0]] 
c3_c = ['C3 Crops',r_c3_c[0],r_c3_win_c[0],r_c3_spr_c[0],r_c3_sum_c[0],r_c3_aut_c[0]]
c4_c = ['C4 Crops',r_c4_c[0],r_c4_win_c[0],r_c4_spr_c[0],r_c4_sum_c[0],r_c4_aut_c[0]]
grass_c = ['C3 Herbaceous',r_grass_c[0],r_grass_win_c[0],r_grass_spr_c[0],r_grass_sum_c[0],r_grass_aut_c[0]]
irr_c = ['Irrigated Crops',r_irr_c[0],r_irr_win_c[0],r_irr_spr_c[0],r_irr_sum_c[0],r_irr_aut_c[0]]

datac = [veg_c,decid_c,conif_c,c3_c,c4_c,grass_c,irr_c]
cortable_c = pd.DataFrame(data=datac,columns=headers)

print('************************************************************\n************************************************************')
print("Matched VODX")
print(cortable_x.round(2))
print('************************************************************')
print("Matchec VODC")
print(cortable_c.round(2))

###################################################################
### Scatter Plots, with seasons, only Vegetation >10% - EXPERIMENTAL
###################################################################
import sys, glob, os, re, math, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
import pylab as pl
import scipy.optimize as optim
from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr, linregress
from numpy.polynomial.polynomial import polyfit
import seaborn as sns


#x = lai_cgls.loc[:,patch_nonveg <perc].to_numpy().flatten()
#y = vodx.loc[:,patch_nonveg <perc].to_numpy().flatten()

#sns_plot = sns.regplot(x,y)

#sns_plot.savefig(imgdir+'LAI_vs_{0}VODX_VEGPatch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)


perc=0.1
plt.rcParams['figure.figsize'] = (11,11)
fig, axs = plt.subplots(1,1)
#fig.suptitle('LAI vs VODX : Vegetation > {0}%'.format(round(perc*100)),fontsize=16)

axs.set_title('LAI vs VODX : Vegetation > {0}%\n{1}'.format(round(perc*100),ecol),fontsize=16)
#x = np.nan_to_num(lai_cgls.loc[:,patch_nonveg <perc].to_numpy().flatten())
#y = np.nan_to_num(vodx.loc[:,patch_nonveg <perc].to_numpy().flatten())

x = lai_cgls.loc[:,patch_nonveg <perc].mean(axis=1)
y = vodx.loc[:,patch_nonveg <perc].mean(axis=1)

x[np.invert(~np.isnan(y))]=np.nan
y[np.invert(~np.isnan(x))]=np.nan

x = x[~np.isnan(x)]
y = y[~np.isnan(y)]


slope, intercept, r_value, p_value, std_err = linregress(x,y)
#axs.scatter(lai_cgls.loc[:,patch_nonveg <perc],vodx.loc[:,patch_nonveg <perc],c='black',s=2)

axs.set_xlim(0,6.5)
axs.set_ylim(0,1.5)
xlims=plt.xlim()
### Attempt 1 at extending the line
'''
xlims=plt.xlim()
pd.concat([pd.Series(xlims[0]),x])
pd.concat([pd.Series(np.polyval(p1,xlims[0])),y])
x.append(pd.Series(xlims[1]))
y.append(pd.Series(np.polyval(p1,xlims[1])))
'''

### Attempt 2 at extending the regression
ext = np.arange(xlims[0],xlims[1],1)
line = [slope*xi + intercept for xi in ext]
axs.set_title('LAI vs VODX : Vegetation > {0}%\n{1}'.format(round(perc*100),ecol),fontsize=16)
axs.scatter(x,y,c='black',s=2)
#plt.plot(x,b+m*x,color='purple',linestyle='--')
plt.plot(ext,line,color='purple',linestyle='--')

#axs.scatter(lai_cgls_win.loc[:,patch_nonveg <perc].mean(axis=1),vodx_win.loc[:,patch_nonveg <perc].mean(axis=1),c='blue',s=4)
#axs.scatter(lai_cgls_spr.loc[:,patch_nonveg <perc].mean(axis=1),vodx_spr.loc[:,patch_nonveg <perc].mean(axis=1),c='green',s=4)
#axs.scatter(lai_cgls_sum.loc[:,patch_nonveg <perc].mean(axis=1),vodx_sum.loc[:,patch_nonveg <perc].mean(axis=1),c='red',s=4)
#axs.scatter(lai_cgls_aut.loc[:,patch_nonveg <perc].mean(axis=1),vodx_aut.loc[:,patch_nonveg <perc].mean(axis=1),c='yellow',s=4)

axs.set_xlabel('LAI [m^2/m^-2]',color='green')
axs.set_ylabel('VODX',color='red')
#plt.show()
plt.savefig(imgdir+'LAI_vs_{0}VODX_VEGPatch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)