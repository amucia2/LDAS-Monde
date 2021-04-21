# Figures for paper LDAS-Monde EnKF for HESS
# author: B. Bonan
# Last modified: March 2021

## patchFrac_figures_Tony.py
import sys, glob, os, re, math, time
os.chdir('/cnrm/vegeo/muciaa/ldas_chain_python/LDAS_v1.3.0')
from ldasPost import *
from ldasMapSet import *
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
from scipy.stats.stats import pearsonr
from numpy.polynomial.polynomial import polyfit
os.chdir('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing')



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
vodx = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2/CONUS/VODX.Pdata')
vodc = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2/CONUS/VODC.Pdata')
lai_cgls = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/CONUS/LAI_CGLS.Pdata')

### Matched Obs
'''
vodx = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2/CONUS/LAIfromVODCAX_V8_2003-2018.Pdata')
vodc = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2/CONUS/LAIfromVODCAC_V8_2003-2018.Pdata')
'''

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

patch_3 = patch_frac[3,:,:].reshape((280*140))
patch_4 = patch_frac[4,:,:].reshape((280*140))
patch_6 = patch_frac[6,:,:].reshape((280*140))
patch_9 = patch_frac[9,:,:].reshape((280*140))

###################################################################
###################################################################
### Scatter Plots, with seasons and 4 Patch types
###################################################################
###################################################################

imgdir = '/cnrm/vegeo/muciaa/Images/US00/VOD_Comparisons/Patch/Experimental/'
perc = 0.5
eco = "ECOSG"
ecol = "ECOCLIMAP SG"
match = ""

###################################################################
### VODX
###################################################################

plt.rcParams['figure.figsize'] = (11,11)
fig, axs = plt.subplots(2,2)
fig.suptitle('LAI vs VODX : Dominant Vegetation > {0}%\n{1}'.format(round(perc*100),ecol),fontsize=16)

axs[0,0].set_title('Deciduous Forests')
axs[0,0].scatter(lai_cgls.loc[:,patch_3 >perc],vodx.loc[:,patch_3 >perc],c='black',s=2)
axs[0,0].scatter(lai_cgls_win.loc[:,patch_3 >perc].mean(axis=1),vodx_win.loc[:,patch_3 >perc].mean(axis=1),c='blue',s=4)
axs[0,0].scatter(lai_cgls_spr.loc[:,patch_3 >perc].mean(axis=1),vodx_spr.loc[:,patch_3 >perc].mean(axis=1),c='green',s=4)
axs[0,0].scatter(lai_cgls_sum.loc[:,patch_3 >perc].mean(axis=1),vodx_sum.loc[:,patch_3 >perc].mean(axis=1),c='red',s=4)
axs[0,0].scatter(lai_cgls_aut.loc[:,patch_3 >perc].mean(axis=1),vodx_aut.loc[:,patch_3 >perc].mean(axis=1),c='yellow',s=4)
axs[0,0].set_xlim(0,6.5)
axs[0,0].set_ylim(0,1.5)
axs[0,0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[0,0].set_ylabel('VODX',color='red')

axs[0,1].set_title('Coniferous Forests')
axs[0,1].scatter(lai_cgls.loc[:,patch_4 >perc],vodx.loc[:,patch_4 >perc],c='black',s=2)
axs[0,1].scatter(lai_cgls_win.loc[:,patch_4 >perc].mean(axis=1),vodx_win.loc[:,patch_4 >perc].mean(axis=1),c='blue',s=4)
axs[0,1].scatter(lai_cgls_spr.loc[:,patch_4 >perc].mean(axis=1),vodx_spr.loc[:,patch_4 >perc].mean(axis=1),c='green',s=4)
axs[0,1].scatter(lai_cgls_sum.loc[:,patch_4 >perc].mean(axis=1),vodx_sum.loc[:,patch_4 >perc].mean(axis=1),c='red',s=4)
axs[0,1].scatter(lai_cgls_aut.loc[:,patch_4 >perc].mean(axis=1),vodx_aut.loc[:,patch_4 >perc].mean(axis=1),c='yellow',s=4)
axs[0,1].set_xlim(0,6.5)
axs[0,1].set_ylim(0,1.5)
axs[0,1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[0,1].set_ylabel('VODX',color='red')

axs[1,0].set_title('C3 Crops')
axs[1,0].scatter(lai_cgls.loc[:,patch_6 >perc],vodx.loc[:,patch_6 >perc],c='black',s=2)
axs[1,0].scatter(lai_cgls_win.loc[:,patch_6 >perc].mean(axis=1),vodx_win.loc[:,patch_6 >perc].mean(axis=1),c='blue',s=4)
axs[1,0].scatter(lai_cgls_spr.loc[:,patch_6 >perc].mean(axis=1),vodx_spr.loc[:,patch_6 >perc].mean(axis=1),c='green',s=4)
axs[1,0].scatter(lai_cgls_sum.loc[:,patch_6 >perc].mean(axis=1),vodx_sum.loc[:,patch_6 >perc].mean(axis=1),c='red',s=4)
axs[1,0].scatter(lai_cgls_aut.loc[:,patch_6 >perc].mean(axis=1),vodx_aut.loc[:,patch_6 >perc].mean(axis=1),c='yellow',s=4)
axs[1,0].set_xlim(0,6.5)
axs[1,0].set_ylim(0,1.0)
axs[1,0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[1,0].set_ylabel('VODX',color='red')

axs[1,1].set_title('Grasslands')
axs[1,1].scatter(lai_cgls.loc[:,patch_9 >perc],vodx.loc[:,patch_9 >perc],c='black',s=2)
axs[1,1].scatter(lai_cgls_win.loc[:,patch_9 >perc].mean(axis=1),vodx_win.loc[:,patch_9 >perc].mean(axis=1),c='blue',s=4)
axs[1,1].scatter(lai_cgls_spr.loc[:,patch_9 >perc].mean(axis=1),vodx_spr.loc[:,patch_9 >perc].mean(axis=1),c='green',s=4)
axs[1,1].scatter(lai_cgls_sum.loc[:,patch_9 >perc].mean(axis=1),vodx_sum.loc[:,patch_9 >perc].mean(axis=1),c='red',s=4)
axs[1,1].scatter(lai_cgls_aut.loc[:,patch_9 >perc].mean(axis=1),vodx_aut.loc[:,patch_9 >perc].mean(axis=1),c='yellow',s=4)
axs[1,1].set_xlim(0,6.5)
axs[1,1].set_ylim(0,1.0)
axs[1,1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[1,1].set_ylabel('VODX',color='red')

#plt.show()
plt.savefig(imgdir+'LAI_vs_{0}VODX_Patch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)

###################################################################
### VODC
###################################################################

plt.rcParams['figure.figsize'] = (11,11)
fig, axs = plt.subplots(2,2)
fig.suptitle('LAI vs VODC : Dominant Vegetation > {0}%\n{1}'.format(round(perc*100),ecol),fontsize=16)

axs[0,0].set_title('Deciduous Forests')
axs[0,0].scatter(lai_cgls.loc[:,patch_3 >perc],vodc.loc[:,patch_3 >perc],c='black',s=2)
axs[0,0].scatter(lai_cgls_win.loc[:,patch_3 >perc].mean(axis=1),vodc_win.loc[:,patch_3 >perc].mean(axis=1),c='blue',s=4)
axs[0,0].scatter(lai_cgls_spr.loc[:,patch_3 >perc].mean(axis=1),vodc_spr.loc[:,patch_3 >perc].mean(axis=1),c='green',s=4)
axs[0,0].scatter(lai_cgls_sum.loc[:,patch_3 >perc].mean(axis=1),vodc_sum.loc[:,patch_3 >perc].mean(axis=1),c='red',s=4)
axs[0,0].scatter(lai_cgls_aut.loc[:,patch_3 >perc].mean(axis=1),vodc_aut.loc[:,patch_3 >perc].mean(axis=1),c='yellow',s=4)
axs[0,0].set_xlim(0,6.5)
axs[0,0].set_ylim(0,1.5)
axs[0,0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[0,0].set_ylabel('VODC',color='red')

axs[0,1].set_title('Coniferous Forests')
axs[0,1].scatter(lai_cgls.loc[:,patch_4 >perc],vodc.loc[:,patch_4 >perc],c='black',s=2)
axs[0,1].scatter(lai_cgls_win.loc[:,patch_4 >perc].mean(axis=1),vodc_win.loc[:,patch_4 >perc].mean(axis=1),c='blue',s=4)
axs[0,1].scatter(lai_cgls_spr.loc[:,patch_4 >perc].mean(axis=1),vodc_spr.loc[:,patch_4 >perc].mean(axis=1),c='green',s=4)
axs[0,1].scatter(lai_cgls_sum.loc[:,patch_4 >perc].mean(axis=1),vodc_sum.loc[:,patch_4 >perc].mean(axis=1),c='red',s=4)
axs[0,1].scatter(lai_cgls_aut.loc[:,patch_4 >perc].mean(axis=1),vodc_aut.loc[:,patch_4 >perc].mean(axis=1),c='yellow',s=4)
axs[0,1].set_xlim(0,6.5)
axs[0,1].set_ylim(0,1.5)
axs[0,1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[0,1].set_ylabel('VODC',color='red')

axs[1,0].set_title('C3 Crops')
axs[1,0].scatter(lai_cgls.loc[:,patch_6 >perc],vodc.loc[:,patch_6 >perc],c='black',s=2)
axs[1,0].scatter(lai_cgls_win.loc[:,patch_6 >perc].mean(axis=1),vodc_win.loc[:,patch_6 >perc].mean(axis=1),c='blue',s=4)
axs[1,0].scatter(lai_cgls_spr.loc[:,patch_6 >perc].mean(axis=1),vodc_spr.loc[:,patch_6 >perc].mean(axis=1),c='green',s=4)
axs[1,0].scatter(lai_cgls_sum.loc[:,patch_6 >perc].mean(axis=1),vodc_sum.loc[:,patch_6 >perc].mean(axis=1),c='red',s=4)
axs[1,0].scatter(lai_cgls_aut.loc[:,patch_6 >perc].mean(axis=1),vodc_aut.loc[:,patch_6 >perc].mean(axis=1),c='yellow',s=4)
axs[1,0].set_xlim(0,6.5)
axs[1,0].set_ylim(0,1.5)
axs[1,0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[1,0].set_ylabel('VODC',color='red')

axs[1,1].set_title('Grasslands')
axs[1,1].scatter(lai_cgls.loc[:,patch_9 >perc],vodc.loc[:,patch_9 >perc],c='black',s=2)
axs[1,1].scatter(lai_cgls_win.loc[:,patch_9 >perc].mean(axis=1),vodc_win.loc[:,patch_9 >perc].mean(axis=1),c='blue',s=4)
axs[1,1].scatter(lai_cgls_spr.loc[:,patch_9 >perc].mean(axis=1),vodc_spr.loc[:,patch_9 >perc].mean(axis=1),c='green',s=4)
axs[1,1].scatter(lai_cgls_sum.loc[:,patch_9 >perc].mean(axis=1),vodc_sum.loc[:,patch_9 >perc].mean(axis=1),c='red',s=4)
axs[1,1].scatter(lai_cgls_aut.loc[:,patch_9 >perc].mean(axis=1),vodc_aut.loc[:,patch_9 >perc].mean(axis=1),c='yellow',s=4)
axs[1,1].set_xlim(0,6.5)
axs[1,1].set_ylim(0,1.5)
axs[1,1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[1,1].set_ylabel('VODC',color='red')

#plt.show()
plt.savefig(imgdir+'LAI_vs_{0}VODC_Patch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)

###################################################################
### Scatter Plots, with seasons, only Vegetation >10%
###################################################################

perc=0.1
plt.rcParams['figure.figsize'] = (11,11)
fig, axs = plt.subplots(1,1)
#fig.suptitle('LAI vs VODX : Vegetation > {0}%'.format(round(perc*100)),fontsize=16)

axs.set_title('LAI vs VODX : Vegetation > {0}%\n{1}'.format(round(perc*100),ecol),fontsize=16)
axs.scatter(lai_cgls.loc[:,patch_nonveg <perc],vodx.loc[:,patch_nonveg <perc],c='black',s=2)
axs.scatter(lai_cgls_win.loc[:,patch_nonveg <perc].mean(axis=1),vodx_win.loc[:,patch_nonveg <perc].mean(axis=1),c='blue',s=4)
axs.scatter(lai_cgls_spr.loc[:,patch_nonveg <perc].mean(axis=1),vodx_spr.loc[:,patch_nonveg <perc].mean(axis=1),c='green',s=4)
axs.scatter(lai_cgls_sum.loc[:,patch_nonveg <perc].mean(axis=1),vodx_sum.loc[:,patch_nonveg <perc].mean(axis=1),c='red',s=4)
axs.scatter(lai_cgls_aut.loc[:,patch_nonveg <perc].mean(axis=1),vodx_aut.loc[:,patch_nonveg <perc].mean(axis=1),c='yellow',s=4)
axs.set_xlim(0,6.5)
axs.set_ylim(0,1.5)
axs.set_xlabel('LAI [m^2/m^-2]',color='green')
axs.set_ylabel('VODX',color='red')
#plt.show()
plt.savefig(imgdir+'LAI_vs_{0}VODX_VEGPatch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)


fig, axs = plt.subplots(1,1)
#fig.suptitle('LAI vs VODC : Vegetation > {0}%'.format(round(perc*100)),fontsize=16)
axs.set_title('LAI vs VODC : Vegetation > {0}%\n{1}'.format(round(perc*100),ecol),fontsize=16)
axs.scatter(lai_cgls.loc[:,patch_nonveg <perc],vodc.loc[:,patch_nonveg <perc],c='black',s=2)
axs.scatter(lai_cgls_win.loc[:,patch_nonveg <perc].mean(axis=1),vodc_win.loc[:,patch_nonveg <perc].mean(axis=1),c='blue',s=4)
axs.scatter(lai_cgls_spr.loc[:,patch_nonveg <perc].mean(axis=1),vodc_spr.loc[:,patch_nonveg <perc].mean(axis=1),c='green',s=4)
axs.scatter(lai_cgls_sum.loc[:,patch_nonveg <perc].mean(axis=1),vodc_sum.loc[:,patch_nonveg <perc].mean(axis=1),c='red',s=4)
axs.scatter(lai_cgls_aut.loc[:,patch_nonveg <perc].mean(axis=1),vodc_aut.loc[:,patch_nonveg <perc].mean(axis=1),c='yellow',s=4)
axs.set_xlim(0,6.5)
axs.set_ylim(0,1.5)
axs.set_xlabel('LAI [m^2/m^-2]',color='green')
axs.set_ylabel('VODC',color='red')

#plt.show()
plt.savefig(imgdir+'LAI_vs_{0}VODC_VEGPatch_{1}perc_{2}.png'.format(match,round(perc*100),eco),format='png',dpi=300)


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