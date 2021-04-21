import sys, glob, os, re, math, time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
import pylab as pl
import scipy.optimize as optim
from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr
#from mpl_toolkits.axes_grid1 import make_axes_locatable

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts

#data = readAllData(obs, analysis_period, assim_hour, patch_frac, out_all_dir,post_from_pickle, mod_ana_dir, ['Model','Analysis','Obs'],IPNT=IPNT, map_set=map_set, filter_nan=not('SIF' in obs or 'E_ACTUAL' in obs or 'FLGPP' in obs))

###################################################################
### FOR RUNNING IN LDAS POST (PDB) 
###################################################################

### Run ldasPost with pdb.set_trace() after postObs

### Model = LAI ISBA (ol or ekf)
### Analysis = VOD (X or C) renamed to LAI in the correct folder
### Obs = LAI CGLS
datap = data.to_pandas()

lai_isba = datap['Model']
vodx = datap['Analysis']
lai_cgls = datap['Obs']
vodc = datap['Analysis']

lai_isba.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/NE/LAI_ISBA.Pdata')
vodx.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/NE/VODX.Pdata')
lai_cgls.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/NE/LAI_CGLS.Pdata')
vodc.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/NE/VODC.Pdata')

lai_isba.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/NEast/LAI_ISBA.Pdata')
vodx.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/NEast/VODX.Pdata')
lai_cgls.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/NEast/LAI_CGLS.Pdata')
vodc.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/NEast/VODC.Pdata')

lai_isba.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/California/LAI_ISBA.Pdata')
vodx.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/California/VODX.Pdata')
lai_cgls.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/California/LAI_CGLS.Pdata')
vodc.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/California/VODC.Pdata')

lai_isba.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/SPlains/LAI_ISBA.Pdata')
vodx.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/SPlains/VODX.Pdata')
lai_cgls.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/SPlains/LAI_CGLS.Pdata')
vodc.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/SPlains/VODC.Pdata')

lai_isba.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/Midwest/LAI_ISBA.Pdata')
vodx.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/Midwest/VODX.Pdata')
lai_cgls.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/Midwest/LAI_CGLS.Pdata')
vodc.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/Midwest/VODC.Pdata')

lai_isba.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/Nebraska/LAI_ISBA.Pdata')
vodx.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/Nebraska/VODX.Pdata')
lai_cgls.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/Nebraska/LAI_CGLS.Pdata')
vodc.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/Nebraska/VODC.Pdata')

lai_isba.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/CONUS/LAI_ISBA.Pdata')
vodx.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/CONUS/VODX.Pdata')
lai_cgls.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/CONUS/LAI_CGLS.Pdata')
vodc.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/CONUS/VODC.Pdata')

###################################################################
### FOR RUNNING SEPERATELY TO PLOT
###################################################################

lai_isba = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/NE/LAI_ISBA.Pdata')
vodx = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/NE/VODX.Pdata')
lai_cgls = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/NE/LAI_CGLS.Pdata')

#C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\NE

lai_isba = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\NE\LAI_ISBA.Pdata')
vodx = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\NE\VODX.Pdata')
vodc = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\NE\VODC.Pdata')
lai_cgls = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\NE\LAI_CGLS.Pdata')

lai_isba = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\NEast\LAI_ISBA.Pdata')
vodx = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\NEast\VODX.Pdata')
vodc = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\NEast\VODC.Pdata')
lai_cgls = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\NEast\LAI_CGLS.Pdata')

lai_isba = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\California\LAI_ISBA.Pdata')
vodx = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\California\VODX.Pdata')
vodc = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\California\VODC.Pdata')
lai_cgls = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\California\LAI_CGLS.Pdata')

lai_isba = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\SPlains\LAI_ISBA.Pdata')
vodx = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\SPlains\VODX.Pdata')
vodc = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\SPlains\VODC.Pdata')
lai_cgls = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\SPlains\LAI_CGLS.Pdata')

lai_isba = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\Midwest\LAI_ISBA.Pdata')
vodx = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\Midwest\VODX.Pdata')
vodc = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\Midwest\VODC.Pdata')
lai_cgls = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\Midwest\LAI_CGLS.Pdata')

lai_isba = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\Nebraska\LAI_ISBA.Pdata')
vodx = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\Nebraska\VODX.Pdata')
vodc = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\Nebraska\VODC.Pdata')
lai_cgls = pd.read_pickle(r'C:\Users\antho\Documents\SXVGO1\ldas_curr\US00\observations\Zones\Nebraska\LAI_CGLS.Pdata')

lai_cgls_12 = lai_cgls['2012-01-01':'2012-12-31']
vodx_12 = vodx['2012-01-01':'2012-12-31']
vodc_12 = vodc['2012-01-01':'2012-12-31']

lai_cgls_10 = lai_cgls['2010-01-01':'2010-12-31']
vodx_10 = vodx['2010-01-01':'2010-12-31']
vodc_10 = vodc['2010-01-01':'2010-12-31']


lai_cgls_09 = lai_cgls['2009-01-01':'2009-12-31']
vodx_09 = vodx['2009-01-01':'2009-12-31']
vodc_09 = vodc['2009-01-01':'2009-12-31']

###################################################################
### All Years Time Series
###################################################################
fig, ax1 = plt.subplots()

#plt.title('LAI vs VOD over Nebraska')
plt.title('LAI vs VOD over SPlains Region')
color='tab:green'
ax1.set_ylabel('LAI [m^2/m^-2]',color=color)
ax1.set_xlabel('Date',color='black')
ax1.plot(lai_cgls.mean(axis=1),color='black',label='LAI CGLS')

ax2 = ax1.twinx()
color='tab:red'
ax2.set_ylabel('VOD',color=color)
ax2.plot(vodx.mean(axis=1),color='green',label='VODX')
ax2.plot(vodc.mean(axis=1),color='red',label='VODC')

#plt.legend()
fig.legend()
fig.tight_layout()
plt.show()

###################################################################
### One year Time Series
###################################################################
fig, ax1 = plt.subplots()

#plt.title('LAI vs VOD over Nebraska - 2012')
plt.title('LAI vs VOD over SPlains - 2012')
color='tab:green'
ax1.set_ylabel('LAI [m^2/m^-2]',color=color)
ax1.set_xlabel('Date',color='black')
ax1.plot(lai_cgls_12.mean(axis=1),color='black',label='LAI CGLS')

ax2 = ax1.twinx()
color='tab:red'
ax2.set_ylabel('VOD',color=color)
ax2.plot(vodx_12.mean(axis=1),color='green',label='VODX')
ax2.plot(vodc_12.mean(axis=1),color='red',label='VODC')

#plt.legend()
fig.legend()
fig.tight_layout()
plt.show()

fig, ax1 = plt.subplots()
#plt.title('LAI vs VOD over Nebraska - 2009')
plt.title('LAI vs VOD over SPlains - 2009')
color='tab:green'
ax1.set_ylabel('LAI [m^2/m^-2]',color=color)
ax1.set_xlabel('Date',color='black')
ax1.plot(lai_cgls_09.mean(axis=1),color='black',label='LAI CGLS')

ax2 = ax1.twinx()
color='tab:red'
ax2.set_ylabel('VOD',color=color)
ax2.plot(vodx_09.mean(axis=1),color='green',label='VODX')
ax2.plot(vodc_09.mean(axis=1),color='red',label='VODC')

#plt.legend()
fig.legend()
fig.tight_layout()
plt.show()

###################################################################
### Scatter Plots
###################################################################


plt.title('LAI CGLS versus VODX')
plt.scatter(lai_cgls.mean(axis=1),vodx.mean(axis=1),c='black')
plt.xlabel('LAI [m^2/m^-2]',color='green')
plt.ylabel('VODX',color='red')
plt.show()

lai_cgls.index.difference(vodx.index)
lai_isba.index.difference(lai_cgls.index)

idx = pd.DatetimeIndex(lai_cgls.index)
vodc = vodc.reindex(idx, fill_value=np.nan)

idx = pd.DatetimeIndex(lai_isba.index)
lai_cgls = lai_cgls.reindex(idx, fill_value=np.nan)

plt.title('LAI CGLS versus VODC')
plt.scatter(lai_cgls.mean(axis=1),vodc.mean(axis=1),c='black')
plt.xlabel('LAI [m^2/m^-2]',color='green')
plt.ylabel('VODC',color='red')
plt.show()



plt.rcParams['figure.figsize'] = (15.5,5.5)
fig, axs = plt.subplots(1,2)
axs[0].set_title('LAI CGLS vs VODX')
axs[0].scatter(lai_cgls.mean(axis=1),vodx.mean(axis=1),c='black')
axs[0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[0].set_ylabel('VODX',color='red')

axs[1].set_title('LAI ISBA vs VODX')
axs[1].scatter(lai_isba.mean(axis=1),vodx.mean(axis=1),c='black')
axs[1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[1].set_ylabel('VODX',color='red')
plt.show()


idx = pd.DatetimeIndex(lai_cgls.index)
vodc = vodc.reindex(idx, fill_value=np.nan)

fig, axs = plt.subplots(1,2)
axs[0].set_title('LAI CGLS vs VODC')
axs[0].scatter(lai_cgls.mean(axis=1),vodc.mean(axis=1),c='black')
axs[0].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[0].set_ylabel('VODC',color='red')

axs[1].set_title('LAI ISBA vs VODC')
axs[1].scatter(lai_isba.mean(axis=1),vodc.mean(axis=1),c='black')
axs[1].set_xlabel('LAI [m^2/m^-2]',color='green')
axs[1].set_ylabel('VODC',color='red')
plt.show()
