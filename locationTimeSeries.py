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


numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts


file = sorted(glob.glob('observations/LAI_V2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),pd.read_pickle(file[18])]
lai_cgls_ = pd.concat(df)

#date = pd.date_range(start="2003-01-01 09:00:00",end='2018-12-31 09:00:00',freq='D')
#lai_cgls = lai_cgls_.reindex(date,fill_value=np.nan)

### Take only values where LAI observations are found (both time and space)
lai_vodx_ = lai_vodx[~np.isnan(lai_cgls_)]

### Take Time series over entire period, averaged over entire domain
lai_cgls = lai_cgls_.mean(axis=1)

#lai_vodx_old = pd.read_pickle('LAIfromVODCAX_TS.PData')
lai_vodx = pd.read_pickle('LAIfromVODX_V8_TS_whereObs.PData').dropna()
lai_vodc = pd.read_pickle('LAIfromVODC_V8_TS_whereObs.PData').dropna()
lai_vodx_int = pd.read_pickle('LAIfromVODCAX_V7_2003-2018_interpolated_TSwhereObs.PData').dropna()

#lai_vodx_old = lai_vodx_old['2003-01-01':'2018-12-31']
#lai_vodx_old = lai_vodx_old[~np.isnan(lai_cgls)].dropna()
#lai_cgls = lai_cgls.dropna()

plt.title("LAI from VOD over CONUS : 2003-2018")
plt.ylabel('LAI [m^2 m^-2]')
plt.plot(lai_vodx,color='red',label='Matched VODX')
#plt.plot(lai_vodc,color='blue',label='Matched VODC')
plt.plot(lai_vodx_int,color='blue',linestyle='--',label='Matched VODX with Interpolated LAI Obs')
plt.plot(lai_cgls,color='green',label='CGLS LAI Observations')
plt.legend()
plt.show()


### Take time sereies over only 2012, averaged over entire domain
lai_vodx2012 = lai_vodx['2012-01-01':'2012-12-31']
lai_vodxold2012 = lai_vodx_old['2012-01-01':'2012-12-31']
lai_vodc2012 = lai_vodc['2012-01-01':'2012-12-31']
obs2012 = lai_cgls['2012-01-01':'2012-12-31']
lai_vodx_int2012 = lai_vodx_int['2012-01-01':'2012-12-31']

plt.title("LAI from VOD over CONUS : 2012")
plt.ylabel('LAI [m^2 m^-2]')
plt.plot(lai_vodx2012,color='red',label='Matched VODX')
#plt.plot(lai_vodc2012,color='blue',label='Matched VODC')
#plt.plot(lai_vodxold2012,color='red',linestyle='--',label='OLD Matched VODX')
plt.plot(lai_vodx_int2012,color='blue',linestyle='--',label='Matched VODX with Interpolated LAI Obs')
plt.plot(obs2012,color='green',label='CGLS LAI Observations')
plt.legend()
plt.show()

### Select individual gridcell (or subsets of gridcells) for plotting
### Will need to be done on server
### To verify I have selected the right area, I can try plotting it as a map (maybe with all other values as Nan, to keep the place in the matrix)

l1 = lai_cgls.mean(axis=0).values.reshape((140,280))
vx1 = lai_vodx.mean(axis=0).values.reshape((140,280))
vc1 = lai_vodc.mean(axis=0).values.reshape((140,280))

l1 = lai_cgls.values.reshape((140,280))
vx1 = lai_vodx.values.reshape((140,280))
vc1 = lai_vodc.values.reshape((140,280))

### Mapping
fig, axes = plt.subplots(1,1)

map = Basemap(llcrnrlon=-130, llcrnrlat=20, urcrnrlon=-60, urcrnrlat=55,lat_0=40,lon_0=-98,resolution='l',ax=axes[0])
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.5)
map.drawcountries(1.5)
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
cs = map.imshow(a1,interpolation='none',cmap='coolwarm_r',vmin=0,vmax=1)
cbar = map.colorbar(cs,location='bottom',pad="10%")
cbar.set_label("SWI [%]")
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()

plt.tight_layout()
plt.show()