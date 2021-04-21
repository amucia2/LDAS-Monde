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


numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts

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


