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

###################################################################
### Specify Data  
###################################################################

exp = 'Nebraska_LDAS/observations/sfx-trip/'
lai = np.load('Nebraska_LDAS/results/ol/LAI_2017-01-01_2017-12-31.PData')
lai_an = np.load('Nebraska_LDAS/results/ekf/LAI_2017-01-01_2017-12-31.PData')
lai_obs = np.load('Nebraska_LDAS/observations/sfx-trip/LAI_V2_2017-01-01_2017-12-31.PData')
lai1 = np.load('Nebraska_LDAS/results/ol/LAI_2012-01-01_2012-12-31.PData')
lai1_an = np.load('Nebraska_LDAS/results/ekf/LAI_2012-01-01_2012-12-31.PData')
lai1_obs = np.load('Nebraska_LDAS/observations/sfx-trip/LAI_V2_2012-01-01_2012-12-31.PData')

###################################################################
### Mapping Funcitons
###################################################################
'''
plt.plot(lai.mean(axis=1)) ; plt.plot(lai_an.mean(axis=1)) ; plt.plot(lai_obs.mean(axis=1),marker='*',linewidth=0) ; plt.show()
v = lai.mean(axis=0).values.reshape((16,40))
ax = plt.gca()
m.imshow(v[::-1],interpolation='none')
plt.show()
'''
###################################################################
### Experimentation Functions
###################################################################

###################################################################
### Data Processing 
### Only takes data from the same time as observations 
###################################################################

tmp = lai.copy()*np.nan 
tmp.ix[lai_obs.index] = lai_obs 

tmp1 = lai1.copy()*np.nan
tmp1.ix[lai1_obs.index] = lai1_obs

### Use the following for graphing time-series
d_ol1 = lai[~np.isnan(tmp)].mean(axis=1).dropna()
d_an1 = lai_an[~np.isnan(tmp)].mean(axis=1).dropna()
d_ob1 = lai_obs[~np.isnan(tmp)].mean(axis=1).dropna()

d_ol2 = lai1[~np.isnan(tmp1)].mean(axis=1).dropna()
d_an2 = lai1_an[~np.isnan(tmp1)].mean(axis=1).dropna()
d_ob2 = lai1_obs[~np.isnan(tmp1)].mean(axis=1).dropna()

### Use the following for mapping yearly results 
s_ol1 = lai[~np.isnan(tmp)].dropna()
s_an1 = lai_an[~np.isnan(tmp)].dropna()
s_ob1 = lai_obs[~np.isnan(tmp)].dropna()

s_ol2 = lai1[~np.isnan(tmp1)].dropna()
s_an2 = lai1_an[~np.isnan(tmp1)].dropna()
s_ob2 = lai1_obs[~np.isnan(tmp1)].dropna()

###################################################################
### Graphing Functions 
### Only takes data from the same time as observations 
###################################################################

f, axarr = plt.subplots(2, sharey = True)

axarr[0].plot(d_ob1, label = 'Observation') ; axarr[0].plot(d_an1, label = 'Analysis') ; axarr[0].plot(d_ol1, label = 'Open-Loop')
axarr[0].set_title('2017 LAI')
axarr[0].legend(loc ='upper left')

axarr[1].plot(d_ob2, label = 'Observation') ; axarr[1].plot(d_an2, label = 'Analysis') ; axarr[1].plot(d_ol2, label = 'Open-Loop')
axarr[1].set_title('2012 LAI')
axarr[1].legend(loc = 'upper left')
plt.show()
raw_input("Press Enter to continue ...")

sys.exit()

###################################################################
### Mapping Functions 
### To map data from the same time as observations  
###################################################################

### Takes data and reshapes over the grid of interest
### Needs to be changed based on data resolution and extent
v1 = s_ol1.mean(axis=0).values.reshape((16,40))
v2 = s_an1.mean(axis=0).values.reshape((16,40))
v3 = s_ob1.mean(axis=0).values.reshape((16,40))

v4 = s_ol2.mean(axis=0).values.reshape((16,40))
v5 = s_an2.mean(axis=0).values.reshape((16,40))
v6 = s_ob2.mean(axis=0).values.reshape((16,40))

### Create grid of subplots
ax1 = plt.subplot2grid((2,3), (0,0))
ax2 = plt.subplot2grid((2,3), (0,1))
ax3 = plt.subplot2grid((2,3), (0,2))

ax4 = plt.subplot2grid((2,3), (1,0))
ax5 = plt.subplot2grid((2,3), (1,1))
ax6 = plt.subplot2grid((2,3), (1,2))

### Creating individual maps and assigning them to a subplot
ax1.set_title("LAI - Open-Loop - 2017")
map1 = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax1)
map1.drawstates(linewidth=1.8)
map1.drawcoastlines(1.8)
map1.drawmapboundary()
map1.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map1.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map1.drawrivers()
map1.imshow(v1,interpolation='none',cmap='YlGn')

ax2.set_title("LAI - Analysis - 2017")
map2 = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax2)
map2.drawstates(linewidth=1.8)
map2.drawcoastlines(1.8)
map2.drawmapboundary()
map2.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map2.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map2.drawrivers()
map2.imshow(v2,interpolation='none',cmap='YlGn')

ax3.set_title("LAI - Observations - 2017")
map3 = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax3)
map3.drawstates(linewidth=1.8)
map3.drawcoastlines(1.8)
map3.drawmapboundary()
map3.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map3.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map3.drawrivers()
map3.imshow(v3,interpolation='none',cmap='YlGn')

ax4.set_title("LAI - Open-Loop - 2000")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax4)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map.drawrivers()
map.imshow(v4,interpolation='none',cmap='YlGn')

ax5.set_title("LAI - Analysis - 2000")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax5)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map.drawrivers()
map.imshow(v5,interpolation='none',cmap='YlGn')

ax6.set_title("LAI - Observations - 2000")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax6)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map.drawrivers()
map.imshow(v6,interpolation='none',cmap='YlGn')

plt.imshow(v3,cmap='YlGn')

#fig = plt.figure()
#cbaxes = fig.add_axes([0.25,0.5,0.75,0.1])
divider = make_axes_locatable(ax2)
cax = divider.append_axes("bottom", size = "10%",pad = .1)
#cax = divider.append_axes("bottom", size = "10%", pad=0.05)
plt.colorbar(ax=ax2, cax=cax,label='Leaf Area Index [$m^2$/$m^2$]',orientation='horizontal')

##plt.colorbar(ax=ax2, label='Leaf Area Index [$m^2$/$m^2$]',fraction=0.046, pad=0.04,orientation='horizontal')

plt.show()

