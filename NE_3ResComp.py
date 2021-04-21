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

###################################################################
### Import Data  
###################################################################

### Import individual yearly files
v2_025 = np.load('../Nebraska_LDAS/observations/sfx-trip/LAI_V2_2014-01-01_2014-12-31.PData')
v2_010 = np.load('../Nebraska_HR/observations/sfx-trip/LAI_V2_2014-01-01_2014-12-31.PData')
v2_001 = np.load('../Nebraska_VHR/observations/sfx-trip/LAI300_2014-01-01_2014-12-31.PData')
###################################################################
### Data Processing 
### Only takes data from the same time as observations 
###################################################################


s_ol2 = v2_025[~np.isnan(v2_025)].dropna()
s_an2 = v2_010[~np.isnan(v2_010)].dropna()
s_ob2 = v2_001

jan = s_ob2['2014-01-01':'2014-01-31'].mean(axis=0).values.reshape((400,1000))
feb = s_ob2['2014-02-01':'2014-02-28'].mean(axis=0).values.reshape((400,1000))
mar = s_ob2['2014-03-01':'2014-03-31'].mean(axis=0).values.reshape((400,1000))
apr = s_ob2['2014-04-01':'2014-04-30'].mean(axis=0).values.reshape((400,1000))
may = s_ob2['2014-05-01':'2014-05-31'].mean(axis=0).values.reshape((400,1000))
jun = s_ob2['2014-06-01':'2014-06-30'].mean(axis=0).values.reshape((400,1000))
jul = s_ob2['2014-07-01':'2014-07-31'].mean(axis=0).values.reshape((400,1000))
aug = s_ob2['2014-08-01':'2014-08-31'].mean(axis=0).values.reshape((400,1000))
sep = s_ob2['2014-09-01':'2014-09-30'].mean(axis=0).values.reshape((400,1000))
oco = s_ob2['2014-10-01':'2014-10-31'].mean(axis=0).values.reshape((400,1000))
nov = s_ob2['2014-11-01':'2014-11-30'].mean(axis=0).values.reshape((400,1000))
dec = s_ob2['2014-12-01':'2014-12-31'].mean(axis=0).values.reshape((400,1000))

#sys.exit()
###################################################################
### Graphing Functions 
### Only takes data from the same time as observations 
###################################################################
###################################################################
### Mapping Functions 
### To map data from the same time as observations  
###################################################################

### Takes data and reshapes over the grid of interest
### Needs to be changed based on data resolution and extent


v4 = s_ol2.mean(axis=0).values.reshape((16,40))
v5 = s_an2.mean(axis=0).values.reshape((40,100))
v6 = s_ob2.mean(axis=0).values.reshape((400,1000))

### [min_lon, min_lat, max_lon, max_lat]
extent = [-105,39.5,-95,43.5]

### For single image plots use the following
### For creating several subplots use the following
### Create grid of subplots

ax01 = plt.subplot2grid((3,1), (0,0))
ax02 = plt.subplot2grid((3,1), (1,0))
ax03 = plt.subplot2grid((3,1), (2,0))

ax01.set_title("LAI 25km")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax01)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(v4,interpolation='none',cmap='YlGn',vmin=0.2,vmax=2)

ax02.set_title("LAI 10km")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax02)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(v5,interpolation='none',cmap='YlGn',vmin=0.2,vmax=2)

ax03.set_title("LAI 1km")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax03)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(v6,interpolation='none',cmap='YlGn',vmin=0.2,vmax=2)
###################################################################
### Experimental colorbar functions
###################################################################
fig = plt.figure(1)
#ax = fig.add_subplot()
#cores = ax.contourf(v2[:],levels=range(4))
#cores = ax.contourf(v2)
#cbar = plt.colorbar(cores)
#cbar = plt.colorbar(ax=ax3, label='Correlation Difference',fraction=0.046, pad=0.04,orientation='horizontal')
### for 'add_axes' it is [left,bottom,witdth,height]
cbaxes = fig.add_axes([0.92,0.15,0.025,0.7])
#fig.colorbar(im,ax=axes.ravel().tolist())
fig.colorbar(im,cax=cbaxes,label='Leaf Area Index [$m^2$/$m^2$]');
###################################################################
### Final Plot
###################################################################

plt.show()

