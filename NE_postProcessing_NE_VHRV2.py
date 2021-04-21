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
#lai_025 = np.load('../Nebraska_LDAS/observations/sfx-trip/LAI_MC.t08_2014-01-01_2014-12-31.PData')
#lai_010 = np.load('../Nebraska_HR/observations/sfx-trip/LAI_MC.t08_2014-01-01_2014-12-31.PData')
#lai_001 = np.load('../Nebraska_VHR/observations/sfx-trip/LAI300_2014-01-01_2014-12-31.PData')
'''
### Observations files contain up to the present, subtract the necessary number of observations here
lai_obs = lai_obs[:-4]
'''
###################################################################
### Data Processing 
### Only takes data from the same time as observations 
###################################################################

### Use the following for mapping yearly results 
#s_ol1 = lai_025
#s_an1 = lai_010
#s_ob1 = lai_001

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
'''
### Single Graph Plot
plt.plot(d_ob1, label = 'Observation'); plt.plot(d_an1, label = 'Analysis') ; plt.plot(d_ol1, label = 'Open-Loop')
plt.legend(loc='upper left')
plt.title('2017 LAI')

### Multiple Graph Plots

f, axarr = plt.subplots(2, sharey = True)

axarr[0].plot(d_ob1, label = 'Observation'); axarr[0].plot(d_an1, label = 'Analysis') ; axarr[0].plot(d_ol1, label = 'Open-Loop')
axarr[0].set_title('LAI')
axarr[0].legend(loc ='upper left')

axarr[1].plot(d_ob2, label = 'Observation') ; axarr[1].plot(d_an2, label = 'Analysis') ; axarr[1].plot(d_ol2, label = 'Open-Loop')
axarr[1].set_title('2012 LAI')
axarr[1].legend(loc = 'upper left')

plt.show()
raw_input("Press Enter to continue ...")

sys.exit()
'''
###################################################################
### Mapping Functions 
### To map data from the same time as observations  
###################################################################

### Takes data and reshapes over the grid of interest
### Needs to be changed based on data resolution and extent

#v1 = s_ol1.mean(axis=0).values.reshape((16,40))
#v2 = s_an1.mean(axis=0).values.reshape((40,100))
#v3 = s_ob1.mean(axis=0).values.reshape((400,1000))

v4 = s_ol2.mean(axis=0).values.reshape((16,40))
v5 = s_an2.mean(axis=0).values.reshape((40,100))
v6 = s_ob2.mean(axis=0).values.reshape((400,1000))

### [min_lon, min_lat, max_lon, max_lat]
extent = [-105,39.5,-95,43.5]

### For single image plots use the following
'''
plt.title("LAI - Correlation - 2000-2018")
map1 = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h')
map1.drawstates(linewidth=1.8)
map1.drawcoastlines(1.8)
map1.drawmapboundary()
map1.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map1.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map1.drawrivers()
map1.imshow(v,interpolation='none',cmap='YlGn',extent=[0,1,0,1],aspect='auto')
plt.show()
'''
### For creating several subplots use the following
### Create grid of subplots

ax01 = plt.subplot2grid((3,1), (0,0))
ax02 = plt.subplot2grid((3,1), (1,0))
ax03 = plt.subplot2grid((3,1), (2,0))

'''
ax1 = plt.subplot2grid((3,4), (0,0))
ax2 = plt.subplot2grid((3,4), (0,1))
ax3 = plt.subplot2grid((3,4), (0,2))
ax4 = plt.subplot2grid((3,4), (0,3))

ax5 = plt.subplot2grid((3,4), (1,0))
ax6 = plt.subplot2grid((3,4), (1,1))
ax7 = plt.subplot2grid((3,4), (1,2))
ax8 = plt.subplot2grid((3,4), (1,3))

ax9 = plt.subplot2grid((3,4), (2,0))
ax10 = plt.subplot2grid((3,4), (2,1))
ax11 = plt.subplot2grid((3,4), (2,2))
ax12 = plt.subplot2grid((3,4), (2,3))
'''
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
'''

### Creating individual maps and assigning them to a subplot
ax1.set_title("January")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax1)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(jan,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)

ax2.set_title("February")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax2)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(feb,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)

ax3.set_title("March")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax3)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(mar,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)

ax4.set_title("April")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax4)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(apr,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)

ax5.set_title("May")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax5)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(may,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)

ax6.set_title("June")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax6)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(jun,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)

ax7.set_title("July")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax7)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(jul,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)

ax8.set_title("August")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax8)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(aug,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)

ax9.set_title("September")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax9)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(sep,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)

ax10.set_title("October")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax10)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(oco,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)

ax11.set_title("November")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax11)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(nov,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)

ax12.set_title("December")
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='h',ax = ax12)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(40,43,2),labels=(1,0,0,1))
map.drawrivers()
im = map.imshow(dec,interpolation='none',cmap='YlGn',vmin=0.2,vmax=5)



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
'''
plt.show()

