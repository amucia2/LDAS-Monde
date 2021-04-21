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

swi = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/observations/obs_test/SWI/SWI_2017-01-01_2017-12-31.PData')

swi_nc = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/observations/obs_test/SWI_nc/SWI_nc_2017-01-01_2017-12-31.PData')

laiV1 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/observations/2017_2018_Test/LAI_V1_2017-01-01_2017-12-31.PData')

laiV2 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/observations/2017_2018_Test/LAI_V2_2017-01-01_2017-12-31.PData')

###################################################################
### Data Processing
###################################################################

#SWI = swi.mean(axis=0).values.reshape((40,100))
#SWI_nc = swi_nc.mean(axis=0).values.reshape((40,100))
#WI_V2 = swi_v2.mean(axis=0).values.reshape((40,100))
#awSWI_nc = rawswi_v2.mean(axis=0).values.reshape((40,100))

###################################################################
### Graphing Functions
###################################################################

plt.plot(swi.mean(axis=1),label='SWI',color='red')
plt.plot(swi_nc.mean(axis=1),label='SWI_nc',color='blue')
#plt.plot(swi_nc.mean(axis=1),label='SWI_nc')
plt.legend(loc='upper left')
plt.show()

sys.exit()

plt.plot(laiV1.mean(axis=1),label='LAI_V1',color='red')
plt.plot(laiV2.mean(axis=1),label='LAI_V2',color='blue')
plt.legend(loc='upper left')
plt.show()
sys.exit()


###################################################################
### Mapping Functions
### To map data from the same time as observations
###################################################################

extent = [-105,39.5,-95,43.5]
'''
plt.title('SWI')
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2],urcrnrlat=extent[3],resolution='h')
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map.drawrivers()
map.imshow(SWI,interpolation='none',cmap='YlGn',extent=[0,1,0,1],aspect='auto')
plt.show()
'''
### Multiple Subplots ###

fig, axs = plt.subplots(2,1)

axs[0].set_title('SWI')
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2],urcrnrlat=extent[3],resolution='h',ax=axs[0])
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map.drawrivers()
map.imshow(SWI_V2,interpolation='none',cmap='YlGn',extent=[0,1,0,1])

axs[1].set_title('SWI_nc')
map = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2],urcrnrlat=extent[3],resolution='h',ax=axs[1])
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map.drawrivers()
map.imshow(SWI_nc,interpolation='none',cmap='YlGn',extent=[0,1,0,1])



plt.show()




