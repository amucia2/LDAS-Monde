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
from matplotlib import ticker
from matplotlib.patches import Polygon


numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    return parts

###################################################################
### Read In Data  
###################################################################

file = sorted(glob.glob(r'/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC_IFS/observations/sfx-trip/raw_SWI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ascat = pd.concat(df)

file = sorted(glob.glob(r'/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC_IFS/observations/sfx-trip/LAI_V2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
lai = pd.concat(df)

lat = np.loadtxt('/cnrm/vegeo/albergelc/Project/SMOS/liste_2018.txt')[:,0]
lon = np.loadtxt('/cnrm/vegeo/albergelc/Project/SMOS/liste_2018.txt')[:,1]

a1 = ascat.mean(axis=0).values.reshape((350,700))
l1 = lai.mean(axis=0).values.reshape((350,700))

###################################################################
### Mapping Functions  
###################################################################
fig, axes = plt.subplots(1,2)
#fig, axes = plt.subplots(1)
#'''
#plt.figure(figsize=(8,6),dpi=300,edgecolor='w')
axes[0].set_title("CGLS SWI (ASCAT)",fontsize=16)
map = Basemap(llcrnrlon=-130, llcrnrlat=20, urcrnrlon=-60, urcrnrlat=55,lat_0=40,lon_0=-98,resolution='l',ax=axes[0])
#map.shadedrelief()
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.5)
map.drawcountries(1.5)
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
x,y = map(np.array(lon), np.array(lat))
cs2 = map.scatter(x,y,marker='s',s=15,color='black')
cs2 = map.scatter(x,y,marker='s',s=10,color='grey')
cs = map.imshow(a1,interpolation='none',cmap='coolwarm_r',vmin=0,vmax=1)
cbar = map.colorbar(cs,location='bottom',pad="10%")
cbar.set_label("SWI [%]")
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
#plt.show()
#'''

#plt.figure(figsize=(8,6),dpi=300,edgecolor='w')
axes[1].set_title("CGLS LAI (GEOV2)",fontsize=16)
map = Basemap(llcrnrlon=-130, llcrnrlat=20, urcrnrlon=-60, urcrnrlat=55,lat_0=40,lon_0=-98,resolution='l',ax=axes[1])
#map.shadedrelief()
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.5)
map.drawcountries(1.5)
map.drawstates(0.5)
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)

plt.sca(axes[1])
x1,y1 = map(-105,39.5)
x2,y2 = map(-95,39.5)
x3,y3 = map(-95,43.5)
x4,y4 = map(-105,43.5)
poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],edgecolor='black',fill=False,linewidth=1)
plt.gca().add_patch(poly)

#draw_screen_poly( lats, lons, map1)
cs1 = map.imshow(l1,interpolation='none',cmap='RdYlGn',vmin=0,vmax=4)
cbar2 = map.colorbar(cs1,location='bottom',pad="10%")
cbar2.set_label("LAI [$m^2$/$m^2$]")
tick_locator = ticker.MaxNLocator(nbins=4)
cbar2.locator = tick_locator
cbar2.update_ticks()


axes[0].annotate("A", xy=(0.025, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')
axes[1].annotate("B", xy=(0.025, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')
plt.tight_layout()
plt.show()


'''

lon = np.linspace(-130,-60)
lat = np.linspace(35,35)


fig, axes = plt.subplots(1) 
map = Basemap(llcrnrlon=-130, llcrnrlat=20, urcrnrlon=-60, urcrnrlat=55,lat_0=40,lon_0=-98,resolution='l',ax=axes) 
map.drawcoastlines(1.5) 
map.drawcountries(1.5) 
#map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5) 
#map.drawparallels(np.arange(35,35),labels=[1,0,0,0],linewidth=.5) 
x,y = map(lon,lat)
map.plot(x,y, linewidth=2, color='red',linestyle='--')

plt.show() 


