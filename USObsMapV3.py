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

def plot_rectangle(bmap, lonmin,lonmax,latmin,latmax,color):
    xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
    ys = [latmin,latmin,latmax,latmax,latmin]
    bmap.plot(xs, ys,latlon = True,color=color,linewidth=3)

###################################################################
### This is to map the subdomains, with no data attached
### For other uses, see other versions of "USObsMap"
###################################################################

###################################################################
### Read In Data  
###################################################################

###################################################################
### Timeseries Plotting Functions  
###################################################################

###################################################################
### Mapping Functions  
###################################################################
#fig, axes = plt.subplots(2,1)
fig, axes = plt.subplots(1)
#'''
#plt.figure(figsize=(8,6),dpi=300,edgecolor='w')
axes.set_title("Domain and Subdomains",fontsize=16)
map = Basemap(llcrnrlon=-130, llcrnrlat=20, urcrnrlon=-60, urcrnrlat=55,lat_0=40,lon_0=-98,resolution='l',ax=axes)
map.drawcoastlines(1.5)
map.drawcountries(1.5)
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
### California
plot_rectangle(map,-125,-115,25,42,'black')
### SPlains
plot_rectangle(map,-110,-90,25,37,'black')
### Midwest
plot_rectangle(map,-105,-85,37,50,'black')
### NEast
plot_rectangle(map,-85,-70,37,50,'black')
### NEast
plot_rectangle(map,-105,-95,39.5,43.5,'blue')
plt.show()
#'''



axes[0].annotate("A", xy=(0.025, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')
axes[1].annotate("B", xy=(0.025, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')
plt.tight_layout()
plt.show()
