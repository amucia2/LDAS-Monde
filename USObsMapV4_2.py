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

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/Panels/LAI_V2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),pd.read_pickle(file[18])]
lai_cgls = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00_V3/observations/sfx-trip/*SSM_COMBINED*'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),pd.read_pickle(file[18])]
ssm_esa = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00_V3/observations/sfx-trip/*VODX*'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),pd.read_pickle(file[18])]
lai_vod  = pd.concat(df)

#lai_vod = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/LAIfromVODCAX_V4.PData')

lat = np.loadtxt('/cnrm/vegeo/albergelc/Project/SMOS/liste_2018.txt')[:,0]
lon = np.loadtxt('/cnrm/vegeo/albergelc/Project/SMOS/liste_2018.txt')[:,1]

cgls = lai_cgls.mean(axis=0).values.reshape((140,280))
esa = ssm_esa.mean(axis=0).values.reshape((140,280))
vod = lai_vod.mean(axis=0).values.reshape((140,280))
'''
vod_ts = lai_vod.mean(axis=1)
cgls_ts = lai_cgls.mean(axis=1)
isba_ts = lai_isba.mean(axis=1)

vod_ts2 = vod_ts['2015-01-01 09:00:00':'2018-12-31 09:00:00'] 
cgls_ts2 = cgls_ts['2015-01-01 09:00:00':'2018-12-31 09:00:00']
isba_ts2 = isba_ts['2015-01-01 09:00:00':'2018-12-31 09:00:00']

vod_ts3 = vod_ts['2012-01-01 09:00:00':'2012-12-31 09:00:00']
cgls_ts3 = cgls_ts['2012-01-01 09:00:00':'2012-12-31 09:00:00']
isba_ts3 = isba_ts['2012-01-01 09:00:00':'2012-12-31 09:00:00']

vod_mean = vod_ts.resample('A').mean() 
cgls_mean = cgls_ts.resample('A').mean()
isba_mean = isba_ts.resample('A').mean()

vod_max = vod_ts.resample('A').max()
cgls_max = cgls_ts.resample('A').max()
isba_max = isba_ts.resample('A').max()
'''

###################################################################
### Timeseries Plotting Functions  
###################################################################
'''
import matplotlib.dates as mdates
years_fmt = mdates.DateFormatter('%Y')
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(years_fmt)

ax.plot(vod_ts,c='r',label='LAI from VOD')
ax.plot(cgls_ts,c='g',label='LAI Obs (CGLS)')
ax.set_ylabel("LAI [$m^2$/$m^2$]")
ax.set_title('LAI Observations vs LAI from VOD') 
plt.legend()
plt.show()

fig, ax = plt.subplots()
ax.plot(vod_ts2,c='r',label='LAI from VOD')
ax.plot(cgls_ts2,c='g',label='LAI Obs (CGLS)')
ax.set_ylabel("LAI [$m^2$/$m^2$]")
ax.set_title('LAI Observations vs LAI from VOD') 
ax.set_ylim(top=3.5)
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.locator_params(axis='x', nbins=5)
plt.legend()
plt.show()

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.plot(vod_ts3,c='r',label='LAI from VOD')
ax.plot(cgls_ts3,c='g',label='LAI Obs (CGLS)')
ax.set_ylabel("LAI [$m^2$/$m^2$]")
ax.set_title('LAI Observations vs LAI from VOD')
ax.set_ylim(top=3)
ax.xaxis.set_major_locator(plt.MaxNLocator(12))
ax.locator_params(axis='x', nbins=13)
plt.legend()
plt.show()


import matplotlib.lines as mlines
fig, ax = plt.subplots()

ax.set_title('Annual Mean LAI Obs vs Matched VOD Obs')
ax.set_xlim(1,1.45)
ax.set_ylim(1,1.45)
ax.set_xlabel('Matched VOD Obs')
ax.set_ylabel('LAI Obs')
line = mlines.Line2D([0, 1], [0, 1], color='red')
ax.add_line(line)
transform = ax.transAxes
line.set_transform(transform)
ax.scatter(vod_mean,cgls_mean)
plt.show()

import matplotlib.lines as mlines
fig, ax = plt.subplots()
ax.set_title('Annual Maximum LAI Obs vs Matched VOD Obs')
ax.set_xlim(1,6)
ax.set_ylim(1,6)
ax.set_xlabel('Matched VOD Obs')
ax.set_ylabel('LAI Obs')
line = mlines.Line2D([0, 1], [0, 1], color='red')
ax.add_line(line)
transform = ax.transAxes
line.set_transform(transform)
ax.scatter(vod_max,cgls_max)
plt.show()
'''
###################################################################
### Mapping Functions  
###################################################################
import matplotlib.gridspec as gridspec
imgdir = '/cnrm/vegeo/muciaa/Images/US00/Observations_Comparisons/'

#fig, axes = plt.subplots(2,1)
#fig, axes = plt.subplots(1)
#'''
#plt.rcParams['figure.figsize'] = (11,8)
fig=plt.figure(constrained_layout=True)
gs = gridspec.GridSpec(ncols=4, nrows=4,figure=fig)


ax1 = fig.add_subplot(gs[:2,:2])

ax1.set_title("CGLS LAI V2",fontsize=16)
map = Basemap(llcrnrlon=-130, llcrnrlat=20, urcrnrlon=-60, urcrnrlat=55,lat_0=40,lon_0=-98,resolution='l',ax=ax1)
#map.shadedrelief()
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.5)
map.drawcountries(1.5)
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
cs1 = map.imshow(cgls,interpolation='none',cmap='RdYlGn',vmin=0,vmax=4)
cbar1 = map.colorbar(cs1,location='bottom',pad="10%")
cbar1.set_label("LAI [$m^2$/$m^2$]")
tick_locator = ticker.MaxNLocator(nbins=4)
cbar1.locator = tick_locator
cbar1.update_ticks()
#plt.show()
#'''
ax2 = fig.add_subplot(gs[:2,2:])
ax2.set_title("ESA CCI SM v04.5 COMBINED",fontsize=16)
map = Basemap(llcrnrlon=-130, llcrnrlat=20, urcrnrlon=-60, urcrnrlat=55,lat_0=40,lon_0=-98,resolution='l',ax=ax2)
#map.shadedrelief()
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.5)
map.drawcountries(1.5)
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),linewidth=.5)
cs2 = map.imshow(esa,interpolation='none',cmap='RdYlGn',vmin=0,vmax=0.5)
cbar2 = map.colorbar(cs2,location='bottom',pad="10%")
cbar2.set_label("SSM [$m^3$/$m^3$]")
tick_locator = ticker.MaxNLocator(nbins=4)
cbar2.locator = tick_locator
cbar2.update_ticks()

ax3 = fig.add_subplot(gs[2:4,1:3])
ax3.set_title("VODCA X-Band VOD",fontsize=16)
map = Basemap(llcrnrlon=-130, llcrnrlat=20, urcrnrlon=-60, urcrnrlat=55,lat_0=40,lon_0=-98,resolution='l',ax=ax3)
#map.shadedrelief()
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.5)
map.drawcountries(1.5)
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
cs2 = map.imshow(vod,interpolation='none',cmap='RdYlGn',vmin=0,vmax=1)
cbar2 = map.colorbar(cs2,location='bottom',pad="10%")
cbar2.set_label("Optical Depth")
tick_locator = ticker.MaxNLocator(nbins=4)
cbar2.locator = tick_locator
cbar2.update_ticks()
'''
plt.sca(axes[1])
x1,y1 = map(-105,39.5)
x2,y2 = map(-95,39.5)
x3,y3 = map(-95,43.5)
x4,y4 = map(-105,43.5)
poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],edgecolor='black',fill=False,linewidth=1)
plt.gca().add_patch(poly)
'''

ax1.annotate("A", xy=(0.025, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')
ax2.annotate("B", xy=(0.025, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')
ax3.annotate("C", xy=(0.025, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')
#plt.tight_layout(rect=[0,0,1,0.95])
plt.tight_layout()
#plt.show()
plt.savefig(imgdir+'LAI_SSM_VOD_Comparison.png',format='png',dpi=300)



###################################################################
### Mapping Functions  
###################################################################


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