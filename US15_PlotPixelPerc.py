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
import scipy.io as spio


def save(path, ext='png', close=True, verbose=True):
    """Save a figure from pyplot.
    Parameters
    ----------
    path : string
    The path (and filename, without the extension) to save the
    figure to.
    ext : string (default='png')
    The file extension. This must be supported by the active
    matplotlib backend (see matplotlib.backends module).  Most
    backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.
    close : boolean (default=True)
    Whether to close the figure after saving.  If you want to save
    the figure multiple times (e.g., to multiple formats), you
    should NOT close it in between saves or you will have to
    re-plot it.
    verbose : boolean (default=True)
    whether to print information about when and where the image
    has been saved.
    """
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
       directory = '.'
    #If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)
    # The final path to save to
    savepath = os.path.join(directory, filename)
    if verbose:
        print("Saving figure to '%s'..." % savepath),
    # Actually save the figure
    plt.savefig(savepath)
    # Close it
    if close:
        plt.close()
    if verbose:
        print("Done")

def filter_nan(s,o):
    """
    this functions removed the data  from simulated and observed data
    whereever the observed data contains nan
    this is used by all other functions, otherwise they will produce nan as 
    output
    """
    data = np.array([s.flatten(),o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]
    #data = data[~np.isnan(data)]
    return data[:,0],data[:,1]
    """
    Nash Sutcliffe efficiency coefficient
    input:
    s: simulated
    o: observed
    output:
    ns: Nash Sutcliffe efficient coefficient
    """
    s,o = filter_nan(s,o)
    return 1 - sum((s-o)**2)/sum((o-np.mean(o))**2)
    #return 1 - sum((np.log(s)-np.log(o))**2)/sum((np.log(o)-np.mean(np.log(o)))**2)

def rmse(s,o):
    """
    Root Mean Squared Error
    input:
    s: simulated
    o: observed
    output: 
    rmses: root mean squared error
    """
    s,o = filter_nan(s,o)
    return np.sqrt(np.mean((s-o)**2))

def correlation(s,o):
    """
    correlation coefficient
    input:
    s: simulated
    o: observed
    output:
    correlation: correlation coefficient
    """
    s,o = filter_nan(s,o)
    if s.size == 0:
        corr = np.NaN
    else:
        #corr = np.corrcoef(o, s)[0,1]	
        corr = pearsonr(o, s)					        
    return corr

numbers = re.compile(r'(\d+)')

def numericalSort(value):
    parts = numbers.split(value)
    return parts

eva_R_i = 14,13,12,11,11,11,10,10,10,10,9,9,9,8
eva_R_d = 4,4,4,3,3,3,2,2,2,2,2,2,1,1
eva_R_n = 81,83,84,85,86,87,87,88,88,89,89,90,90,91

eva_RMSD_i = 7,7,6,6,5,5,5,4,4,3,3,3,3,2
eva_RMSD_d = 4,3,3,3,2,2,2,2,1,1,1,1,1,1
eva_RMSD_n = 90,90,91,91,92,93,93,94,95,95,96,96,97,97

lai_R_i = 33,33,32,32,31,31,31,34,34,33,23,21,20,20
lai_R_d =4,4,4,4,4,3,3,3,3,3,3,3,3,3
lai_R_n = 63,63,65,65,65,65,65,63,63,65,74,76,77,77

lai_RMSD_i = 56,56,56,56,56,57,56,57,57,57,43,38,35,35
lai_RMSD_d = 0,0,0,0,0,0,0,0,0,0,0,0,0,0
lai_RMSD_n = 43,43,44,44,44,43,43,43,43,43,57,62,65,65

ssm_R_i = 27,18,13,10,9,8,7,5,5,4,4,4,4,4
ssm_R_d = 1,0,1,1,1,1,2,2,3,3,3,3,3,3
ssm_R_n = 72,82,86,89,90,91,92,93,92,93,93,93,92,93

ssm_RMSD_i = 47,34,24,17,13,10,8,6,6,5,5,5,6,5
ssm_RMSD_d = 1,2,2,2,2,2,2,2,1,1,1,1,1,1
ssm_RMSD_n = 52,64,74,81,85,88,90,92,93,93,94,94,93,94

### Create Graphs ###
fig = plt.figure(figsize=(8,2.5),dpi=300,edgecolor='w')
ax1 = fig.add_subplot(111)
ticks = ['24','48','72','96','120','144','168','192','216','240','264','288','312','336']
plt.xticks(np.arange(15),ticks,fontsize=4)

ax1.set_xlabel('fc Lead Time (Hours)',fontsize=5)
ax1.set_xlim(-1,14)

import matplotlib.transforms as transforms
offset = transforms.ScaledTranslation(0.05, 0,fig.dpi_scale_trans)
trans = ax1.transData-offset

plt.axvline(x=0.5,color='black',linestyle='--')
#plt.axvline(x=4.5,color='black',linestyle='--')

ax1.set_ylim(0,100)
ax1.set_ylabel('R [% Change]')
### EVAP
#ax1.set_ylim(0.70,0.75)
plt.title('Evapotranspiration : Impact of Initialization on Forecast R',fontsize=7)
l1 = ax1.plot(eva_R_i,color='b',marker='.',linewidth=1,markersize=4,label='Improved')[0]
l2 = ax1.plot(eva_R_d,color='r',marker='.',linewidth=1,markersize=4,label='Degraded')[0]
l3 = ax1.plot(eva_R_n,color='orange',marker='.',linewidth=1,markersize=4,label='Neutral')[0]
#### SSM
#plt.title('SSM : Impact of Initialization on Forecast R',fontsize=7)
#l1 = ax1.plot(ssm_R_i,color='b',marker='.',linewidth=1,markersize=4,label='Improved')[0]
#l2 = ax1.plot(ssm_R_d,color='r',marker='.',linewidth=1,markersize=4,label='Degraded')[0]
#l3 = ax1.plot(ssm_R_n,color='orange',marker='.',linewidth=1,markersize=4,label='Neutral')[0]
#### LAI
#plt.title('LAI : Impact of Initialization on Forecast R',fontsize=7)
#l1 = ax1.plot(lai_R_i,color='b',marker='.',linewidth=1,markersize=4,label='Improved')[0]
#l2 = ax1.plot(lai_R_d,color='r',marker='.',linewidth=1,markersize=4,label='Degraded')[0]
#l3 = ax1.plot(lai_R_n,color='orange',marker='.',linewidth=1,markersize=4,label='Neutral')[0]
plt.show()

fig = plt.figure(figsize=(8,2.5),dpi=300,edgecolor='w')
ax1 = fig.add_subplot(111)
ticks = ['24','48','72','96','120','144','168','192','216','240','264','288','312','336']
plt.xticks(np.arange(15),ticks,fontsize=4)

ax1.set_xlabel('fc Lead Time (Hours)',fontsize=5)
ax1.set_xlim(-1,14)

import matplotlib.transforms as transforms
offset = transforms.ScaledTranslation(0.05, 0,fig.dpi_scale_trans)
trans = ax1.transData-offset
ax1.set_ylim(0,100)
plt.axvline(x=0.5,color='black',linestyle='--')
#plt.axvline(x=4.5,color='black',linestyle='--')
ax1.set_ylabel('RMSD [% Change]')
### EVAP
#ax1.set_ylim(0.70,0.75)
plt.title('Evapotranspiration : Impact of Initialization on Forecast RMSD',fontsize=7)
l1 = ax1.plot(eva_RMSD_i,color='b',marker='.',linewidth=1,markersize=4,label='Improved')[0]
l2 = ax1.plot(eva_RMSD_d,color='r',marker='.',linewidth=1,markersize=4,label='Degraded')[0]
l3 = ax1.plot(eva_RMSD_n,color='orange',marker='.',linewidth=1,markersize=4,label='Neutral')[0]
### SSM
#ax1.set_ylim(0.0,1)
#plt.title('SSM : Impact of Initialization on Forecast RMSD',fontsize=7)
#l1 = ax1.plot(ssm_RMSD_i,color='b',marker='.',linewidth=1,markersize=4,label='Improved')[0]
#l2 = ax1.plot(ssm_RMSD_d,color='r',marker='.',linewidth=1,markersize=4,label='Degraded')[0]
#l3 = ax1.plot(ssm_RMSD_n,color='orange',marker='.',linewidth=1,markersize=4,label='Neutral')[0]
##### LAI
#plt.title('LAI : Impact of Initialization on Forecast RMSD',fontsize=7)
#l1 = ax1.plot(lai_RMSD_i,color='b',marker='.',linewidth=1,markersize=4,label='Improved')[0]
#l2 = ax1.plot(lai_RMSD_d,color='r',marker='.',linewidth=1,markersize=4,label='Degraded')[0]
#l3 = ax1.plot(lai_RMSD_n,color='orange',marker='.',linewidth=1,markersize=4,label='Neutral')[0]
plt.show()
sys.exit()
