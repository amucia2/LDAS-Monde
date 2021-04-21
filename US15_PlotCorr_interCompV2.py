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

lai = 0.998, 0.997, 0.996, 0.995, 0.994, 0.993, 0.993, 0.992, 0.991, 0.99, 0.988, 0.986, 0.984
rzsm = 0.992, 0.991, 0.99, 0.989, 0.987, 0.985, 0.983, 0.981, 0.979, 0.976, 0.974, 0.971, 0.969
wg2 = 0.974, 0.964, 0.953, 0.939, 0.922, 0.903, 0.883, 0.862, 0.843, 0.825, 0.81, 0.799, 0.789
wg3 = 0.985, 0.978, 0.97, 0.96, 0.949, 0.936, 0.922, 0.907, 0.893, 0.878, 0.867, 0.856, 0.847
runoff = 0.893, 0.835, 0.772, 0.697, 0.618, 0.549, 0.481, 0.403, 0.347, 0.306, 0.307, 0.286, 0.266
drain = 0.994, 0.99, 0.983, 0.977, 0.968, 0.959, 0.95, 0.938, 0.925, 0.915, 0.906, 0.897, 0.886
evap = 0.964, 0.949, 0.933, 0.915, 0.896, 0.877, 0.862, 0.848, 0.836, 0.825, 0.817, 0.805, 0.801



laie    = 0.067, 0.08, 0.091, 0.101, 0.109, 0.117, 0.124, 0.13, 0.137, 0.149, 0.162, 0.173, 0.184
rzsme   = 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.013, 0.014, 0.015, 0.016, 0.016, 0.017
wg2e    = 0.02, 0.024, 0.028, 0.032, 0.036, 0.04, 0.044, 0.048, 0.051, 0.054, 0.056, 0.058, 0.059
wg3e    = 0.015, 0.018, 0.021, 0.025, 0.028, 0.032, 0.035, 0.038, 0.041, 0.044, 0.046, 0.048, 0.049
runoffe = 0.884, 1.067, 1.238, 1.414, 1.571, 1.702, 1.818, 1.939, 2.011, 2.075, 2.067, 2.088, 2.133
draine  = 0.029, 0.045, 0.065, 0.082, 0.1, 0.119, 0.135, 0.152, 0.168, 0.182, 0.195, 0.207, 0.219
evape   = 0.41, 0.473, 0.534, 0.593, 0.651, 0.702, 0.74, 0.777, 0.805, 0.832, 0.85, 0.875, 0.885

### Create Graphs ###
#fig = plt.figure(figsize=(8,2.5),dpi=300,edgecolor='w')
fig, axes = plt.subplots(2,1,sharex = True)
#ticks = ['FC2','FC3','FC4','FC5','FC6','FC7','FC8','FC9','FC10','FC11','FC12','FC13','FC14']
ticks = ['2','3','4','5','6','7','8','9','10','11','12','13','14']
plt.xticks(np.arange(15),ticks,fontsize=8)

axes[1].set_xlabel('fc Lead Time (Days)')
axes[0].set_xlim(-1,13)
axes[1].set_xlim(-1,13)

axes[0].set_ylim(0.0,1.05)
#ax1.set_yticks((0.0,0.2,0.4,0.6,0.8,1.0))
axes[0].set_yticks((0.0,0.25,0.5,0.75,1.0))
axes[0].set_ylabel('Correlation')
### EVAP
#ax1.set_ylim(0.70,0.75)
axes[0].set_title('Correlation: SEKF vs SEKF Initialized Forecast')
l1 = axes[0].plot(lai,linewidth=1,marker='.',markersize=10,label='LAI')[0]
l2 = axes[0].plot(rzsm,linewidth=1,marker='.',markersize=10,label='RZSM')[0]
l3 = axes[0].plot(wg2,linewidth=1,marker='.',markersize=10,label='WG2')[0]
#l4 = axes[0].plot(wg3,linewidth=1,marker='.',markersize=10,label='WG3')[0]
l5 = axes[0].plot(runoff,linewidth=1,marker='.',markersize=10,label='RUNOFF')[0]
l6 = axes[0].plot(drain,linewidth=1,marker='.',markersize=10,label='DRAIN')[0]
l7 = axes[0].plot(evap,linewidth=1,marker='.',markersize=10,label='EVAP')[0]
#plt.legend(loc='best',fontsize=10,numpoints=1)
#plt.show()

#### RMSD
le, re, we, rue, de, ee = ([] for i in range(6))
for i in range(len(laie)):
    le.append(laie[i]*100)
    re.append(rzsme[i]*100)
    we.append(wg2e[i]*100)
    rue.append(runoffe[i]*100)
    de.append(draine[i]*100)
    ee.append(evape[i]*100)


#ax1.set_ylim(-0.2,2.5)
axes[1].set_ylim(-20,250)
#ax1.set_ylabels(
#ax1.set_yticks((0.2,0.4,0.6,0.8,1.0))
axes[1].set_ylabel('RMSD Change [%]')
### EVAP
#ax1.set_ylim(0.70,0.75)
axes[1].set_title('RMSD Change : SEKF vs SEKF Initialized Forecast')
r1 = axes[1].plot(le,linewidth=1,marker='.',markersize=10,label='LAI')[0]
r2 = axes[1].plot(re,linewidth=1,marker='.',markersize=10,label='RZSM')[0]
r3 = axes[1].plot(we,linewidth=1,marker='.',markersize=10,label='WG2')[0]
#r4 = axes[1].plot(wg3e,linewidth=1,marker='.',markersize=10,label='WG3')[0]
r5 = axes[1].plot(rue,linewidth=1,marker='.',markersize=10,label='RUNOFF')[0]
r6 = axes[1].plot(de,linewidth=1,marker='.',markersize=10,label='DRAIN')[0]
r7 = axes[1].plot(ee,linewidth=1,marker='.',markersize=10,label='EVAP')[0]
#plt.legend(loc='best',fontsize=4,numpoints=1)

axes[0].annotate("A", xy=(0.025, 0.9), xycoords="axes fraction",fontsize=14,weight='bold')
axes[1].annotate("B", xy=(0.025, 0.9), xycoords="axes fraction",fontsize=14,weight='bold')

handles, labels = axes[0].get_legend_handles_labels()
axes[0].legend(handles,labels,loc='best',fontsize=10,numpoints=1)
plt.tight_layout()
plt.show()


sys.exit()
