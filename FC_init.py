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


def haversine(lon1,lat1,lon2,lat2,yy_lon,yy_lat):
    """
    calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal to radians
    lon1b, lat1b, lon2b, lat2b = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2b-lon1b
    dlat = lat2b-lat1b

    a = sin(dlat/2)**2 + cos(lat1b)*cos(lat2b)*sin(dlon/2)**2
    c = 2*asin(sqrt(a))
    r = 6371. # radius of Earth in km

    return c*r,lat2,lon2,yy_lat,yy_lon,lat1,lon1

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

def NS(s,o):
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

def ubrmse(s,o):
    """
    Unbiased root Mean Squared Error
    input:
    s: simulated
    o: observed
    output: 
    ubrmses: unbiased root mean squared error
    """
    s,o = filter_nan(s,o)
    n = len(o)
    if (n!=0):
        o_moy = np.mean(o) 
        s_moy = np.mean(s)
        somme = 0.
        for i in range(n):
            somme = somme + ((s[i]-s_moy)-(o[i]-o_moy))**2
        return np.sqrt(somme/n)
    else:
        return np.nan

def bias(s, o):
    """
        Bias
        input:
        s: simulated
        o: observed
        output:
        bias: bias
        """
    s,o = filter_nan(s,o)
    return np.mean(s-o)

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts

# Load Data
evap_imp_rmsd = [33.1,33.6,29.8,27.8,23.9,20.9,19.4,18.7]
evap_deg_rmsd = [6.3,5.5,5.2,2.9,4.5,4.7,4.6,4.3]
evap_neu_rmsd = [60.6,60.9,65.0,69.3,71.6,74.5,76.0,77.0]

lai_imp_rmsd = [74.2,58.9,59.5,60.0,60.2,60.7,37.2,37.0]
lai_deg_rmsd = [0.2,0.1,0.1,0.1,0.0,0.0,0.1,0.1]
lai_neu_rmsd = [25.6,41.0,40.4,39.9,39.7,39.2,62.7,62.9]

ssm_imp_rmsd = [49.2,26.9,10.7,4.5,2.7,1.8,1.4,1.2]
ssm_deg_rmsd = [2.8,1.4,1.9,1.8,1.5,1.4,1.4,1.1]
ssm_neu_rmsd = [48.0,71.7,87.4,93.7,95.8,96.9,97.2,97.7]

evap_imp_r = [19.6,18.8,16.9,16.1,14.5,13.1,12.3,12.1]
evap_deg_r = [3.3,2.7,2.6,2.4,1.7,1.6,1.5,1.3]
evap_neu_r = [77.2,78.5,80.6,81.5,83.8,85.3,86.2,86.6]

lai_imp_r = [54.7,31.0,30.5,32.3,31.1,31.3,20.2,19.5]
lai_deg_r = [3.8,4.2,3.9,3.7,3.2,2.9,2.8,2.8]
lai_neu_r = [41.5,64.9,65.6,64.0,65.7,65.8,77.0,77.7]

ssm_imp_r = [29.0,12.1,5.7,3.6,3.4,2.6,1.7,1.7]
ssm_deg_r = [1.0,0.1,0.4,0.8,1.3,1.4,1.1,1.2]
ssm_neu_r = [70.1,87.7,94.0,95.5,95.3,96.1,97.2,97.1]


# Plot Data
#fig, ax1= plt.subplots()
fig = plt.figure(figsize=(8,2),dpi=300,edgecolor='w')
ax1 = fig.add_subplot(111)
#ax2 = ax1.twinx()

ticks = ['NO FC','FC2','FC4', 'FC6', 'FC8', 'FC10', 'FC12', 'FC14']

#dif = np.array(ekf)-np.array(ol)
ax1.set_xlim(-1,8)
#plt.title('In-Situ USCRN SM {0}m vs ISBA {1} Mean Correlation for OL and EKF Forecasts'.format(dep, response),fontsize=7)
plt.xticks(np.arange(8),ticks)

#plt.title('EKF vs OL Initialization: Evapotranspiration RMSD [Kg/m2/s] % Pixels Improved, Degraded, or Neutral',fontsize=7)
#ax1.set_ylabel('[%]')
#l1 = ax1.plot(evap_imp_rmsd,color='b',marker='h',linewidth=1,markersize=4,label='Improved')[0]
#l2 = ax1.plot(evap_deg_rmsd,color='r',marker='h',linewidth=1,markersize=4,label='Degraded')[0]
#l3 = ax1.plot(evap_neu_rmsd,color='y',marker='h',linewidth=1,markersize=4,label='Neutral')[0]

#plt.title('EKF vs OL Initialization: LAI RMSD [m2/m2] % Pixels Improved, Degraded, or Neutral',fontsize=7)
#ax1.set_ylabel('[%]')
#l1 = ax1.plot(lai_imp_rmsd,color='b',marker='h',linewidth=1,markersize=4,label='Improved')[0]
#l2 = ax1.plot(lai_deg_rmsd,color='r',marker='h',linewidth=1,markersize=4,label='Degraded')[0]
#l3 = ax1.plot(lai_neu_rmsd,color='y',marker='h',linewidth=1,markersize=4,label='Neutral')[0]

#plt.title('EKF vs OL Initialization: SSM RMSD [m3/m3] % Pixels Improved, Degraded, or Neutral',fontsize=7)
#ax1.set_ylabel('[%]')
#l1 = ax1.plot(ssm_imp_rmsd,color='b',marker='h',linewidth=1,markersize=4,label='Improved')[0]
#l2 = ax1.plot(ssm_deg_rmsd,color='r',marker='h',linewidth=1,markersize=4,label='Degraded')[0]
#l3 = ax1.plot(ssm_neu_rmsd,color='y',marker='h',linewidth=1,markersize=4,label='Neutral')[0]

#plt.title('EKF vs OL Initialization: Evapotranspiration R % Pixels Improved, Degraded, or Neutral',fontsize=7)
#ax1.set_ylabel('[%]')
#l1 = ax1.plot(evap_imp_r,color='b',marker='h',linewidth=1,markersize=4,label='Improved')[0]
#l2 = ax1.plot(evap_deg_r,color='r',marker='h',linewidth=1,markersize=4,label='Degraded')[0]
#l3 = ax1.plot(evap_neu_r,color='y',marker='h',linewidth=1,markersize=4,label='Neutral')[0]

plt.title('EKF vs OL Initialization: LAI R % Pixels Improved, Degraded, or Neutral',fontsize=7)
ax1.set_ylabel('[%]')
l1 = ax1.plot(lai_imp_r,color='b',marker='h',linewidth=1,markersize=4,label='Improved')[0]
l2 = ax1.plot(lai_deg_r,color='r',marker='h',linewidth=1,markersize=4,label='Degraded')[0]
l3 = ax1.plot(lai_neu_r,color='y',marker='h',linewidth=1,markersize=4,label='Neutral')[0]

#plt.title('EKF vs OL Initialization: SSM R % Pixels Improved, Degraded, or Neutral',fontsize=7)
#ax1.set_ylabel('[%]')
#l1 = ax1.plot(ssm_imp_r,color='b',marker='h',linewidth=1,markersize=4,label='Improved')[0]
#l2 = ax1.plot(ssm_deg_r,color='r',marker='h',linewidth=1,markersize=4,label='Degraded')[0]
#l3 = ax1.plot(ssm_neu_r,color='y',marker='h',linewidth=1,markersize=4,label='Neutral')[0]

#ax2.set_ylabel('R Diff (EKF-OL)',fontsize=8)
#l7 = ax2.plot(dif,color='black',marker='.',linewidth=0,markersize=4)[0]
#ax2.set_ylim(-0.0025,0.0025)
ax1.set_ylim(0,100)
#ax2.axhline(0,-1,8,color='black')
#ax2.tick_params(labelsize=8)
#ax2.fontsize(8)

#plt.figlegend((l1,l2,l3),('OL WG3','EKF WG3','OL WG6','EKF WG6','OL WG8','EKF WG8'),loc=(0.80,0.60),fontsize=4)
plt.legend(loc='center right',fontsize=6)
plt.rcParams.update({'font.size':8})
#plt.savefig('images/USFC_SMN_cor_difV5.png',format='png',dpi=600)
plt.show()
sys.exit()



######################################


