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
from scipy.stats import zscore
from scipy import stats

###################################################################
### Useful Functions
###################################################################

def pearsonr_ci(x,y,alpha=0.05):
    ''' calculate Pearson correlation along with the confidence interval using scipy and numpy
    Parameters
    ----------
    x, y : iterable object such as a list or np.array
      Input for correlation calculation
    alpha : float
      Significance level. 0.05 by default
    Returns
    -------
    r : float
      Pearson's correlation coefficient
    pval : float
      The corresponding p value
    lo, hi : float
      The lower and upper bound of confidence intervals
    '''
    r, p = stats.pearsonr(x,y)
    r_z = np.arctanh(r)
    se = 1/np.sqrt(x.size-3)
    z = stats.norm.ppf(1-alpha/2)
    lo_z, hi_z = r_z-z*se, r_z+z*se
    lo, hi = np.tanh((lo_z, hi_z))
    return r, p, lo, hi

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

def ZscoreM(x):
    for ii in [1,2,3,4,5,6,7,8,9,10,11,12]:
        for name in x.columns:
            season = (x[name].index.month == ii)
            x[name].ix[season] = (x[name].ix[season]-x[name].ix[season].mean())/x[name].ix[season].std()
    return x

def ZscoreY(x):
    for ii in [1,2,3,4,5,6,7,8,9,10,11,12]:
        season = (x.index.month == ii)
        x.ix[season] = (x.ix[season]-x.ix[season].mean())/x.ix[season].std()
    return x

#def ts_std(*argv):
    #std_list = np.Series()
    #for i in len(argv[1]):
        #std_list.append(argv[1]
###################################################################
### Import Data  
###################################################################

### OL LAI
ol_lai_1 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS1/results/ol/LAI_2018-01-01_2018-12-31.PData')
ol_lai_2 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS2/results/ol/LAI_2018-01-01_2018-12-31.PData')
ol_lai_3 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS3/results/ol/LAI_2018-01-01_2018-12-31.PData')
ol_lai_4 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS4/results/ol/LAI_2018-01-01_2018-12-31.PData')


### OL SWI
ol_swi_1 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS1/results/ol/WG2_2018-01-01_2018-12-31.PData')
ol_swi_2 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS2/results/ol/WG2_2018-01-01_2018-12-31.PData')
ol_swi_3 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS3/results/ol/WG2_2018-01-01_2018-12-31.PData')
ol_swi_4 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS4/results/ol/WG2_2018-01-01_2018-12-31.PData')

### EKF LAI
ekf_lai_1 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS1/results/ekf/LAI_2018-01-01_2018-12-31.PData')
ekf_lai_2 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS2/results/ekf/LAI_2018-01-01_2018-12-31.PData')
ekf_lai_3 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS3/results/ekf/LAI_2018-01-01_2018-12-31.PData')
ekf_lai_4 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS4/results/ekf/LAI_2018-01-01_2018-12-31.PData')

### EKF SWI
ekf_swi_1 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS1/results/ekf/WG2_2018-01-01_2018-12-31.PData')
ekf_swi_2 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS2/results/ekf/WG2_2018-01-01_2018-12-31.PData')
ekf_swi_3 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS3/results/ekf/WG2_2018-01-01_2018-12-31.PData')
ekf_swi_4 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/ENS/US_ENS4/results/ekf/WG2_2018-01-01_2018-12-31.PData')

###################################################################
### Data Processing 
###################################################################

ol_lai_1 = ol_lai_1.resample('M').mean()
ol_lai_2 = ol_lai_2.resample('M').mean()
ol_lai_3 = ol_lai_3.resample('M').mean()
ol_lai_4 = ol_lai_4.resample('M').mean()

ol_swi_1 = ol_swi_1.resample('M').mean()
ol_swi_2 = ol_swi_2.resample('M').mean()
ol_swi_3 = ol_swi_3.resample('M').mean()
ol_swi_4 = ol_swi_4.resample('M').mean()

ekf_lai_1 = ekf_lai_1.resample('M').mean()
ekf_lai_2 = ekf_lai_2.resample('M').mean()
ekf_lai_3 = ekf_lai_3.resample('M').mean()
ekf_lai_4 = ekf_lai_4.resample('M').mean()

ekf_swi_1 = ekf_swi_1.resample('M').mean()
ekf_swi_2 = ekf_swi_2.resample('M').mean()
ekf_swi_3 = ekf_swi_3.resample('M').mean()
ekf_swi_4 = ekf_swi_4.resample('M').mean()

#sys.exit()
###################################################################
### Computer Z Scores or STD
###################################################################

    
    

###################################################################
### Computer Correlations
###################################################################



#sys.exit()
###################################################################
### Graphing Functions 
raw_input("Press Enter to continue ...")

'''
plt.title('LAI from ECO_II vs ECO_SG with Observed Corn Yield over Nebraska')
plt.plot(lai_Y['Yield_b/a'],label='Corn Yield - USDA',marker='o',linestyle='--',markersize=9,color='green')

plt.plot(lai_Y['Obs'],label='LAI Observations',color='black',marker='*',linestyle='--',markersize=9)
#plt.plot(co2yz['CO2'],label='CO2',color='black',linewidth =1)

plt.plot(lai_Y['Model'],label='OL ECO_II LAI',color='blue')
plt.plot(lai_Y['Analysis'],label='EKF ECO_II LAI',color='red')
plt.plot(lai_sg_olz.mean(axis=1),label='OL ECO_SG LAI',color='green')
plt.plot(lai_sg_ekfz.mean(axis=1),label='EKF ECO_SG LAI',color='orange')
plt.axhline(color='black',linewidth=.5)
plt.legend(loc='upper left')
'''
'''
plt.title('ENS 1 vs ENS 2 LAI over CONUS')
plt.plot(ol_lai_1.mean(axis=1),label='ENS1 OL LAI',color='b')
plt.plot(ol_lai_2.mean(axis=1),label='ENS2 OL LAI',color='c')

plt.plot(ekf_lai_1.mean(axis=1),label='ENS1 EKF LAI',color='r')
plt.plot(ekf_lai_2.mean(axis=1),label='ENS2 EKF LAI',color='orange')
'''

#plt.title('ENS 1 vs ENS 2 SWI over CONUS')
#plt.plot(ol_swi_1.mean(axis=1),label='ENS1 OL SWI',color='b',linestyle='--')
#plt.plot(ol_swi_2.mean(axis=1),label='ENS2 OL SWI',color='c',linestyle='--')
#plt.plot(ol_swi_3.mean(axis=1),label='ENS3 OL SWI',color='g',linestyle='--')
#plt.plot(ol_swi_4.mean(axis=1),label='ENS4 OL SWI',color='r',linestyle='--')

#plt.plot(ekf_swi_1.mean(axis=1),label='ENS1 EKF SWI',color='b')
#plt.plot(ekf_swi_2.mean(axis=1),label='ENS2 EKF SWI',color='c')
#plt.plot(ekf_swi_3.mean(axis=1),label='ENS3 EKF SWI',color='g')
#plt.plot(ekf_swi_4.mean(axis=1),label='ENS4 EKF SWI',color='r')


plt.title('ENS 1 vs ENS 2 LAI over CONUS')
plt.plot(ol_lai_1.mean(axis=1),label='ENS1 OL LAI',color='b',linestyle='--')
plt.plot(ol_lai_2.mean(axis=1),label='ENS2 OL LAI',color='c',linestyle='--')
plt.plot(ol_lai_3.mean(axis=1),label='ENS3 OL LAI',color='g',linestyle='--')
plt.plot(ol_lai_4.mean(axis=1),label='ENS4 OL LAI',color='r',linestyle='--')

plt.plot(ekf_lai_1.mean(axis=1),label='ENS1 EKF LAI',color='b')
plt.plot(ekf_lai_2.mean(axis=1),label='ENS2 EKF LAI',color='c')
plt.plot(ekf_lai_3.mean(axis=1),label='ENS3 EKF LAI',color='g')
plt.plot(ekf_lai_4.mean(axis=1),label='ENS4 EKF LAI',color='r')


plt.legend(loc='upper left')


plt.rcParams.update({'font.size': 10})
#fig.tight_layout()
#save('../vegeoFigures/Yield-RainF', ext='ps', close=True, verbose=True)
plt.show()

