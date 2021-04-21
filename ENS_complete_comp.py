###################################################################
### Import Python Libraries
###################################################################
import sys, glob, os, re, math, time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
import pylab as pl
import scipy.optimize as optim
import scipy as sp
from sklearn.metrics import mean_squared_error
from math import sqrt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import zscore
from scipy import stats
from tqdm import tqdm


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

### Import entire folders of yearly data
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    return parts

def rmse(y_actual, y_predicted):
    return sqrt(mean_squared_error(y_actual,y_predicted))


###################################################################
### Import Data  
###################################################################
        
#rain_1 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens1/results/ol/RAINF_ISBA_2018-01-01_2018-12-31.PData')
#rain_2 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens2/results/ol/RAINF_ISBA_2018-01-01_2018-12-31.PData')
#rain_3 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens3/results/ol/RAINF_ISBA_2018-01-01_2018-12-31.PData')
#rain_4 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens4/results/ol/RAINF_ISBA_2018-01-01_2018-12-31.PData')
#rain_5 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens5/results/ol/RAINF_ISBA_2018-01-01_2018-12-31.PData')
#rain_6 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens6/results/ol/RAINF_ISBA_2018-01-01_2018-12-31.PData')
#rain_7 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens7/results/ol/RAINF_ISBA_2018-01-01_2018-12-31.PData')
#rain_8 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens8/results/ol/RAINF_ISBA_2018-01-01_2018-12-31.PData')
#rain_9 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens9/results/ol/RAINF_ISBA_2018-01-01_2018-12-31.PData')
#rain_10 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens10/results/ol/RAINF_ISBA_2018-01-01_2018-12-31.PData')

numbers = re.compile(r'(\d+)')                                                                 
def numericalSort(value):                         
    parts = numbers.split(value)
    return parts


rain={}
lai={}
swi={}


#rain_era = pd.read_pickle('/cnrm/vegeo/albergelc/ldas_chain_python/LDAS_v1.2.3/MIDW_025/')
lai_era = pd.read_pickle('/cnrm/vegeo/albergelc/ldas_chain_python/LDAS_v1.2.3/MIDW_025/model/LAI_2018-01-01_2018-12-31.PData')
swi_era = pd.read_pickle('/cnrm/vegeo/albergelc/ldas_chain_python/LDAS_v1.2.3/MIDW_025/model/WG2_2018-01-01_2018-12-31.PData')

lai_obs = pd.read_pickle('/cnrm/vegeo/albergelc/ldas_chain_python/LDAS_v1.2.3/MIDW_025/observations_e5/sfx-trip/LAI_V1_2018-01-01_2018-12-31.PData')
lai_obsV2 = pd.read_pickle('/cnrm/vegeo/albergelc/ldas_chain_python/LDAS_v1.2.3/MIDW_025/observations_e5/sfx-trip/LAI_V2_2018-01-01_2018-12-31.PData')
swi_obs = pd.read_pickle('/cnrm/vegeo/albergelc/ldas_chain_python/LDAS_v1.2.3/MIDW_025/observations_e5/sfx-trip/raw_SWI_nc_2018-01-01_2018-12-31.PData')

rain0 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens0/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData')
lai0 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens0/results/ol/LAI_2018-01-01_2018-12-31.PData')
swi0 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens0/results/ol/WG2_2018-01-01_2018-12-31.PData')

### For Future - This will import all data directories into list - Would still need to read_pickle ###
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens*/results/ol/WG2_2018-01-01_2018-12-31.PData'),key=numericalSort)
'''
#rain = [pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens1/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens2/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens3/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens4/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens5/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens6/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens7/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens8/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens9/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens10/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens11/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens12/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens13/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens14/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens15/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens16/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens17/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens18/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens19/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens20/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens21/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens22/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens23/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens24/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens25/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens26/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens27/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens28/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens29/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens30/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens31/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens32/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens33/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens34/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens35/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens36/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens37/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens38/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens39/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens40/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens41/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens42/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens43/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens44/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens45/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens46/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens47/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens48/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens49/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens50/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData')]

#lai = [pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens1/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens2/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens3/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens4/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens5/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens6/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens7/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens8/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens9/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens10/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens11/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens12/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens13/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens14/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens15/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens16/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens17/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens18/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens19/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens20/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens21/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens22/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens23/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens24/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens25/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens26/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens27/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens28/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens29/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens30/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens31/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens32/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens33/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens34/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens35/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens36/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens37/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens38/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens39/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens40/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens41/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens42/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens43/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens44/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens45/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens46/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens47/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens48/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens49/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens50/results/ol/LAI_2018-01-01_2018-12-31.PData')]

#swi = [pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens1/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens2/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens3/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens4/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens5/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens6/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens7/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens8/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens9/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens10/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens11/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens12/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens13/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens14/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens15/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens16/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens17/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens18/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens19/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens20/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens21/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens22/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens23/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens24/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens25/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens26/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens27/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens28/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens29/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens30/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens31/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens32/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens33/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens34/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens35/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens36/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens37/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens38/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens39/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens40/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens41/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens42/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens43/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens44/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens45/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens46/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens47/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens48/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens49/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens50/results/ol/WG2_2018-01-01_2018-12-31.PData')]

#gpp = [pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens1/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens2/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens3/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens4/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens5/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens6/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens7/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens8/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens9/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens10/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens11/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens12/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens13/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens14/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens15/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens16/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens17/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens18/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens19/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens20/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens21/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens22/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens23/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens24/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens25/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens26/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens27/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens28/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens29/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens30/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens31/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens32/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens33/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens34/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens35/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens36/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens37/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens38/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens39/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens40/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens41/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens42/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens43/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens44/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens45/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens46/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens47/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens48/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens49/results/ol/GPPC_2018-01-01_2018-12-31.PData'),\
        #pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens50/results/ol/GPPC_2018-01-01_2018-12-31.PData')]

#eva = [pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens1/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens2/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens3/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens4/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens5/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens6/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens7/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens8/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens9/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens10/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens11/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens12/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens13/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens14/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens15/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens16/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens17/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens18/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens19/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens20/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens21/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens22/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens23/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens24/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens25/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens26/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens27/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens28/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens29/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens30/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens31/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens32/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens33/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens34/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens35/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens36/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens37/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens38/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens39/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens40/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens41/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens42/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens43/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens44/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens45/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens46/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens47/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens48/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens49/results/ol/EVAPC_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens50/results/ol/EVAPC_2018-01-01_2018-12-31.PData')]
'''
lai_d = [pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens5/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
       pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens10/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens15/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens20/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens25/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens30/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens35/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens40/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens45/results/ol/LAI_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens50/results/ol/LAI_2018-01-01_2018-12-31.PData')]

swi_d = [pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens5/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens10/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens15/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens20/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens25/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens30/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens35/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens40/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens45/results/ol/WG2_2018-01-01_2018-12-31.PData'),\
       pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens50/results/ol/WG2_2018-01-01_2018-12-31.PData')]

rain_d = [pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens5/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
       pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens10/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens15/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens20/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens25/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens30/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens35/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens40/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens45/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData'),\
        pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens50/results/ol/RAINFC_ISBA_2018-01-01_2018-12-31.PData')]

r = {}
s = {}
l = {}
g = {}
e = {}
_4 = {}
_6 = {}
_8 = {}

for iii in  range(1,51):
    rname = "/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens"+str(iii)+"/results/ol/RAINFC*"
    sname = "/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens"+str(iii)+"/results/ol/WG2*"
    lname = "/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens"+str(iii)+"/results/ol/LAI*"
    gname = "/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens"+str(iii)+"/results/ol/GPPC*"
    ename = "/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens"+str(iii)+"/results/ol/EVAPC*"
    name4 = "/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens"+str(iii)+"/results/ol/WG4*"
    name6 = "/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens"+str(iii)+"/results/ol/WG6*"
    name8 = "/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens"+str(iii)+"/results/ol/WG8*"
    
    
    for file in glob.glob(rname):
        name = file
        r[name] = pd.DataFrame(pd.read_pickle(file))
    for file in glob.glob(sname):
        name = file
        s[name] = pd.DataFrame(pd.read_pickle(file))
    for file in glob.glob(lname):
        name = file
        l[name] = pd.DataFrame(pd.read_pickle(file))
    for file in glob.glob(gname):
        name = file
        g[name] = pd.DataFrame(pd.read_pickle(file))
    for file in glob.glob(ename):
        name = file
        e[name] = pd.DataFrame(pd.read_pickle(file))
    for file in glob.glob(name4):
        name = file
        _4[name] = pd.DataFrame(pd.read_pickle(file))
    for file in glob.glob(name6):
        name = file
        _6[name] = pd.DataFrame(pd.read_pickle(file))
    for file in glob.glob(name8):
        name = file
        _8[name] = pd.DataFrame(pd.read_pickle(file))

rain = pd.Panel(r)
swi = pd.Panel(s)
lai = pd.Panel(l)
gpp = pd.Panel(g)
eva = pd.Panel(e)
wg4 = pd.Panel(_4)
wg6 = pd.Panel(_6)
wg8 = pd.Panel(_8)

era_5lai = pd.read_pickle('Data/LAI_2018-01-01_2018-12-31_mw.PData')
era_5swi = pd.read_pickle('Data/WG2_2018-01-01_2018-12-31_mw.PData')


print "Data Import Complete"

###################################################################
### Data Processing
###################################################################

### Calculates mean, standard deviation, and z (std/mean) using Panels
rstd = rain.std(axis=0)
rmean = rain.mean(axis=0)
rz = rstd/rmean

sstd = swi.std(axis=0)
smean = swi.mean(axis=0)
sz = sstd/smean

lstd = lai.std(axis=0)
lmean = lai.mean(axis=0)
lz = lstd/lmean

gstd = gpp.std(axis=0)
gmean = gpp.mean(axis=0)
gz = gstd/gmean

estd = eva.std(axis=0)
emean = eva.mean(axis=0)
ez = estd/emean

_4std = wg4.std(axis=0)
mean4 = wg4.mean(axis=0)
_4z = _4std/mean4

_6std = wg6.std(axis=0)
mean6 = wg6.mean(axis=0)
_6z = _6std/mean6

_8std = wg8.std(axis=0)
mean8 = wg8.mean(axis=0)
_8z = _8std/mean6

### Calculates annual and seasonal data and reshapes it for mapping
rz_win = rz.iloc[334:364].append(rz.iloc[0:59]).mean(axis=0)
rz_spr = rz.iloc[59:151].mean(axis=0)
rz_sum = rz.iloc[151:243].mean(axis=0)
rz_aut = rz.iloc[243:334].mean(axis=0)

sz_win = sz.iloc[334:364].append(sz.iloc[0:59]).mean(axis=0)
sz_spr = sz.iloc[59:151].mean(axis=0)
sz_sum = sz.iloc[151:243].mean(axis=0)
sz_aut = sz.iloc[243:334].mean(axis=0)

lz_win = lz.iloc[334:364].append(lz.iloc[0:59]).mean(axis=0)
lz_spr = lz.iloc[59:151].mean(axis=0)
lz_sum = lz.iloc[151:243].mean(axis=0)
lz_aut = lz.iloc[243:334].mean(axis=0)

gz_win = gz.iloc[334:364].append(gz.iloc[0:59]).mean(axis=0)
gz_spr = gz.iloc[59:151].mean(axis=0)
gz_sum = gz.iloc[151:243].mean(axis=0)
gz_aut = gz.iloc[243:334].mean(axis=0)

ez_win = ez.iloc[334:364].append(ez.iloc[0:59]).mean(axis=0)
ez_spr = ez.iloc[59:151].mean(axis=0)
ez_sum = ez.iloc[151:243].mean(axis=0)
ez_aut = ez.iloc[243:334].mean(axis=0)

_4z_win = _4z.iloc[334:364].append(_4z.iloc[0:59]).mean(axis=0)
_4z_spr = _4z.iloc[59:151].mean(axis=0)
_4z_sum = _4z.iloc[151:243].mean(axis=0)
_4z_aut = _4z.iloc[243:334].mean(axis=0)

_6z_win = _6z.iloc[334:364].append(_6z.iloc[0:59]).mean(axis=0)
_6z_spr = _6z.iloc[59:151].mean(axis=0)
_6z_sum = _6z.iloc[151:243].mean(axis=0)
_6z_aut = _6z.iloc[243:334].mean(axis=0)

_8z_win = _8z.iloc[334:364].append(_8z.iloc[0:59]).mean(axis=0)
_8z_spr = _8z.iloc[59:151].mean(axis=0)
_8z_sum = _8z.iloc[151:243].mean(axis=0)
_8z_aut = _8z.iloc[243:334].mean(axis=0)

ZR1 = rz_win.reshape(52,80);ZR2 = rz_spr.reshape(52,80);ZR3 = rz_sum.reshape(52,80);ZR4 = rz_aut.reshape(52,80)
ZS1 = sz_win.reshape(52,80);ZS2 = sz_spr.reshape(52,80);ZS3 = sz_sum.reshape(52,80);ZS4 = sz_aut.reshape(52,80)   
ZL1 = lz_win.reshape(52,80);ZL2 = lz_spr.reshape(52,80);ZL3 = lz_sum.reshape(52,80);ZL4 = lz_aut.reshape(52,80)
ZG1 = gz_win.reshape(52,80);ZG2 = gz_spr.reshape(52,80);ZG3 = gz_sum.reshape(52,80);ZG4 = gz_aut.reshape(52,80)
ZE1 = ez_win.reshape(52,80);ZE2 = ez_spr.reshape(52,80);ZE3 = ez_sum.reshape(52,80);ZE4 = ez_aut.reshape(52,80)

_4Z1 = _4z_win.reshape(52,80);_4Z2 = _4z_spr.reshape(52,80);_4Z3 = _4z_sum.reshape(52,80);_4Z4 = _4z_aut.reshape(52,80)
_6Z1 = _6z_win.reshape(52,80);_6Z2 = _6z_spr.reshape(52,80);_6Z3 = _6z_sum.reshape(52,80);_6Z4 = _6z_aut.reshape(52,80)
_8Z1 = _8z_win.reshape(52,80);_8Z2 = _8z_spr.reshape(52,80);_8Z3 = _8z_sum.reshape(52,80);_8Z4 = _8z_aut.reshape(52,80)

rz = rz.mean(axis=0)
lz = lz.mean(axis=0)
sz = sz.mean(axis=0)
gz = gz.mean(axis=0)
ez = ez.mean(axis=0)
_4z = _4z.mean(axis=0)
_6z = _6z.mean(axis=0)
_8z = _8z.mean(axis=0)

Z1 = rz.reshape(52,80)
Z2 = lz.reshape(52,80)
Z3 = sz.reshape(52,80)
Z4 = gz.reshape(52,80)
Z5 = ez.reshape(52,80)
Z6 = _4z.reshape(52,80)
Z7 = _6z.reshape(52,80)
Z8 = _8z.reshape(52,80)

print "Data Processing Complete"

'''
#raw_input("Press Enter to continue...")

var = np.empty((rain[0].shape[0],rain[0].shape[1]))
var[:] = np.nan
var = pd.DataFrame(var)

val = np.empty((lai[0].shape[0],lai[0].shape[1]))
val[:] = np.nan
val = pd.DataFrame(val)

vas = np.empty((swi[0].shape[0],swi[0].shape[1]))
vas[:] = np.nan
vas = pd.DataFrame(vas)

vag = np.empty((gpp[0].shape[0],gpp[0].shape[1]))
vag[:] = np.nan
vag = pd.DataFrame(vag)

vae = np.empty((eva[0].shape[0],eva[0].shape[1]))
vae[:] = np.nan
vae = pd.DataFrame(vae)

mar = np.empty((rain[0].shape[0],rain[0].shape[1]))
mar[:] = np.nan
mar = pd.DataFrame(mar)

mal = np.empty((lai[0].shape[0],lai[0].shape[1]))
mal[:] = np.nan
mal = pd.DataFrame(mal)

mas = np.empty((swi[0].shape[0],swi[0].shape[1]))
mas[:] = np.nan
mas = pd.DataFrame(mas)

mag = np.empty((gpp[0].shape[0],gpp[0].shape[1]))
mag[:] = np.nan
mag = pd.DataFrame(mag)

mae = np.empty((eva[0].shape[0],eva[0].shape[1]))
mae[:] = np.nan
mae = pd.DataFrame(mae)

zar = np.empty((rain[0].shape[0],rain[0].shape[1]))
zar[:] = np.nan
zar = pd.DataFrame(zar)

zal = np.empty((lai[0].shape[0],lai[0].shape[1]))
zal[:] = np.nan
zal = pd.DataFrame(zal)

zas = np.empty((swi[0].shape[0],swi[0].shape[1]))
zas[:] = np.nan
zas = pd.DataFrame(zas)

zag = np.empty((gpp[0].shape[0],gpp[0].shape[1]))
zag[:] = np.nan
zag = pd.DataFrame(zag)

zae = np.empty((eva[0].shape[0],eva[0].shape[1]))
zae[:] = np.nan
zae = pd.DataFrame(zae)


while True:
    response = raw_input('Calculate Variance/Standard Deviation? (Y)/(N):  ')
    if response != "Y" and response != "N":
        print 'Error: Input Not Expected. Please enter (Y) or (N):  '
    else:
        break

if response == 'Y':
    print 'Calculating Standard Deviation. This will take a while...'
    for c in tqdm(xrange(rain[0].shape[1])):
        for r in xrange(rain[0].shape[0]):
            #print rain[0][c][r]
            indvar = (rain[0][c][r],rain[1][c][r],rain[2][c][r],rain[3][c][r],rain[4][c][r],rain[5][c][r],rain[6][c][r],rain[7][c][r],rain[8][c][r],rain[9][c][r],rain[10][c][r],rain[11][c][r],rain[12][c][r],rain[13][c][r],rain[14][c][r],rain[15][c][r],rain[16][c][r],rain[17][c][r],rain[18][c][r],rain[19][c][r],rain[20][c][r],rain[21][c][r],rain[22][c][r],rain[23][c][r],rain[24][c][r],rain[25][c][r],rain[26][c][r],rain[27][c][r],rain[28][c][r],rain[29][c][r],rain[30][c][r],rain[31][c][r],rain[32][c][r],rain[33][c][r],rain[34][c][r],rain[35][c][r],rain[36][c][r],rain[37][c][r],rain[38][c][r],rain[39][c][r],rain[40][c][r],rain[41][c][r],rain[42][c][r],rain[43][c][r],rain[44][c][r],rain[45][c][r],rain[46][c][r],rain[47][c][r],rain[48][c][r],rain[49][c][r])
            
            #indval = (lai[0][c][r],lai[1][c][r],lai[2][c][r],lai[3][c][r],lai[4][c][r],lai[5][c][r],lai[6][c][r],lai[7][c][r],lai[8][c][r],lai[9][c][r],lai[10][c][r],lai[11][c][r],lai[12][c][r],lai[13][c][r],lai[14][c][r],lai[15][c][r],lai[16][c][r],lai[17][c][r],lai[18][c][r],lai[19][c][r],lai[20][c][r],lai[21][c][r],lai[22][c][r],lai[23][c][r],lai[24][c][r],lai[25][c][r],lai[26][c][r],lai[27][c][r],lai[28][c][r],lai[29][c][r],lai[30][c][r],lai[31][c][r],lai[32][c][r],lai[33][c][r],lai[34][c][r],lai[35][c][r],lai[36][c][r],lai[37][c][r],lai[38][c][r],lai[39][c][r],lai[40][c][r],lai[41][c][r],lai[42][c][r],lai[43][c][r],lai[44][c][r],lai[45][c][r],lai[46][c][r],lai[47][c][r],lai[48][c][r],lai[49][c][r])
            
            #indvas = (swi[0][c][r],swi[1][c][r],swi[2][c][r],swi[3][c][r],swi[4][c][r],swi[5][c][r],swi[6][c][r],swi[7][c][r],swi[8][c][r],swi[9][c][r],swi[10][c][r],swi[11][c][r],swi[12][c][r],swi[13][c][r],swi[14][c][r],swi[15][c][r],swi[16][c][r],swi[17][c][r],swi[18][c][r],swi[19][c][r],swi[20][c][r],swi[21][c][r],swi[22][c][r],swi[23][c][r],swi[24][c][r],swi[25][c][r],swi[26][c][r],swi[27][c][r],swi[28][c][r],swi[29][c][r],swi[30][c][r],swi[31][c][r],swi[32][c][r],swi[33][c][r],swi[34][c][r],swi[35][c][r],swi[36][c][r],swi[37][c][r],swi[38][c][r],swi[39][c][r],swi[40][c][r],swi[41][c][r],swi[42][c][r],swi[43][c][r],swi[44][c][r],swi[45][c][r],swi[46][c][r],swi[47][c][r],swi[48][c][r],swi[49][c][r])
            
            #indvag = (gpp[0][c][r],gpp[1][c][r],gpp[2][c][r],gpp[3][c][r],gpp[4][c][r],gpp[5][c][r],gpp[6][c][r],gpp[7][c][r],gpp[8][c][r],gpp[9][c][r],gpp[10][c][r],gpp[11][c][r],gpp[12][c][r],gpp[13][c][r],gpp[14][c][r],gpp[15][c][r],gpp[16][c][r],gpp[17][c][r],gpp[18][c][r],gpp[19][c][r],gpp[20][c][r],gpp[21][c][r],gpp[22][c][r],gpp[23][c][r],gpp[24][c][r],gpp[25][c][r],gpp[26][c][r],gpp[27][c][r],gpp[28][c][r],gpp[29][c][r],gpp[30][c][r],gpp[31][c][r],gpp[32][c][r],gpp[33][c][r],gpp[34][c][r],gpp[35][c][r],gpp[36][c][r],gpp[37][c][r],gpp[38][c][r],gpp[39][c][r],gpp[40][c][r],gpp[41][c][r],gpp[42][c][r],gpp[43][c][r],gpp[44][c][r],gpp[45][c][r],gpp[46][c][r],gpp[47][c][r],gpp[48][c][r],gpp[49][c][r])
            
            #indvae = (eva[0][c][r],eva[1][c][r],eva[2][c][r],eva[3][c][r],eva[4][c][r],eva[5][c][r],eva[6][c][r],eva[7][c][r],eva[8][c][r],eva[9][c][r],eva[10][c][r],eva[11][c][r],eva[12][c][r],eva[13][c][r],eva[14][c][r],eva[15][c][r],eva[16][c][r],eva[17][c][r],eva[18][c][r],eva[19][c][r],eva[20][c][r],eva[21][c][r],eva[22][c][r],eva[23][c][r],eva[24][c][r],eva[25][c][r],eva[26][c][r],eva[27][c][r],eva[28][c][r],eva[29][c][r],eva[30][c][r],eva[31][c][r],eva[32][c][r],eva[33][c][r],eva[34][c][r],eva[35][c][r],eva[36][c][r],eva[37][c][r],eva[38][c][r],eva[39][c][r],eva[40][c][r],eva[41][c][r],eva[42][c][r],eva[43][c][r],eva[44][c][r],eva[45][c][r],eva[46][c][r],eva[47][c][r],eva[48][c][r],eva[49][c][r])
        
            var[c][r] = np.nanstd(indvar,ddof=1)
            #val[c][r] = np.nanstd(indval)
            #vas[c][r] = np.nanstd(indvas)
            #vag[c][r] = np.nanstd(indvag)
            #vae[c][r] = np.nanstd(indvae)
            
            #mar[c][r] = np.nanmean(indvar)
            #mal[c][r] = np.nanmean(indval)
            #mas[c][r] = np.nanmean(indvas)
            #mag[c][r] = np.nanmean(indvag)
            #mae[c][r] = np.nanmean(indvae)
            
            #zar[c][r] = var[c][r]/mar[c][r]
            #zal[c][r] = val[c][r]/mal[c][r]
            #zas[c][r] = vas[c][r]/mas[c][r]
            #zag[c][r] = vag[c][r]/mag[c][r]
            #zae[c][r] = vae[c][r]/mae[c][r]
            
            
    print 'Variables Successfully Calculated'
    print 'Saving Files...'
    var.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/RainStd_dof1.PData')
    #val.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/LAIStd.PData')
    #vas.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/SWIStd.PData')
    #vag.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/GPPStd.PData')
    #vae.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/EVAStd.PData')
    #mar.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/RainMean.PData')
    #mal.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/LAIMean.PData')
    #mas.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/SWIMean.PData')
    #mag.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/GPPMean.PData')
    #mae.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/EVAMean.PData')
    #zar.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/RainZ.PData')
    #zal.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/LAIZ.PData')
    #zas.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/SWIZ.PData')
    #zag.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/GPPZ.PData')
    #zae.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/EVAZ.PData')
    print 'Variables Successfully Saved'
    
else:
    print 'Trying to Read in Saved files...'
    var = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/RainStd.PData')
    #val = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/LAIStd.PData')
    #vas = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/SWIStd.PData')
    #vag = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/GPPStd.PData')
    #vae = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/EVAStd.PData')
    #mar = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/RainMean.PData')
    #mal = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/LAIMean.PData')
    #mas = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/SWIMean.PData')
    #mag = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/GPPMean.PData')
    #mae = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/EVAMean.PData')
    #zar = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/RainZ.PData')
    #zal = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/LAIZ.PData')
    #zas = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/SWIZ.PData')
    #zag = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/GPPZ.PData')
    #zae = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/EVAZ.PData')
    print 'Variables Successfully Read'
        
#raw_input("Press Enter to continue...")
'''

###################################################################
### Compute Correlations
###################################################################


###################################################################
### Graphing
###################################################################

while True:
    response = raw_input('Would you like Rain,SWI,LAI,GPP,EVAP Time Series? (Y)/(N):  ')
    if response != "Y" and response != "N":
        print 'Error: Input Not Expected. Please enter (Y) or (N):  '
    elif response == 'Y':
        from matplotlib.dates import DateFormatter
        myFmt = DateFormatter("%b")
        fig, ax = plt.subplots()
        plt.title('ENS Rainfall',fontsize=24)
        for i in xrange(len(rain)):
            ax.plot(rain.iloc[i].mean(axis=1),label='ENS'+str(i+1))
        ax.plot(rain0.mean(axis=1),label='ENS CTRL',color='black',linewidth=2)
        #plt.legend(loc='upper right',fontsize=8)
        plt.yticks(fontsize=22)
        plt.xticks(fontsize=22)
        plt.ylabel('Rainfall [mm/day]',fontsize=24)
        ax.xaxis.set_major_formatter(myFmt) 
        ax.locator_params(axis='y', nbins=4)
        plt.show()
        
        fig, ax = plt.subplots()
        plt.title('ENS SSM OL',fontsize=24)
        myFmt = DateFormatter("%b")
        for i in xrange(len(swi)):
            ax.plot(swi.iloc[i].mean(axis=1),label='ENS'+str(i+1))
        ax.plot(swi0.mean(axis=1),label='ENS CTRL',color='black',linestyle='--')
        ax.plot(swi_era.mean(axis=1),label='ERA5',color='r',linewidth=2)
        #plt.plot(swi_obs.mean(axis=1),label='Obs',color='green',linewidth=2,linestyle=':')
        #plt.legend(loc='upper right',fontsize=8)
        plt.yticks(fontsize=22)
        plt.xticks(fontsize=22)
        plt.ylabel('4cm SSM [%]',fontsize=24)
        ax.xaxis.set_major_formatter(myFmt) 
        ax.locator_params(axis='y', nbins=4)
        plt.show()

        fig, ax = plt.subplots()
        plt.title('ENS LAI OL',fontsize=24)
        myFmt = DateFormatter("%b")
        for i in xrange(len(swi)):
            plt.plot(lai.iloc[i].mean(axis=1),label='ENS'+str(i+1))
        plt.plot(lai0.mean(axis=1),label='ENS CTRL',color='black',linestyle='--')
        plt.plot(lai_era.mean(axis=1),label='ERA5',color='r',linewidth=2)
        plt.plot(lai_obs.mean(axis=1),label='Obs',color='g',linewidth=3,linestyle=':')
        plt.plot(lai_obsV2.mean(axis=1),label='Obs',color='b',linewidth=3,linestyle=':')
        #plt.legend(loc='upper right',fontsize=8)
        plt.yticks(fontsize=22)
        plt.xticks(fontsize=22)
        plt.ylabel('Leaf Area Index [$m^2$/$m^2$]',fontsize=24)
        ax.xaxis.set_major_formatter(myFmt) 
        ax.locator_params(axis='y', nbins=4)
        plt.show()
        
        fig, ax = plt.subplots()
        plt.title('ENS GPP OL',fontsize=24)
        myFmt = DateFormatter("%b")
        for i in xrange(len(swi)):
            plt.plot(gpp.iloc[i].mean(axis=1),label='ENS'+str(i+1))
        #plt.legend(loc='upper right',fontsize=8)
        plt.yticks(fontsize=22)
        plt.xticks(fontsize=22)
        plt.ylabel('Gross Primary Production [mg C/$m^2$/day]',fontsize=24)
        ax.xaxis.set_major_formatter(myFmt)
        plt.ylim([0,0.05])
        ax.locator_params(axis='y', nbins=4)
        plt.show()
        
        fig, ax = plt.subplots()
        plt.title('ENS EVAP OL',fontsize=24)
        myFmt = DateFormatter("%b")
        for i in xrange(len(swi)):
            plt.plot(eva.iloc[i].mean(axis=1),label='ENS'+str(i+1))
        #plt.legend(loc='upper right',fontsize=8)
        plt.yticks(fontsize=22)
        plt.xticks(fontsize=22)
        plt.ylabel('Evapotranspiration [mm/day]',fontsize=24)
        ax.xaxis.set_major_formatter(myFmt)
        ax.locator_params(axis='y', nbins=4)
        plt.show()
        
        
        break
    else:
        break
    
while True:
    response = raw_input('Would you like WG2, WG4, WG6, WG8 Time Series? (Y)/(N) - CURRENTLY NOT IMPLEMENTED:  ')
    if response != "Y" and response != "N":
        print 'Error: Input Not Expected. Please enter (Y) or (N):  '
    elif response == 'Y':      
        
        break
    else:
        break

while True:
    response = raw_input('Would you like ENS0, ENS5, ENS10.... Time Series? (Y)/(N):  ')
    if response != "Y" and response != "N":
        print 'Error: Input Not Expected. Please enter (Y) or (N):  '
    elif response == 'Y':
        plt.title('ENS Rainfall OL',fontsize=16)
        for i in xrange(len(swi_d)):
            plt.plot(rain_d[i].mean(axis=1),label='ENS'+str((i+1)*5))
        #plt.plot(swi0.mean(axis=1),label='ENS CTRL',color='black',linestyle='--')
        #plt.plot(swi_era.mean(axis=1),label='ERA5',color='r',linewidth=2)
        #plt.plot(swi_obs.mean(axis=1),label='Obs',color='green',linewidth=2,linestyle=':')
        axes = plt.gca()
        axes.set_ylim([0,18])
        plt.legend(loc='upper right',fontsize=8)
        plt.show()        
        
        plt.title('ENS SWI OL',fontsize=16)
        for i in xrange(len(swi_d)):
            plt.plot(swi_d[i].mean(axis=1),label='ENS'+str((i+1)*5))
        plt.plot(swi0.mean(axis=1),label='ENS CTRL',color='black',linestyle='--')
        plt.plot(swi_era.mean(axis=1),label='ERA5',color='r',linewidth=2)
        plt.plot(swi_obs.mean(axis=1),label='Obs',color='green',linewidth=2,linestyle=':')
        plt.legend(loc='upper right',fontsize=8)
        plt.show()

        plt.title('ENS LAI OL',fontsize=16)
        for i in xrange(len(swi_d)):
            plt.plot(lai_d[i].mean(axis=1),label='ENS'+str((i+1)*5))
        plt.plot(lai0.mean(axis=1),label='ENS CTRL',color='black',linestyle='--')
        plt.plot(lai_era.mean(axis=1),label='ERA5',color='r',linewidth=2)
        plt.plot(lai_obs.mean(axis=1),label='Obs',color='g',linewidth=3,linestyle=':')
        plt.plot(lai_obsV2.mean(axis=1),label='Obs',color='b',linewidth=3,linestyle=':')
        plt.legend(loc='upper right',fontsize=8)
        plt.show()
        
        break
    else:
        break
    



###################################################################
### Mapping
###################################################################
while True:
    response = raw_input('Would you like Seasonal Standard Deviation Maps? (Y)/(N):  ')
    if response != "Y" and response != "N":
        print 'Error: Input Not Expected. Please enter (Y) or (N):  '
    elif response == 'Y':
        
        ### Rainfall
        #plt.figure(figsize=(8,6),dpi=300,edgecolor='w')
        fig, axes = plt.subplots(2, 2)
        axes[0,0].set_title("Winter")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZR1,interpolation='none',cmap='YlGn',vmin=0.20,vmax=1.75)
        
        axes[0,1].set_title("Spring")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZR2,interpolation='none',cmap='YlGn',vmin=0.20,vmax=1.75)
        
        axes[1,0].set_title("Summer")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZR3,interpolation='none',cmap='YlGn',vmin=0.20,vmax=1.75)
        
        axes[1,1].set_title("Autumn")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZR4,interpolation='none',cmap='YlGn',vmin=0.20,vmax=1.75)
        
        cbar_ax = fig.add_axes([0.9, 0.15, 0.025, 0.7])
        #fig.tight_layout()
        fig.colorbar(im,cax=cbar_ax)
        fig.suptitle('Rainfall (Standard Deviation / Mean)- 25km',fontsize=16)
        #plt.savefig('images/ENS_Winter_rainfall.png',format='png',dpi=300)
        plt.show()
        
        ### SSM
        fig, axes = plt.subplots(2, 2)
        axes[0,0].set_title("Winter")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZS1,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[0,1].set_title("Spring")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZS2,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[1,0].set_title("Summer")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZS3,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[1,1].set_title("Autumn")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZS4,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        cbar_ax = fig.add_axes([0.9, 0.15, 0.025, 0.7])
        fig.colorbar(im,cax=cbar_ax)
        fig.suptitle('SSM (Standard Deviation / Mean)- 25km',fontsize=16) 
        plt.show()
        
        ### LAI
        fig, axes = plt.subplots(2, 2)
        axes[0,0].set_title("Winter")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZL1,interpolation='none',cmap='YlGn',vmin=0.0,vmax=0.35)
        
        axes[0,1].set_title("Spring")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZL2,interpolation='none',cmap='YlGn',vmin=0.0,vmax=0.35)
        
        axes[1,0].set_title("Summer")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZL3,interpolation='none',cmap='YlGn',vmin=0.0,vmax=0.35)
        
        axes[1,1].set_title("Autumn")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZL4,interpolation='none',cmap='YlGn',vmin=0.0,vmax=0.35)
        
        cbar_ax = fig.add_axes([0.9, 0.15, 0.025, 0.7])
        fig.colorbar(im,cax=cbar_ax)
        fig.suptitle('LAI (Standard Deviation/Mean)- 25km',fontsize=16) 
        plt.show()

        break
    else:
        break
    
    
while True:
    response = raw_input('Would you like Seasonal Standard Deviation Maps for WG2, WG4, WG6, WG8? (Y)/(N):  ')
    if response != "Y" and response != "N":
        print 'Error: Input Not Expected. Please enter (Y) or (N):  '
    elif response == 'Y':
        ### WG2
        fig, axes = plt.subplots(2, 2)
        axes[0,0].set_title("Winter")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZS1,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[0,1].set_title("Spring")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZS2,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[1,0].set_title("Summer")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZS3,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[1,1].set_title("Autumn")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(ZS4,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        cbar_ax = fig.add_axes([0.9, 0.15, 0.025, 0.7])
        fig.colorbar(im,cax=cbar_ax)
        fig.suptitle('WG2 (Standard Deviation / Mean)- 25km',fontsize=16) 
        plt.show()
        
        ### WG4
        fig, axes = plt.subplots(2, 2)
        axes[0,0].set_title("Winter")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_4Z1,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[0,1].set_title("Spring")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_4Z2,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[1,0].set_title("Summer")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_4Z3,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[1,1].set_title("Autumn")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_4Z4,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        cbar_ax = fig.add_axes([0.9, 0.15, 0.025, 0.7])
        fig.colorbar(im,cax=cbar_ax)
        fig.suptitle('WG4 (Standard Deviation / Mean)- 25km',fontsize=16) 
        plt.show()
        
        ### WG6
        fig, axes = plt.subplots(2, 2)
        axes[0,0].set_title("Winter")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_6Z1,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[0,1].set_title("Spring")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_6Z2,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[1,0].set_title("Summer")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_6Z3,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[1,1].set_title("Autumn")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_6Z4,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        cbar_ax = fig.add_axes([0.9, 0.15, 0.025, 0.7])
        fig.colorbar(im,cax=cbar_ax)
        fig.suptitle('WG6 (Standard Deviation / Mean)- 25km',fontsize=16) 
        plt.show()
        
        ### WG8
        fig, axes = plt.subplots(2, 2)
        axes[0,0].set_title("Winter")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_8Z1,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[0,1].set_title("Spring")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[0,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_8Z2,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[1,0].set_title("Summer")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,0])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_8Z3,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        
        axes[1,1].set_title("Autumn")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i',ax=axes[1,1])
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(_8Z4,interpolation='none',cmap='YlGn',vmin=0.003,vmax=0.15)
        cbar_ax = fig.add_axes([0.9, 0.15, 0.025, 0.7])
        fig.colorbar(im,cax=cbar_ax)
        fig.suptitle('WG8 (Standard Deviation / Mean)- 25km',fontsize=16) 
        plt.show()
        
        break
    else:
        break
    
while True:
    response = raw_input('Would you like Annual Standard Deviation Maps? (Y)/(N):  ')
    if response != "Y" and response != "N":
        print 'Error: Input Not Expected. Please enter (Y) or (N):  '
    elif response == 'Y':
        #var_m = var.mean(axis=0)
        #val_m = val.mean(axis=0)
        #vas_m = vas.mean(axis=0)
        #vag_m = vag.mean(axis=0)
        #vae_m = vae.mean(axis=0)
        
        #mar_m = mar.mean(axis=0)
        #mal_m = mal.mean(axis=0)
        #mas_m = mas.mean(axis=0)
        #mag_m = mag.mean(axis=0)
        #mae_m = mae.mean(axis=0)
        
        #zar_m = zar.mean(axis=0)
        #zal_m = zal.mean(axis=0)
        #zas_m = zas.mean(axis=0)
        #zag_m = zag.mean(axis=0)
        #zae_m = zae.mean(axis=0)        

        #V1 = var_m.reshape(52,80)
        #V2 = val_m.reshape(52,80)
        #V3 = vas_m.reshape(52,80)
        #V4 = vag_m.reshape(52,80)
        #V5 = vae_m.reshape(52,80)
        
        #M1 = mar_m.reshape(52,80)
        #M2 = mal_m.reshape(52,80)
        #M3 = mas_m.reshape(52,80)
        #M4 = mag_m.reshape(52,80)
        #M5 = mae_m.reshape(52,80)
        
        #Z1 = zar_m.reshape(52,80)
        #Z2 = zal_m.reshape(52,80)
        #Z3 = zas_m.reshape(52,80)
        #Z4 = zag_m.reshape(52,80)
        #Z5 = zae_m.reshape(52,80)
        

        plt.title("Annual Rainfall (Standard Deviation / Mean)- 25km")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i')
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(Z1,interpolation='none',cmap='YlGn')
        plt.colorbar()
        #plt.colorbar(V1,location='bottom',pad='5%')
        #cbar.set_label('Standard Deviation')
        plt.show()

        plt.title("Annual LAI (Standard Deviation / Mean) - 25km") 
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i')
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(Z2,interpolation='none',cmap='YlGn')
        plt.colorbar()
        plt.show()

        plt.title("Annual SSM (Standard Deviation / Mean) - 25km")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i')
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(Z3,interpolation='none',cmap='YlGn')
        plt.colorbar()
        plt.show()
        
        plt.title("Annual GPP (Standard Deviation / Mean) - 25km")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i')
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(Z4,interpolation='none',cmap='YlGn')
        plt.colorbar()
        plt.show()
        
        plt.title("Annual Evapotranspiration (Standard Deviation / Mean) - 25km")
        map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=37, urcrnrlon=-85, urcrnrlat=50,resolution='i')
        map.drawstates(linewidth=1.8)
        map.drawcoastlines(1.8)
        map.drawmapboundary()
        map.drawcountries(1.8)
        map.drawrivers()
        im = map.imshow(Z5,interpolation='none',cmap='YlGn')
        plt.colorbar()
        plt.show()
        
        break
    else:
        break


###################################################################
### MatPlotLib Finishing Touches
###################################################################
### for 'add_axes' it is [left,bottom,witdth,height], plt.figure(#) is seperate image, or on the same plot
#fig = plt.figure(0)
#cbaxes = fig.add_axes([0.08,0.85,0.7,0.025])
#fig.colorbar(im,ax=axes.ravel().tolist())
#cbaxes = matplotlib.colorbar.make_axes(location='bottom')
#cbar = fig.colorbar(im,cax=cbaxes,orientation='horizontal',ticks=[-.75,0,.75])
#cbar.ax.set_xticklabels(['-0.75','0','0.75'])

#save('../vegeoFigures/Yield-RainF', ext='ps', close=True, verbose=True)
plt.show()

print 'Program Execution Complete'




