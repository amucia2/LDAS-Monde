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
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/observations/sfx-trip/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
obs = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ol/LAI_20*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_2/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol2 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_3/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol3 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_4/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol4 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_5/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol5 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_6/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol6 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_7/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol7 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_8/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol8 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_9/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol9 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_10/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol10 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_11/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol11 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_12/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol12 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_13/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol13 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_14/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ol14 = pd.concat(df)


file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekf/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_2/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf2 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_3/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf3 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_4/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf4 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_5/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf5 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_6/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf6 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_7/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf7 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_8/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf8 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_9/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf9 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_10/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf10 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_11/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf11 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_12/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf12 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_13/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf13 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_14/LAI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
ekf14 = pd.concat(df)
print "LAI DATA IMPORT: COMPLETE"

### Testing Raw SWI vs SWI
### Raw SWI makes some crazy bad correlations - Raw SWI is far, far too high
### RAW SWI
#file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/observations/sfx-trip/SWI*.PData'),key=numericalSort)
#df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
#sobs = pd.concat(df)

### CDF Matched SWI
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/observations/sfx-trip/tmp/SWI*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sobs = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ol/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sol = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_2/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sol2 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_4/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sol4 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_6/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sol6 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_8/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sol8 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_10/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sol10 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_12/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sol12 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/olfc_14/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sol14 = pd.concat(df)


file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekf/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sekf = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_2/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sekf2 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_4/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sekf4 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_6/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sekf6 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_8/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sekf8 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_10/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sekf10 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_12/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sekf12 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekffc_14/WG2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
sekf14 = pd.concat(df)
print "SWI DATA IMPORT: COMPLETE"



###################################################################
### Data Processing
####################################################################

tmp = ol2.copy()*np.nan
tmp.ix[obs.index] = obs

#### Use the following for graphing time-series
t_obs = obs.mean(axis=1)

t_ol = ol[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol2 = ol2[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol3 = ol3[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol4 = ol4[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol5 = ol5[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol6 = ol6[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol7 = ol7[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol8 = ol8[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol9 = ol9[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol10 = ol10[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol11 = ol11[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol12 = ol12[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol13 = ol13[~np.isnan(tmp)].mean(axis=1).dropna()
t_ol14 = ol14[~np.isnan(tmp)].mean(axis=1).dropna()

t_ekf = ekf[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf2 = ekf2[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf3 = ekf3[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf4 = ekf4[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf5 = ekf5[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf6 = ekf6[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf7 = ekf7[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf8 = ekf8[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf9 = ekf9[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf10 = ekf10[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf11 = ekf11[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf12 = ekf12[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf13 = ekf13[~np.isnan(tmp)].mean(axis=1).dropna()
t_ekf14 = ekf14[~np.isnan(tmp)].mean(axis=1).dropna()
print "LAI TIME SERIES PROCESSING: COMPLETE"

sobs.index = sobs.index.shift(-1,freq='D')
tmp = sol2.copy()*np.nan
tmp.ix[sobs.index] = sobs

### Use the following for graphing time-series
st_obs = sobs.mean(axis=1)

st_ol = sol[~np.isnan(tmp)].mean(axis=1).dropna()
st_ol2 = sol2[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ol3 = sol3[~np.isnan(tmp)].mean(axis=1).dropna()
st_ol4 = sol4[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ol5 = sol5[~np.isnan(tmp)].mean(axis=1).dropna()
st_ol6 = sol6[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ol7 = sol7[~np.isnan(tmp)].mean(axis=1).dropna()
st_ol8 = sol8[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ol9 = sol9[~np.isnan(tmp)].mean(axis=1).dropna()
st_ol10 = sol10[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ol11 = sol11[~np.isnan(tmp)].mean(axis=1).dropna()
st_ol12 = sol12[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ol13 = sol13[~np.isnan(tmp)].mean(axis=1).dropna()
st_ol14 = sol14[~np.isnan(tmp)].mean(axis=1).dropna()

st_ekf = sekf[~np.isnan(tmp)].mean(axis=1).dropna()
st_ekf2 = sekf2[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ekf3 = sekf3[~np.isnan(tmp)].mean(axis=1).dropna()
st_ekf4 = sekf4[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ekf5 = sekf5[~np.isnan(tmp)].mean(axis=1).dropna()
st_ekf6 = sekf6[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ekf7 = sekf7[~np.isnan(tmp)].mean(axis=1).dropna()
st_ekf8 = sekf8[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ekf9 = sekf9[~np.isnan(tmp)].mean(axis=1).dropna()
st_ekf10 = sekf10[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ekf11 = sekf11[~np.isnan(tmp)].mean(axis=1).dropna()
st_ekf12 = sekf12[~np.isnan(tmp)].mean(axis=1).dropna()
#st_ekf13 = sekf13[~np.isnan(tmp)].mean(axis=1).dropna()
st_ekf14 = sekf14[~np.isnan(tmp)].mean(axis=1).dropna()
print "SWI TIME SERIES PROCESSING: COMPLETE"



### Use the following for mapping yearly results 
m_obs = obs.mean(axis=0)

m_ol = ol[~np.isnan(tmp)].dropna()
m_ol2 = ol2[~np.isnan(tmp)].dropna()
m_ol3 = ol3[~np.isnan(tmp)].dropna()
m_ol4 = ol4[~np.isnan(tmp)].dropna()
m_ol5 = ol5[~np.isnan(tmp)].dropna()
m_ol6 = ol6[~np.isnan(tmp)].dropna()
m_ol7 = ol7[~np.isnan(tmp)].dropna()
m_ol8 = ol8[~np.isnan(tmp)].dropna()
m_ol9 = ol9[~np.isnan(tmp)].dropna()
m_ol10 = ol10[~np.isnan(tmp)].dropna()
m_ol11 = ol11[~np.isnan(tmp)].dropna()
m_ol12 = ol12[~np.isnan(tmp)].dropna()
m_ol13 = ol13[~np.isnan(tmp)].dropna()
m_ol14 = ol14[~np.isnan(tmp)].dropna()

m_ekf = ekf[~np.isnan(tmp)].dropna()
m_ekf2 = ekf2[~np.isnan(tmp)].dropna()
m_ekf3 = ekf3[~np.isnan(tmp)].dropna()
m_ekf4 = ekf4[~np.isnan(tmp)].dropna()
m_ekf5 = ekf5[~np.isnan(tmp)].dropna()
m_ekf6 = ekf6[~np.isnan(tmp)].dropna()
m_ekf7 = ekf7[~np.isnan(tmp)].dropna()
m_ekf8 = ekf8[~np.isnan(tmp)].dropna()
m_ekf9 = ekf9[~np.isnan(tmp)].dropna()
m_ekf10 = ekf10[~np.isnan(tmp)].dropna()
m_ekf11 = ekf11[~np.isnan(tmp)].dropna()
m_ekf12 = ekf12[~np.isnan(tmp)].dropna()
m_ekf13 = ekf13[~np.isnan(tmp)].dropna()
m_ekf14 = ekf14[~np.isnan(tmp)].dropna()
print "MAP SERIES PROCESSING: COMPLETE"

###################################################################
### Compute Correlations
###################################################################

### Convert all Series into DataFrames in order to do 'corrwith' correlations
### That omit missing values from either data source

d_obs = pd.DataFrame(t_obs)
d_ol = pd.DataFrame(t_ol);d_ol2 = pd.DataFrame(t_ol2);d_ol3 = pd.DataFrame(t_ol3)
d_ol4 = pd.DataFrame(t_ol4);d_ol5 = pd.DataFrame(t_ol5);d_ol6 = pd.DataFrame(t_ol6)
d_ol7 = pd.DataFrame(t_ol7);d_ol8 = pd.DataFrame(t_ol8);d_ol9 = pd.DataFrame(t_ol9)
d_ol10 = pd.DataFrame(t_ol10);d_ol11 = pd.DataFrame(t_ol11);d_ol12 = pd.DataFrame(t_ol12)
d_ol13 = pd.DataFrame(t_ol13);d_ol14 = pd.DataFrame(t_ol14)

d_ekf = pd.DataFrame(t_ekf)
d_ekf2 = pd.DataFrame(t_ekf2)
d_ekf3 = pd.DataFrame(t_ekf3)
d_ekf4 = pd.DataFrame(t_ekf4)
d_ekf5 = pd.DataFrame(t_ekf5)
d_ekf6 = pd.DataFrame(t_ekf6)
d_ekf7 = pd.DataFrame(t_ekf7)
d_ekf8 = pd.DataFrame(t_ekf8)
d_ekf9 = pd.DataFrame(t_ekf9)
d_ekf10 = pd.DataFrame(t_ekf10)
d_ekf11 = pd.DataFrame(t_ekf11)
d_ekf12 = pd.DataFrame(t_ekf12)
d_ekf13 = pd.DataFrame(t_ekf13)
d_ekf14 = pd.DataFrame(t_ekf14)

sd_obs = pd.DataFrame(st_obs)
sd_ol = pd.DataFrame(st_ol)
sd_ol2 = pd.DataFrame(st_ol2)
#sd_ol3 = pd.DataFrame(st_ol3)
sd_ol4 = pd.DataFrame(st_ol4)
#sd_ol5 = pd.DataFrame(st_ol5)
sd_ol6 = pd.DataFrame(st_ol6)
#sd_ol7 = pd.DataFrame(st_ol7)
sd_ol8 = pd.DataFrame(st_ol8)
#sd_ol9 = pd.DataFrame(st_ol9)
sd_ol10 = pd.DataFrame(st_ol10)
#sd_ol11 = pd.DataFrame(st_ol11)
sd_ol12 = pd.DataFrame(st_ol12)
#sd_ol13 = pd.DataFrame(st_ol13)
sd_ol14 = pd.DataFrame(st_ol14)

sd_ekf = pd.DataFrame(st_ekf)
sd_ekf2 = pd.DataFrame(st_ekf2)
#sd_ekf3 = pd.DataFrame(st_ekf3)
sd_ekf4 = pd.DataFrame(st_ekf4)
#sd_ekf5 = pd.DataFrame(st_ekf5)
sd_ekf6 = pd.DataFrame(st_ekf6)
#sd_ekf7 = pd.DataFrame(st_ekf7)
sd_ekf8 = pd.DataFrame(st_ekf8)
#sd_ekf9 = pd.DataFrame(st_ekf9)
sd_ekf10 = pd.DataFrame(st_ekf10)
#sd_ekf11 = pd.DataFrame(st_ekf11)
sd_ekf12 = pd.DataFrame(st_ekf12)
#sd_ekf13 = pd.DataFrame(st_ekf13)
sd_ekf14 = pd.DataFrame(st_ekf14)

cor_ol_obs = d_obs.corrwith(d_ol,axis=0,drop=True)
cor_ol2_obs = d_obs.corrwith(d_ol2,axis=0,drop=True)
cor_ol3_obs = d_obs.corrwith(d_ol3,axis=0,drop=True)
cor_ol4_obs = d_obs.corrwith(d_ol4,axis=0,drop=True)
cor_ol5_obs = d_obs.corrwith(d_ol5,axis=0,drop=True)
cor_ol6_obs = d_obs.corrwith(d_ol6,axis=0,drop=True)
cor_ol7_obs = d_obs.corrwith(d_ol7,axis=0,drop=True)
cor_ol8_obs = d_obs.corrwith(d_ol8,axis=0,drop=True)
cor_ol9_obs = d_obs.corrwith(d_ol9,axis=0,drop=True)
cor_ol10_obs = d_obs.corrwith(d_ol10,axis=0,drop=True)
cor_ol11_obs = d_obs.corrwith(d_ol11,axis=0,drop=True)
cor_ol12_obs = d_obs.corrwith(d_ol12,axis=0,drop=True)
cor_ol13_obs = d_obs.corrwith(d_ol13,axis=0,drop=True)
cor_ol14_obs = d_obs.corrwith(d_ol14,axis=0,drop=True)

cor_ekf_obs = d_obs.corrwith(d_ekf,axis=0,drop=True)
cor_ekf2_obs = d_obs.corrwith(d_ekf2,axis=0,drop=True)
cor_ekf3_obs = d_obs.corrwith(d_ekf3,axis=0,drop=True)
cor_ekf4_obs = d_obs.corrwith(d_ekf4,axis=0,drop=True)
cor_ekf5_obs = d_obs.corrwith(d_ekf5,axis=0,drop=True)
cor_ekf6_obs = d_obs.corrwith(d_ekf6,axis=0,drop=True)
cor_ekf7_obs = d_obs.corrwith(d_ekf7,axis=0,drop=True)
cor_ekf8_obs = d_obs.corrwith(d_ekf8,axis=0,drop=True)
cor_ekf9_obs = d_obs.corrwith(d_ekf9,axis=0,drop=True)
cor_ekf10_obs = d_obs.corrwith(d_ekf10,axis=0,drop=True)
cor_ekf11_obs = d_obs.corrwith(d_ekf11,axis=0,drop=True)
cor_ekf12_obs = d_obs.corrwith(d_ekf12,axis=0,drop=True)
cor_ekf13_obs = d_obs.corrwith(d_ekf13,axis=0,drop=True)
cor_ekf14_obs = d_obs.corrwith(d_ekf14,axis=0,drop=True)
##
scor_ol_obs = sd_obs.corrwith(sd_ol,axis=0,drop=True)
scor_ol2_obs = sd_obs.corrwith(sd_ol2,axis=0,drop=True)
#scor_ol3_obs = sd_obs.corrwith(sd_ol3,axis=0,drop=True)
scor_ol4_obs = sd_obs.corrwith(sd_ol4,axis=0,drop=True)
#scor_ol5_obs = sd_obs.corrwith(sd_ol5,axis=0,drop=True)
scor_ol6_obs = sd_obs.corrwith(sd_ol6,axis=0,drop=True)
#scor_ol7_obs = sd_obs.corrwith(sd_ol7,axis=0,drop=True)
scor_ol8_obs = sd_obs.corrwith(sd_ol8,axis=0,drop=True)
#scor_ol9_obs = sd_obs.corrwith(sd_ol9,axis=0,drop=True)
scor_ol10_obs = sd_obs.corrwith(sd_ol10,axis=0,drop=True)
#scor_ol11_obs = sd_obs.corrwith(sd_ol11,axis=0,drop=True)
scor_ol12_obs = sd_obs.corrwith(sd_ol12,axis=0,drop=True)
#scor_ol13_obs = sd_obs.corrwith(sd_ol13,axis=0,drop=True)
scor_ol14_obs = sd_obs.corrwith(sd_ol14,axis=0,drop=True)

scor_ekf_obs = sd_obs.corrwith(sd_ekf,axis=0,drop=True)
scor_ekf2_obs = sd_obs.corrwith(sd_ekf2,axis=0,drop=True)
#scor_ekf3_obs = sd_obs.corrwith(sd_ekf3,axis=0,drop=True)
scor_ekf4_obs = sd_obs.corrwith(sd_ekf4,axis=0,drop=True)
#scor_ekf5_obs = sd_obs.corrwith(sd_ekf5,axis=0,drop=True)
scor_ekf6_obs = sd_obs.corrwith(sd_ekf6,axis=0,drop=True)
#scor_ekf7_obs = sd_obs.corrwith(sd_ekf7,axis=0,drop=True)
scor_ekf8_obs = sd_obs.corrwith(sd_ekf8,axis=0,drop=True)
#scor_ekf9_obs = sd_obs.corrwith(sd_ekf9,axis=0,drop=True)
scor_ekf10_obs = sd_obs.corrwith(sd_ekf10,axis=0,drop=True)
#scor_ekf11_obs = sd_obs.corrwith(sd_ekf11,axis=0,drop=True)
scor_ekf12_obs = sd_obs.corrwith(sd_ekf12,axis=0,drop=True)
#scor_ekf13_obs = sd_obs.corrwith(sd_ekf13,axis=0,drop=True)
scor_ekf14_obs = sd_obs.corrwith(sd_ekf14,axis=0,drop=True)

### LAI and SWI 
### [LAI_Correlation, LAI_RMSD, SWI_Correlation, SWI_RMSD]
### EVAP
cor_ol_obs = [0.686,1.187,0.864,0.045]
cor_ol2_obs = [0.683,1.193,0.856,0.046]
cor_ol3_obs = [0.685,1.196,0.849,0.047]
cor_ol4_obs = [0.685,1.196,0.841,0.048]
cor_ol5_obs = [0.689,1.188,0.828,0.050]
cor_ol6_obs = [0.682,1.202,0.815,0.052]
cor_ol7_obs = [0.690,1.169,0.805,0.053]
cor_ol8_obs = [0.685,1.210,0.792,0.055]
cor_ol9_obs = [0.702,1.180,0.777,0.057]
cor_ol10_obs = [0.689,1.207,0.766,0.058]
cor_ol11_obs = [0.686,1.215,0.757,0.059]
cor_ol12_obs = [0.683,1.219,0.746,0.061]
cor_ol13_obs = [0.685,1.210,0.740,0.061]
cor_ol14_obs = [0.689,1.213,0.736,0.062]

cor_ekf_obs = [0.870,0.716,0.879,0.042]
cor_ekf2_obs = [0.811,0.858,0.864,0.044]
cor_ekf3_obs = [0.811,0.864,0.855,0.046]
cor_ekf4_obs = [0.814,0.864,0.845,0.047]
cor_ekf5_obs = [0.814,0.861,0.831,0.049]
cor_ekf6_obs = [0.811,0.860,0.818,0.051]
cor_ekf7_obs = [0.813,0.870,0.807,0.053]
cor_ekf8_obs = [0.815,0.866,0.793,0.054]
cor_ekf9_obs = [0.822,0.847,0.778,0.056]
cor_ekf10_obs = [0.816,0.857,0.768,0.058]
cor_ekf11_obs = [0.781,0.945,0.758,0.059]
cor_ekf12_obs = [0.772,0.961,0.747,0.060]
cor_ekf13_obs = [0.771,0.965,0.741,0.061]
cor_ekf14_obs = [0.776,0.965,0.736,0.061]

ecor_ol_obs = [0.686,1.187,0.864,0.045]
ecor_ol2_obs = [0.683,1.193,0.856,0.046]
#ecor_ol3_obs = [0.685,1.196,0.849,0.047]
ecor_ol4_obs = [0.685,1.196,0.841,0.048]
#ecor_ol5_obs = [0.689,1.188,0.828,0.050]
ecor_ol6_obs = [0.682,1.202,0.815,0.052]
#ecor_ol7_obs = [0.690,1.169,0.805,0.053]
ecor_ol8_obs = [0.685,1.210,0.792,0.055]
#ecor_ol9_obs = [0.702,1.180,0.777,0.057]
ecor_ol10_obs = [0.689,1.207,0.766,0.058]
#ecor_ol11_obs = [0.686,1.215,0.757,0.059]
ecor_ol12_obs = [0.683,1.219,0.746,0.061]
#ecor_ol13_obs = [0.685,1.210,0.740,0.061]
ecor_ol14_obs = [0.689,1.213,0.736,0.062]

ecor_ekf_obs = [0.870,0.716,0.879,0.042]
ecor_ekf2_obs = [0.811,0.858,0.864,0.044]
#ecor_ekf3_obs = [0.811,0.864,0.855,0.046]
ecor_ekf4_obs = [0.814,0.864,0.845,0.047]
#ecor_ekf5_obs = [0.814,0.861,0.831,0.049]
ecor_ekf6_obs = [0.811,0.860,0.818,0.051]
#ecor_ekf7_obs = [0.813,0.870,0.807,0.053]
ecor_ekf8_obs = [0.815,0.866,0.793,0.054]
#ecor_ekf9_obs = [0.822,0.847,0.778,0.056]
ecor_ekf10_obs = [0.816,0.857,0.768,0.058]
#ecor_ekf11_obs = [0.781,0.945,0.758,0.059]
ecor_ekf12_obs = [0.772,0.961,0.747,0.060]
#ecor_ekf13_obs = [0.771,0.965,0.741,0.061]
ecor_ekf14_obs = [0.776,0.965,0.736,0.061]




print("LAI Correlations")
print("------------------------------------------------------------------")
print("Correlation - OL vs Obs: " + str(round(cor_ol_obs[0],4)))
print("Correlation - OL 2 Day FC vs Obs: " + str(round(cor_ol2_obs[0],4)))
print("Correlation - OL 3 Day FC vs Obs: " + str(round(cor_ol3_obs[0],4)))
print("Correlation - OL 4 Day FC vs Obs: " + str(round(cor_ol4_obs[0],4)))
print("Correlation - OL 5 Day FC vs Obs: " + str(round(cor_ol5_obs[0],4)))
print("Correlation - OL 6 Day FC vs Obs: " + str(round(cor_ol6_obs[0],4)))
print("Correlation - OL 7 Day FC vs Obs: " + str(round(cor_ol7_obs[0],4)))
print("Correlation - OL 8 Day FC vs Obs: " + str(round(cor_ol8_obs[0],4)))
print("Correlation - OL 9 Day FC vs Obs: " + str(round(cor_ol9_obs[0],4)))
print("Correlation - OL 10 Day FC vs Obs: " + str(round(cor_ol10_obs[0],4)))
print("Correlation - OL 11 Day FC vs Obs: " + str(round(cor_ol11_obs[0],4)))
print("Correlation - OL 12 Day FC vs Obs: " + str(round(cor_ol12_obs[0],4)))
print("Correlation - OL 13 Day FC vs Obs: " + str(round(cor_ol13_obs[0],4)))
print("Correlation - OL 14 Day FC vs Obs: " + str(round(cor_ol14_obs[0],4)))
print("------------------------------------------------------------------")
print("Correlation - EKF vs Obs: " + str(round(cor_ekf_obs[0],4)))
print("Correlation - EKF 2 Day FC vs Obs: " + str(round(cor_ekf2_obs[0],4)))
print("Correlation - EKF 3 Day FC vs Obs: " + str(round(cor_ekf3_obs[0],4)))
print("Correlation - EKF 4 Day FC vs Obs: " + str(round(cor_ekf4_obs[0],4)))
print("Correlation - EKF 5 Day FC vs Obs: " + str(round(cor_ekf5_obs[0],4)))
print("Correlation - EKF 6 Day FC vs Obs: " + str(round(cor_ekf6_obs[0],4)))
print("Correlation - EKF 7 Day FC vs Obs: " + str(round(cor_ekf7_obs[0],4)))
print("Correlation - EKF 8 Day FC vs Obs: " + str(round(cor_ekf8_obs[0],4)))
print("Correlation - EKF 9 Day FC vs Obs: " + str(round(cor_ekf9_obs[0],4)))
print("Correlation - EKF 10 Day FC vs Obs: " + str(round(cor_ekf10_obs[0],4)))
print("Correlation - EKF 11 Day FC vs Obs: " + str(round(cor_ekf11_obs[0],4)))
print("Correlation - EKF 12 Day FC vs Obs: " + str(round(cor_ekf12_obs[0],4)))
print("Correlation - EKF 13 Day FC vs Obs: " + str(round(cor_ekf13_obs[0],4)))
print("Correlation - EKF 14 Day FC vs Obs: " + str(round(cor_ekf14_obs[0],4)))
print("------------------------------------------------------------------")

print("SWI Correlations")
print("------------------------------------------------------------------")
print("Correlation - OL vs Obs: " + str(round(scor_ol_obs[0],4)))
print("Correlation - OL 2 Day FC vs Obs: " + str(round(scor_ol2_obs[0],4)))
#print("Correlation - OL 3 Day FC vs Obs: " + str(round(scor_ol3_obs[0],4)))
print("Correlation - OL 4 Day FC vs Obs: " + str(round(scor_ol4_obs[0],4)))
#print("Correlation - OL 5 Day FC vs Obs: " + str(round(scor_ol5_obs[0],4)))
print("Correlation - OL 6 Day FC vs Obs: " + str(round(scor_ol6_obs[0],4)))
#print("Correlation - OL 7 Day FC vs Obs: " + str(round(scor_ol7_obs[0],4)))
print("Correlation - OL 8 Day FC vs Obs: " + str(round(scor_ol8_obs[0],4)))
#print("Correlation - OL 9 Day FC vs Obs: " + str(round(scor_ol9_obs[0],4)))
print("Correlation - OL 10 Day FC vs Obs: " + str(round(scor_ol10_obs[0],4)))
#print("Correlation - OL 11 Day FC vs Obs: " + str(round(scor_ol11_obs[0],4)))
print("Correlation - OL 12 Day FC vs Obs: " + str(round(scor_ol12_obs[0],4)))
#print("Correlation - OL 13 Day FC vs Obs: " + str(round(scor_ol13_obs[0],4)))
print("Correlation - OL 14 Day FC vs Obs: " + str(round(scor_ol14_obs[0],4)))
print("------------------------------------------------------------------")
print("Correlation - EKF vs Obs: " + str(round(scor_ekf_obs[0],4)))
print("Correlation - EKF 2 Day FC vs Obs: " + str(round(scor_ekf2_obs[0],4)))
#print("Correlation - EKF 3 Day FC vs Obs: " + str(round(scor_ekf3_obs[0],4)))
print("Correlation - EKF 4 Day FC vs Obs: " + str(round(scor_ekf4_obs[0],4)))
#print("Correlation - EKF 5 Day FC vs Obs: " + str(round(scor_ekf5_obs[0],4)))
print("Correlation - EKF 6 Day FC vs Obs: " + str(round(scor_ekf6_obs[0],4)))
#print("Correlation - EKF 7 Day FC vs Obs: " + str(round(scor_ekf7_obs[0],4)))
print("Correlation - EKF 8 Day FC vs Obs: " + str(round(scor_ekf8_obs[0],4)))
#print("Correlation - EKF 9 Day FC vs Obs: " + str(round(scor_ekf9_obs[0],4)))
print("Correlation - EKF 10 Day FC vs Obs: " + str(round(scor_ekf10_obs[0],4)))
#print("Correlation - EKF 11 Day FC vs Obs: " + str(round(scor_ekf11_obs[0],4)))
print("Correlation - EKF 12 Day FC vs Obs: " + str(round(scor_ekf12_obs[0],4)))
#print("Correlation - EKF 13 Day FC vs Obs: " + str(round(scor_ekf13_obs[0],4)))
print("Correlation - EKF 14 Day FC vs Obs: " + str(round(scor_ekf14_obs[0],4)))
print("------------------------------------------------------------------")

#lai_ol=[cor_ol_obs[0], cor_ol2_obs[0], cor_ol4_obs[0],cor_ol6_obs[0],cor_ol8_obs[0],cor_ol10_obs[0],cor_ol12_obs[0],cor_ol14_obs[0]]
#lai_ekf=[cor_ekf_obs[0], cor_ekf2_obs[0], cor_ekf4_obs[0],cor_ekf6_obs[0],cor_ekf8_obs[0],cor_ekf10_obs[0],cor_ekf12_obs[0],cor_ekf14_obs[0]]

#swi_ol=[scor_ol_obs[0], scor_ol2_obs[0], scor_ol4_obs[0],scor_ol6_obs[0],scor_ol8_obs[0],scor_ol10_obs[0],scor_ol12_obs[0],scor_ol14_obs[0]]
#swi_ekf=[scor_ekf_obs[0], scor_ekf2_obs[0], scor_ekf4_obs[0],scor_ekf6_obs[0],scor_ekf8_obs[0],scor_ekf10_obs[0],scor_ekf12_obs[0],scor_ekf14_obs[0]]

lai_ol=[cor_ol_obs[0], cor_ol2_obs[0], cor_ol4_obs[0],cor_ol6_obs[0],cor_ol8_obs[0],cor_ol10_obs[0],cor_ol12_obs[0],cor_ol14_obs[0]]
lai_ekf=[cor_ekf_obs[0], cor_ekf2_obs[0], cor_ekf4_obs[0],cor_ekf6_obs[0],cor_ekf8_obs[0],cor_ekf10_obs[0],cor_ekf12_obs[0],cor_ekf14_obs[0]]

swi_ol=[cor_ol_obs[2], cor_ol2_obs[2], cor_ol4_obs[2],cor_ol6_obs[2],cor_ol8_obs[2],cor_ol10_obs[2],cor_ol12_obs[2],cor_ol14_obs[2]]
swi_ekf=[cor_ekf_obs[2], cor_ekf2_obs[2], cor_ekf4_obs[2],cor_ekf6_obs[2],cor_ekf8_obs[2],cor_ekf10_obs[2],cor_ekf12_obs[2],cor_ekf14_obs[2]]

rlai_ol=[cor_ol_obs[1], cor_ol2_obs[1], cor_ol4_obs[1],cor_ol6_obs[1],cor_ol8_obs[1],cor_ol10_obs[1],cor_ol12_obs[1],cor_ol14_obs[1]]
rlai_ekf=[cor_ekf_obs[1], cor_ekf2_obs[1], cor_ekf4_obs[1],cor_ekf6_obs[1],cor_ekf8_obs[1],cor_ekf10_obs[1],cor_ekf12_obs[1],cor_ekf14_obs[1]]

rswi_ol=[cor_ol_obs[3], cor_ol2_obs[3], cor_ol4_obs[3],cor_ol6_obs[3],cor_ol8_obs[3],cor_ol10_obs[3],cor_ol12_obs[3],cor_ol14_obs[3]]
rswi_ekf=[cor_ekf_obs[3], cor_ekf2_obs[3], cor_ekf4_obs[3],cor_ekf6_obs[3],cor_ekf8_obs[3],cor_ekf10_obs[3],cor_ekf12_obs[3],cor_ekf14_obs[3]]

#    adfdfdkkkk1212121212asdf
#    asdfhjknmm3434343434asdf

###################################################################
### Graphing
###################################################################
'''
plt.title('LAI over CONUS',fontsize=24)
plt.plot(t_obs,label='Observations',marker='*',color='green',linewidth=0,markersize=10)

plt.plot(t_ol,label='OL',color='blue',linewidth=1)
#plt.plot(t_ol2,label='OL 2 Day FC',linestyle='--',linewidth=2,color='cyan')
#plt.plot(t_ol4,label='OL 4 Day FC',linestyle='--')
#plt.plot(t_ol6,label='OL 6 Day FC',linestyle='--',linewidth=2,color='yellow')
#plt.plot(t_ol8,label='OL 8 Day FC',linestyle='--')
plt.plot(t_ol10,label='OL 10 Day FC',linestyle='--',linewidth=1,color='blue')
#plt.plot(t_ol12,label='OL 12 Day FC',linestyle='--')
#plt.plot(t_ol14,label='OL 14 Day FC',linestyle='--',linewidth=2,color='black')

plt.plot(t_ekf,label='EKF',color='red',linewidth=1)
#plt.plot(t_ekf2,label='EKF 2 Day FC',linestyle=':',linewidth=2,color='cyan')
#plt.plot(t_ekf4,label='EKF 4 Day FC',linestyle=':')
#plt.plot(t_ekf6,label='EKF 6 Day FC',linestyle=':',linewidth=2,color='yellow')
#plt.plot(t_ekf8,label='EKF 8 Day FC',linestyle=':')
plt.plot(t_ekf10,label='EKF 10 Day FC',linestyle=':',linewidth=1,color='red')
#plt.plot(t_ekf12,label='EKF 12 Day FC',linestyle=':')
#plt.plot(t_ekf14,label='EKF 14 Day FC',linestyle=':',linewidth=2,color='black')

#plt.axhline(color='black',linewidth=.5)
plt.legend(loc='upper left')
#plt.ylim([0,1])
plt.ylabel('Leaf Area Index [$m^2$/$m^2$]',fontsize=24)
plt.yticks(fontsize=22)
plt.xticks(fontsize=22)
plt.locator_params(axis='y', nbins=4)


#plt.ylabel('',fontsize=24)

#plt.rcParams.update({'font.size': 10})
plt.show()


plt.title('SWI over CONUS',fontsize=24)
plt.plot(st_obs,label='Observations',marker='*',color='green',linewidth=0,markersize=10)

plt.plot(st_ol,label='OL',color='blue',linewidth=1)
#plt.plot(st_ol2,label='OL 2 Day FC',linestyle='--',linewidth=2,color='cyan')
#plt.plot(st_ol4,label='OL 4 Day FC',linestyle='--')
#plt.plot(st_ol6,label='OL 6 Day FC',linestyle='--',linewidth=2,color='yellow')
#plt.plot(st_ol8,label='OL 8 Day FC',linestyle='--')
plt.plot(st_ol10,label='OL 10 Day FC',linestyle='--',linewidth=1,color='blue')
#plt.plot(st_ol12,label='OL 12 Day FC',linestyle='--')
#plt.plot(st_ol14,label='OL 14 Day FC',linestyle='--',linewidth=2,color='black')

plt.plot(st_ekf,label='EKF',color='red',linewidth=1)
#plt.plot(st_ekf2,label='EKF 2 Day FC',linestyle=':',linewidth=2,color='cyan')
#plt.plot(st_ekf4,label='EKF 4 Day FC',linestyle=':')
#plt.plot(st_ekf6,label='EKF 6 Day FC',linestyle=':',linewidth=2,color='yellow')
#plt.plot(st_ekf8,label='EKF 8 Day FC',linestyle=':')
plt.plot(st_ekf10,label='EKF 10 Day FC',linestyle=':',linewidth=1,color='red')
#plt.plot(st_ekf12,label='EKF 12 Day FC',linestyle=':')
#plt.plot(st_ekf14,label='EKF 14 Day FC',linestyle=':',linewidth=2,color='black')

#plt.axhline(color='black',linewidth=.5)
plt.legend(loc='upper left')
plt.ylim([0.2,0.35])
#plt.ylabel('Leaf Area Index [$m^2$/$m^2$]',fontsize=24)
plt.ylabel('SWI',fontsize=24)
#plt.rcParams.update({'font.size': 10})
plt.yticks(fontsize=22)
plt.xticks(fontsize=22)
plt.locator_params(axis='y', nbins=4)
plt.show()
'''

#ollai = ['OL','OL FC2','OL FC4', 'OL FC6', 'OL FC8', 'OL FC10', 'OL FC12', 'OL FC14']
#ekflai = ['EKF','EKF FC2','EKF FC4', 'EKF FC6', 'EKF FC8', 'EKF FC10', 'EKF FC12', 'EKF FC14']
ticks = ['NO FC','FC2','FC4', 'FC6', 'FC8', 'FC10', 'FC12', 'FC14']

fig, ax1 = plt.subplots()

plt.title('LAI Forecast Satellite Correlations over CONUS',fontsize=24)
ax1.set_ylabel('R',fontsize=24)
ax1.set_xlabel('Forecast Day',fontsize=24)
plt.xticks(np.arange(8),ticks,fontsize=20)
plt.yticks(fontsize=20)
#ax1.yaxis.set_major_locator(plt.MaxNLocator(5))
#plt.locator_params(axis='y', nbins=4)
ax1.set_yticks([0.7,0.75,0.8,0.85,0.9])
ax1.plot(lai_ol,markersize=20,marker='.',linewidth=1,label='OL',color='b',linestyle='--')
ax1.plot(lai_ekf,markersize=20,marker='.',linewidth=1,label='EKF',color='r',linestyle='--')

#ax1.plot(swi_ol,markersize=20,marker='.',linewidth=1,label='OL',color='b',linestyle='--')
#ax1.plot(swi_ekf,markersize=20,marker='.',linewidth=1,label='EKF',color='r',linestyle ='--')

ax1.margins(0.05)
ax1.legend(loc='upper right',fontsize=24)
plt.show()

#ax2=ax1.twinx()
#ax2.set_ylabel('RMSD',fontsize=16)
#ax2.plot(rlai_ol,markersize=12,marker='*',linewidth=0.5,label='OL LAI RMSD',color='b',linestyle='-')
#ax2.plot(rlai_ekf,markersize=12,marker='*',linewidth=0.5,label='EKF LAI RMSD',color='r',linestyle='-')
##ax2.plot(rswi_ol,markersize=12,marker='*',linewidth=0.5,label='OL SWI RMSD',color='c',linestyle='--')
##ax2.plot(rswi_ekf,markersize=12,marker='*',linewidth=0.5,label='EKF SWI RMSD',color='m',linestyle='--')
#ax2.margins(0.05)
#plt.legend(loc='upper right')



#plt.title('SWI Correlations',fontsize=16)
#plt.plot(scor_ol_obs,label='OL',color='blue',linewidth=1)





###################################################################
### Mapping
###################################################################




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
#plt.show()


