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
'''
### Import individual yearly files
lai = np.load('Nebraska_HR/results/ol/LAI_2017-01-01_2017-12-31.PData')
lai_an = np.load('Nebraska_HR/results/ekf/LAI_2017-01-01_2017-12-31.PData')
lai_obs = np.load('Nebraska_HR/observations/sfx-trip/LAI_V2_2017-01-01_2017-12-31.PData')

lai1 = np.load('Nebraska_LDAS/results/ol/LAI_2012-01-01_2012-12-31.PData')
lai1_an = np.load('Nebraska_LDAS/results/ekf/LAI_2012-01-01_2012-12-31.PData')
lai1_obs = np.load('Nebraska_LDAS/observations/sfx-trip/LAI_V2_2012-01-01_2012-12-31.PData')
'''
### Import entire folders of yearly data
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    return parts

### It is required to have as many 'read_pickle's as years in the folder 
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_LDAS/observations/sfx-trip/LAI_MC.t08*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15])]
lai_c4 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_LDAS/results/ol/LAI_20*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
lai = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_LDAS/results/ekf/LAI_20*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
lai_an = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_LDAS/observations/sfx-trip/LAI_V2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
lai_obs = pd.concat(df)

### Neb ECO SG



### Observations files contain up to the present, subtract the necessary number of observations here
lai_obs = lai_obs[:-4]

###################################################################
### Data Processing 
### Only takes data from the same time as observations 
###################################################################

tmp = lai.copy()*np.nan 
tmp.ix[lai_obs.index] = lai_obs 

### Use the following for graphing time-series
d_ol1 = lai[~np.isnan(tmp)].mean(axis=1).dropna()
d_an1 = lai_an[~np.isnan(tmp)].mean(axis=1).dropna()
d_ob1 = lai_obs[~np.isnan(tmp)].mean(axis=1).dropna()

### Use the following for mapping yearly results 
s_ol1 = lai[~np.isnan(tmp)].dropna()
s_an1 = lai_an[~np.isnan(tmp)].dropna()
s_ob1 = lai_obs[~np.isnan(tmp)].dropna()

'''
### Correlation Testing 
ol_cor = []
an_cor = []
for row in range(len(s_ol1.axes[1])):
    ### The '[0]' at the end only returns the first variable, correlation
    cor = sp.stats.pearsonr(s_ob1[row], s_ol1[row])[0]
    ol_cor.append(cor)
for row in range(len(s_an1.axes[1])):
    ### The '[0]' at the end only returns the first variable, correlation
    cor = sp.stats.pearsonr(s_ob1[row], s_an1[row])[0]
    an_cor.append(cor)

dif_cor = np.subtract(an_cor,ol_cor)
#sys.exit()
'''

lai_c4M = lai_c4.resample('M').mean()
lai_M = lai_obs.resample('M').mean()
lai_anM = d_an1.resample('M').mean()
lai_olM = d_ol1.resample('M').mean()


for ii in [1,2,3,4,5,6,7,8,9,10,11,12]:                                                        
    season = (lai_M.index.month == ii)
    lai_M.ix[season] = (lai_M.ix[season]-lai_M.ix[season].mean())/lai_M.ix[season].std()
    lai_anM.ix[season] = (lai_anM.ix[season]-lai_anM.ix[season].mean())/lai_anM.ix[season].std()
    lai_olM.ix[season] = (lai_olM.ix[season]-lai_olM.ix[season].mean())/lai_olM.ix[season].std()

for ii in [1,2,3,4,5,6,7,8,9,10,11,12]:
    season = (lai_c4M.index.month == ii)    
    lai_c4M.ix[season] = (lai_c4M.ix[season]-lai_c4M.ix[season].mean())/lai_c4M.ix[season].std()


###################################################################
### Graphing Functions 
### Only takes data from the same time as observations 
###################################################################
'''
### Single Graph Plot
plt.plot(lai_M.mean(axis=1),label='Monthly Anomally - All Obs')
plt.plot(lai_c4M.mean(axis=1),label='Monthly Anomally - C4 Obs')
plt.plot(lai_olM,label='Open Loop')
plt.plot(lai_anM,label='Analysis')


plt.legend(loc='upper left')
plt.title('LAI Monthly Anomally')
'''
### Multiple Graph Plots

f, axarr = plt.subplots(2, sharey = True)

axarr[0].plot(lai_M.mean(axis=1),label='Monthly Anomally - All Obs',linestyle="",marker="*"); axarr[0].plot(lai_c4M.mean(axis=1),label='Monthly Anomally - C4 Obs',linestyle="",marker="*",color="red"); axarr[0].plot(lai_olM,label='Open Loop') ; axarr[0].plot(lai_anM,label='Analysis')
axarr[0].set_title('All Observations')
axarr[0].legend(loc ='upper left')

axarr[1].plot(lai_c4M.mean(axis=1),label='Monthly Anomally - C4 Obs',linestyle="",marker="*")
axarr[1].plot(lai_M.mean(axis=1),label='Monthly Anomally - All Obs',linestyle="",marker="+",color="red")
axarr[1].set_title('C4 Observations')
axarr[1].legend(loc ='upper left')

plt.show()
raw_input("Press Enter to continue ...")

sys.exit()

###################################################################
### Final Plot
###################################################################

plt.show()

