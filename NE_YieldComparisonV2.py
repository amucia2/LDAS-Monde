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

###################################################################
### Import Data  
###################################################################

### Import individual yearly files
lai_isba = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/Nebraska/LAI_ISBA.Pdata')
lai_cgls = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2/Nebraska/LAI_CGLS.Pdata')
vodx = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/Nebraska/VODX.Pdata')
vodc = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/Nebraska/VODC.Pdata')


### Annual Corn Yield (2000-2018) read as a CSV
corn_y = pd.read_csv('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/Nebraska/NE_Corn_Production_2000_2018_bu-acre.csv')['Value'][3:]
date = pd.date_range(start="2002-12-31",end='2017-12-31',freq='Y')
corn_y = pd.DataFrame(corn_y.values,index=date)

### Import entire folders of yearly data
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    return parts

###################################################################
### Data Processing 
###################################################################
### Taking Mean after resampling gives different values than if you take the mean before
### To be consistent, I am taking the mean first, then resampling
### Yearly Resampling - mean before
lai_cgls_A = lai_cgls.mean(axis=1).resample('A',label='left').mean()
lai_isba_A = lai_isba.mean(axis=1).resample('A',label='left').mean()
vodx_A = vodx.mean(axis=1).resample('A',label='left').mean()
vodc_A = vodc.mean(axis=1).resample('A',label='left').mean()

lai_cgls_A = lai_cgls.mean(axis=1).resample('A',label='left').max()
lai_isba_A = lai_isba.mean(axis=1).resample('A',label='left').max()
vodx_A = vodx.mean(axis=1).resample('A',label='left').max()
vodc_A = vodc.mean(axis=1).resample('A',label='left').max()

### Year Resampling - mean after
#lai_isba_A = lai_isba.resample('A',label='left').mean()
#lai_cgls_A = lai_cgls.resample('A',label='left').mean()
#vodx_A = vodx.resample('A',label='left').mean()
#vodc_A = vodc.resample('A',label='left').mean()

###################################################################
### Computer Z Scores
###################################################################
corn_y_Z = corn_y.apply(zscore)
corn_y_Z = pd.Series(corn_y_Z[0].values,index=pd.date_range(start="2002-12-31",end='2017-12-31',freq='Y'))

#lai_isba_Z = lai_isba_A.apply(zscore).mean(axis=1)
#lai_cgls_Z = lai_cgls_A.apply(zscore).mean(axis=1)
#vodx_Z = vodx_A.apply(zscore).mean(axis=1)
#vodc_Z = vodc_A.apply(zscore).mean(axis=1)

lai_isba_Z = pd.DataFrame(lai_isba_A).apply(zscore)
lai_cgls_Z = pd.DataFrame(lai_cgls_A).apply(zscore)
vodx_Z = pd.DataFrame(vodx_A).apply(zscore)
vodc_Z = pd.DataFrame(vodc_A).apply(zscore)


###################################################################
### Compute Correlations
###################################################################

cor_cgls_vodx = np.corrcoef(lai_cgls_Z[0],vodx_Z[0])
cor_cgls_vodc = np.corrcoef(lai_cgls_Z[0],vodc_Z[0])
cor_cgls_isba = np.corrcoef(lai_cgls_Z[0],lai_isba_Z[0])

cor_corn_cgls = np.corrcoef(corn_y_Z,lai_cgls_Z[0])
cor_corn_vodx = np.corrcoef(corn_y_Z,vodx_Z[0])
cor_corn_vodc = np.corrcoef(corn_y_Z,vodc_Z[0])
cor_corn_isba = np.corrcoef(corn_y_Z,lai_isba_Z[0])

print("Annual Observed LAI Anomaly vs Corn Yield Anomaly: " + str(round(cor_corn_cgls[0,1],4)))
print("------------------------------------------------------------------")
print("Annual ISBA LAI Anomaly vs Corn Yield Anomaly: " + str(round(cor_corn_isba[0,1],4)))
print("Annual VODX Anomaly vs Corn Yield Anomaly: " + str(round(cor_corn_vodx[0,1],4)))
print("Annual VODC Anomaly vs Corn Yield Anomaly: " +str(round(cor_corn_vodc[0,1],4)))
print("------------------------------------------------------------------")


###################################################################
### Graphing Functions 
###################################################################

### Debugging
'''
lai_iir = pd.read_pickle('Data/Nebraska/LAI_2000_2018_Nebraska.PData')['2003':'2017']['Obs']
lai_ii_A = lai_iir.resample('A',label='left').mean()
lai_ii = pd.DataFrame(lai_ii_A).apply(zscore)

lai_cgls = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2/Nebraska/LAI_CGLS.Pdata')['2003':'2017']
lai_cgls_A = lai_cgls.mean(axis=1).resample('A',label='left').mean()
lai_cgls_Z = pd.DataFrame(lai_cgls_A).apply(zscore)

lai_cgls = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2/Nebraska/LAI_CGLS.Pdata')
lai_cgls_A = lai_cgls.mean(axis=1).resample('A',label='left').mean()
lai_cgls_Zx = pd.DataFrame(lai_cgls_A).apply(zscore)

cor_corn_cgls = np.corrcoef(corn_y_Z,lai_cgls_Z[0])
cor_corn_cglsx = np.corrcoef(corn_y_Z,lai_cgls_Zx[0])

plt.plot(lai_iir['Obs'],color='black');plt.plot(lai_cgls.mean(axis=1),color='green');plt.show()
plt.plot(lai_iir['Obs'],color='black');plt.plot(lai_cgls_nan.mean(axis=1),color='green');plt.show()

plt.plot(lai_ii,color='black');plt.plot(lai_cgls_Z,color='green');plt.show()
'''

### Plot Yearly Results ###
#fig = plt.figure(figsize=(8,2),dpi=300,edgecolor='w')
#ax1 = fig.add_subplot(111)


plt.rcParams['figure.figsize'] = (15.5,5.5)
fig, ax1 = plt.subplots()
plt.title('Interannual Variability of LAI, VOD, and Corn Yield over Nebraska - Maximum',fontsize=16)
#'''
### ONLY YIELD AND OBS
ax1.plot(corn_y_Z,label='Corn Yield - USDA',marker='o',linestyle='--',markersize=8,color='green')
ax1.plot(lai_cgls_Z,label='LAI Observations',marker='*',linestyle='-',markersize=8,color='black')
ax1.plot(vodx_Z,label='VODX Observations',marker='^',linestyle='-',markersize=8,color='blue')
ax1.plot(vodc_Z,label='VODC Observations',marker='v',linestyle='-',markersize=8,color='red')
### Set Date Ticks and labels
import matplotlib.dates as mdates
years = mdates.YearLocator()  
years_fmt = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_locator(years)
ax1.xaxis.set_major_formatter(years_fmt)
datemin = np.datetime64(vodx_Z.index[0], 'Y')
datemax = np.datetime64(vodx_Z.index[-1], 'Y') + np.timedelta64(2, 'Y')
ax1.set_xlim(datemin, datemax)

ax1.format_xdata = mdates.DateFormatter('%Y-%m-%d')
ax1.format_ydata = lambda x: '$%1.2f' % x  # format the price.

fig.autofmt_xdate()

ax1.axhline(color='black',linewidth=1)

plt.legend(loc='upper left',numpoints=1)
plt.ylim([-2.5,2.5])
#plt.ylim([-1,1])
plt.ylabel('Anomaly',fontsize=16)

#plt.xlim('2002-01-31','2017-12-31')

#plt.rcParams.update({'font.size': 8})
fig.tight_layout()
#save('../vegeoFigures/Yield-RainF', ext='ps', close=True, verbose=True)
plt.show()



