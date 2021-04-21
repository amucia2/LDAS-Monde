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
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]
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

'''
### Read in data from original Pickles
var = input("Select variable (WG2-8, LAI, EVAP)? : ")

### Import files from all 6 exps
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/ol/{0}*.PData'.format(var)),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
ol = pd.concat(df)


file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/ekf_swi/{0}*.PData'.format(var)),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
ekf_ssm = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/ekf_lai/{0}*.PData'.format(var)),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
ekf_lai = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/ekf_lai_swi/{0}*.PData'.format(var)),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
ekf_laissm = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/ekf_vod/{0}*.PData'.format(var)),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
ekf_vod = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/ekf_vod_swi/{0}*.PData'.format(var)),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
ekf_vodssm = pd.concat(df)

### Or saved in concatenated version (same size)
ol = pd.read_pickle('Data/US00_concat/ol*{0}*.PData'.format(var))
ekf_lai = pd.read_pickle('Data/US00_concat/ekf_lai*{0}*.PData'.format(var))
ekf_vod = pd.read_pickle('Data/US00_concat/ekf_vod*{0}*.PData'.format(var))
ekf_vodssm = pd.read_pickle('Data/US00_concat/ekf_vodssm*{0}*.PData'.format(var))

if var == 'EVAP':
    ob = 'ET_04'
elif var == 'LAI':
    ob = 'LAI_V2'
elif var == 'WG2':
    ob = 'SSM'

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/sfx-trip/{0}*.PData'.format(ob)),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
obs = pd.concat(df)

'''
### Option to open the daily mean files (far smaller)
ol_m = pd.read_pickle('Data/US00_concat/TS_ol_LAI.PData')
#ekf_ssm_m = pd.read_pickle('Data/US00_concat/TS_ekfssm_LAI.PData')
ekf_lai_m = pd.read_pickle('Data/US00_concat/TS_ekflai_LAI.PData')
#ekf_laissm_m = pd.read_pickle('Data/US00_concat/TS_ekflaissm_LAI.PData')
ekf_vod_m = pd.read_pickle('Data/US00_concat/TS_ekfvod_LAI.PData')
ekf_vodssm_m = pd.read_pickle('Data/US00_concat/TS_ekfvodssm_LAI.PData')

obs_m = pd.read_pickle('Data/US00_concat/TS_obs_LAI.PData')
vod = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/LAIfromVODCAX_V4.PData')
vod_m = vod.mean(axis=1)

ol_m.index = ol_m.index.rename('date')
#ekf_ssm_m.index = ekf_ssm_m.index.rename('date')
ekf_lai_m.index = ekf_lai_m.index.rename('date')
#ekf_laissm_m.index = ekf_laissm_m.index.rename('date')
ekf_vod_m.index = ekf_vod_m.index.rename('date')
ekf_vodssm_m.index = ekf_vodssm_m.index.rename('date')

#obs.index = obs.index.normalize()
#obs = obs.mean(axis=1)

#begin_date = pd.to_datetime('20000101090000')
#end_date = pd.to_datetime('20181231090000')

#pd_period = pd.date_range(begin_date, end_date,freq='D')
#obs = obs.reindex(pd_period)


#ol_m = ol.mean(axis=1)
#ekf_ssm_m = ekf_ssm.mean(axis=1)
#ekf_lai_m = ekf_lai.mean(axis=1)
#ekf_laissm_m = ekf_laissm.mean(axis=1)
#ekf_vod_m = ekf_vod.mean(axis=1)
#ekf_vodssm_m = ekf_vodssm.mean(axis=1)

sys.exit()
'''
year = input("Select Year to Display : ")
    
### Plot time series

fig = plt.figure(figsize=(8,2.5),dpi=300,edgecolor='w')
ax1 = fig.add_subplot(111)
plt.title('Assimilation Scheme Differences',fontsize=7)
ax1.set_ylabel('{0}'.format(var))
ax1.set_xlabel('Date')

l1 = ax1.plot(ol_m[year],color='b',linewidth=1,label='ol',linestyle='--')[0]
#l2 = ax1.plot(ekf_ssm_m[year],color='orange',linewidth=1,label='ekf_ssm',linestyle='--')[0]
l3 = ax1.plot(ekf_lai_m[year],color='cyan',linewidth=1,label='ekf_lai',linestyle='--')[0]
#l4 = ax1.plot(ekf_laissm_m[year],color='r',linewidth=1,label='ekf_laissm',linestyle='--')[0]
l5 = ax1.plot(ekf_vod_m[year],color='indigo',linewidth=1,label='ekf_vod',linestyle='--')[0]
l6 = ax1.plot(ekf_vodssm_m[year],color='m',linewidth=1,label='ekf_vodssm',linestyle='--')[0]

l7 = ax1.plot(vod_m[year],color='lime',linewidth=1,label='VOD Obs')[0]
#l8 = ax1.plot(obs_m[year],color='green',linewidth=1,label='LAI Obs')[0]

plt.legend()
plt.show()
'''

### Plot Annual Time Series
### Averaging each month/day over the 19 years to get a single year 
'''
begin_date = pd.to_datetime('20000101090000')
end_date = pd.to_datetime('20181231090000')

pd_period = pd.date_range(begin_date, end_date,freq='D')
obs = obs_m.reindex(pd_period)
obs_A = obs_m.groupby([obs_m.index.month,obs_m.index.day],dropna=False).mean()
'''

ol_A = ol_m.groupby([ol_m.index.month,ol_m.index.day]).mean()
#ekf_ssm_A = ekf_ssm_m.groupby([ekf_ssm_m.index.month,ekf_ssm_m.index.day]).mean()
ekf_lai_A = ekf_lai_m.groupby([ekf_lai_m.index.month,ekf_lai_m.index.day]).mean()
#ekf_laissm_A = ekf_laissm_m.groupby([ekf_laissm_m.index.month,ekf_laissm_m.index.day]).mean()
ekf_vod_A = ekf_vod_m.groupby([ekf_vod_m.index.month,ekf_vod_m.index.day]).mean()
ekf_vodssm_A = ekf_vodssm_m.groupby([ekf_vodssm_m.index.month,ekf_vodssm_m.index.day]).mean()
vod_A = vod_m.groupby([vod_m.index.month,vod_m.index.day]).mean()
obs_A = obs_m.groupby([obs_m.index.month,obs_m.index.day]).mean()

ol_A = ol_A.reset_index().iloc[:,2]
#ekf_ssm_A = ekf_ssm_A.reset_index().iloc[:,2]
ekf_lai_A = ekf_lai_A.reset_index().iloc[:,2]
#ekf_laissm_A = ekf_laissm_A.reset_index().iloc[:,2]
ekf_vod_A = ekf_vod_A.reset_index().iloc[:,2]
ekf_vodssm_A = ekf_vodssm_A.reset_index().iloc[:,2]
vod_A = vod_A.reset_index().iloc[:,2]
obs_A = obs_A.reset_index().iloc[:,2]

ol_r = ol_A.rolling(7).mean()
#ekf_ssm_r = ekf_ssm_A.rolling(7).mean()
ekf_lai_r = ekf_lai_A.rolling(7).mean()
#ekf_laissm_r = ekf_laissm_A.rolling(7).mean()
ekf_vod_r = ekf_vod_A.rolling(7).mean()
ekf_vodssm_r = ekf_vodssm_A.rolling(7).mean()
vod_r = vod_A.rolling(7).mean()
obs_r = obs_A.rolling(7).mean()

fig = plt.figure(figsize=(8,2.5),dpi=300,edgecolor='w')
ax1 = fig.add_subplot(111)
plt.title('Assimilation Scheme Differences',fontsize=7)
ax1.set_ylabel('{0}'.format(var))
ax1.set_xlabel('Date')

l1 = ax1.plot(ol_A,color='b',linewidth=1,label='ol',linestyle='--')[0]
#l2 = ax1.plot(ekf_ssm_A,color='orange',linewidth=1,label='ekf_ssm',linestyle='--')[0]
l3 = ax1.plot(ekf_lai_A,color='cyan',linewidth=1,label='ekf_lai',linestyle='--')[0]
#l4 = ax1.plot(ekf_laissm_A,color='r',linewidth=1,label='ekf_laissm',linestyle='--')[0]
l5 = ax1.plot(ekf_vod_A,color='indigo',linewidth=1,label='ekf_vod',linestyle='--')[0]
l6 = ax1.plot(ekf_vodssm_A,color='m',linewidth=1,label='ekf_vodssm',linestyle='--')[0]

l7 = ax1.plot(vod_A,color='lime',linewidth=1,label='VOD Obs')[0]
l8 = ax1.plot(obs_A,color='green',linewidth=1,label='{0} Obs'.format(ob))[0]


plt.legend()
plt.show()

### Plot Annual Time Series with Rolling Mean
fig = plt.figure(figsize=(8,2.5),dpi=300,edgecolor='w')
#fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.title('Assimilation Scheme Differences',fontsize=7)
ax1.set_ylabel('{0} [m2/m2]'.format(var))
ax1.set_xlabel('Day of Year')

l1 = ax1.plot(ol_r,color='b',linewidth=1,label='ol',linestyle='--')[0]
#l2 = ax1.plot(ekf_ssm_r,color='orange',linewidth=1,label='ekf_ssm',linestyle='--')[0]
l3 = ax1.plot(ekf_lai_r,color='cyan',linewidth=1,label='ekf_lai',linestyle='--')[0]
#l4 = ax1.plot(ekf_laissm_r,color='r',linewidth=1,label='ekf_laissm',linestyle='--')[0]
l5 = ax1.plot(ekf_vod_r,color='indigo',linewidth=1,label='ekf_vod',linestyle='--')[0]
l6 = ax1.plot(ekf_vodssm_r,color='m',linewidth=1,label='ekf_vodssm',linestyle='--')[0]
l7 = ax1.plot(vod_r,color='lime',linewidth=1,label='VOD Obs')[0]
#l8 = ax1.plot(obs_r,color='green',linewidth=1,label='{0} Obs'.format(ob))[0]

plt.legend()
plt.show()
