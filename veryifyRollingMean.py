import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#from mpl_toolkits.basemap import Basemap
import pandas as pd
#from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr
import pylab as pl
import sys, glob, os, re

'''
Model : LAI ISBA
Analysis : LAI CGLS
Obs : VOD
'''

#A = pd.read_pickle('LAIfromVODCAX_V4.PData') 

#B = pd.rolling_mean(A,window=30,min_periods=5,center=True)
#B = A.rolling(window=30,min_periods=5,center=True).mean()


rollVOD = pd.read_pickle('LAIfromVODCAX_V5_rollingMean_2000.PData')
VOD = pd.read_pickle('MatchedVOD/LAI_2000-01-01_2000-12-31.PData')
obs = pd.read_pickle('observations/LAI_V2_2000-01-01_2000-12-31.PData')

rVOD = VOD.rolling(window=30,min_periods=5,center=True).mean()
r2VOD = rVOD.rolling(window=30,min_periods=5,center=True).mean()

plt.plot(rollVOD.mean(axis=1),color='r',label='VOD Rolling Mean - newly produced')
plt.plot(VOD.mean(axis=1),color='b',label='VOD Assimilated')
plt.plot(rVOD.mean(axis=1),color='b',linestyle='--',label='VOD Rolling applied')
#plt.plot(r2VOD.mean(axis=1),color='r',linestyle='--',label='VOD Rolling applied twice')
plt.plot(obs.mean(axis=1),color='black',label='LAI Observations')
plt.legend()
plt.show()



vod_roll = pd.read_pickle('LAIfromVODCAX_TS.PData')

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts


file = sorted(glob.glob('observations/LAI_V2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),pd.read_pickle(file[18])]
lai_cgls_ = pd.concat(df)

date = pd.date_range(start="2003-01-01 09:00:00",end='2018-12-31 09:00:00',freq='D')
lai_cgls = lai_cgls_.reindex(date,fill_value=np.nan)

laifromVOD = vod_roll['2003-01-01':'2018-12-31']
vod2012 = laifromVOD['2012-01-01':'2012-12-31']
obs2012 = lai_cgls_['2012-01-01':'2012-12-31']

plt.title("LAI from VOD (rolling mean) over CONUS : 2003-2018")
plt.ylabel('LAI [m^2/m^2]')
plt.plot(laifromVOD,color='red',label='LAI from VOD w/ rolling Mean')
plt.plot(lai_cgls_.mean(axis=1),color='green',label='CGLS LAI Observations')
plt.legend()
plt.show()

plt.title("LAI from VOD (rolling mean) over CONUS : 2012")
plt.ylabel('LAI [m^2/m^2]')
plt.plot(vod2012,color='red',label='LAI from VOD w/ rolling Mean')
plt.plot(obs2012.mean(axis=1),color='green',label='CGLS LAI Observations')
plt.legend()
plt.show()

### Take only values where LAI observations are found (both time and space)

lai_vodx_ = lai_vodx[~np.isnan(lai_cgls_)]

lai_cgls = lai_cgls_.mean(axis=1)

lai_vodx_old = pd.read_pickle('LAIfromVODCAX_TS.PData')
lai_vodx = pd.read_pickle('LAIfromVODX_V8_TS_whereObs.PData').dropna()
lai_vodc = pd.read_pickle('LAIfromVODC_V8_TS_whereObs.PData').dropna()
lai_vodx_int = pd.read_pickle('LAIfromVODCAX_V7_2003-2018_interpolated_TSwhereObs.PData').dropna()
#mask = lai_vodx[0].isnan().groupby(lai_vodx.index.normalize()).transform('any')


lai_vodx_old = lai_vodx_old['2003-01-01':'2018-12-31']
lai_vodx_old = lai_vodx_old[~np.isnan(lai_cgls)].dropna()
lai_cgls = lai_cgls.dropna()

plt.title("LAI from VOD over CONUS : 2003-2018")
plt.ylabel('LAI [m^2 m^-2]')
plt.plot(lai_vodx,color='red',label='Matched VODX')
#plt.plot(lai_vodc,color='blue',label='Matched VODC')
plt.plot(lai_vodx_int,color='blue',linestyle='--',label='Matched VODX with Interpolated LAI Obs')
plt.plot(lai_cgls,color='green',label='CGLS LAI Observations')
plt.legend()
plt.show()

lai_vodx2012 = lai_vodx['2012-01-01':'2012-12-31']
lai_vodxold2012 = lai_vodx_old['2012-01-01':'2012-12-31']
lai_vodc2012 = lai_vodc['2012-01-01':'2012-12-31']
obs2012 = lai_cgls['2012-01-01':'2012-12-31']
lai_vodx_int2012 = lai_vodx_int['2012-01-01':'2012-12-31']

plt.title("LAI from VOD over CONUS : 2012")
plt.ylabel('LAI [m^2 m^-2]')
plt.plot(lai_vodx2012,color='red',label='Matched VODX')
#plt.plot(lai_vodc2012,color='blue',label='Matched VODC')
#plt.plot(lai_vodxold2012,color='red',linestyle='--',label='OLD Matched VODX')
plt.plot(lai_vodx_int2012,color='blue',linestyle='--',label='Matched VODX with Interpolated LAI Obs')
plt.plot(obs2012,color='green',label='CGLS LAI Observations')
plt.legend()
plt.show()