import numpy as np
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sys, glob, os, re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import pylab as pl
import pandas as pd
from netCDF4 import Dataset
import scipy.stats
from datetime import datetime
from datetime import date, timedelta



### Import files of interest
f = '/cnrm/vegeo/muciaa/NO_SAVE/US_FC15/time_series/OL/all/USCRN_LDAS_OL_48.3100_-105.100_WG3.DAT'

#d = np.loadtxt(f,dtype=str)[:,0]
#h = np.loadtxt(f,dtype=str)[:,1]

### Other way to load data

data = pd.read_csv(f,sep='\t',header=None)
t = data[0].str.split(' ',1,expand=True)
t[1] = t[1].replace('24','23')
t[1] = t[1].str.zfill(2)
t[2] = t[0]+' '+t[1]
t[2] = pd.to_datetime(t[2])

### Removes Hours Minutes Seconds

t[2] = pd.DatetimeIndex(t[2]).normalize()

#dt_f = list()
#date = list()

#for i in range(len(data)):
#    if h[i] == '24':
#        h[i] = '23'
#    date.append(d[i] + ' ' + h[i].zfill(2))
#    dt_f.append(datetime.strptime(t.at[i,2],'%Y-%m-%d %H'))

#date_set = set(t.at[0,2] + timedelta(x) for x in range((t.at[-1,2]-t.at[0,2]).days))
#date_set = set(dt_f[0] + timedelta(x) for x in range((dt_f[-1]-dt_f[0]).days))

#missing = sorted(date_set - set(t[2]))

missing = pd.date_range(start='2017-01-01',end='2018-12-31').difference(t[2])

print missing

