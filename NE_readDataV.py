import pandas as pd
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import sys

forcing_dir = '/cnrm/vegeo/muciaa/NO_SAVE/Nebraska_ECOSG/forcings/'

### Select Years ###
begin_date = pd.to_datetime('20110101')
end_date = pd.to_datetime('20111231')

pd_period = pd.date_range(begin_date, end_date + pd.DateOffset(hours=24),freq='H') + pd.DateOffset(hours=9)
pd_day = pd.date_range(begin_date,end_date,freq='D')

filename = forcing_dir+'FORCING_'+begin_date.strftime("%Y%m%d%H")+'.nc'
ncfile = nc.Dataset(filename,'r')
nbtime, nblat, nblon = ncfile.variables['Rainf'][:].shape
ncfile.close()

var_period = pd.DataFrame(index=pd_period,columns=range(0,nblat*nblon))

for dat in pd_day:
    print(dat)
    
    pd_hours = pd.date_range(dat + pd.DateOffset(hours=9), dat + pd.DateOffset(hours=24+9),freq='H')
    
    filename = forcing_dir + 'FORCING_' + dat.strftime("%Y%m%d") + '09.nc'
    ncfile = nc.Dataset(filename,'r')
    
    var_aux = ncfile.variables['Rainf'][:].data
    lats = ncfile.variables['LAT'][:].data
    lons = ncfile.variables['LON'][:].data
    var_aux[var_aux > 10000.0] = np.nan
    
    #var_period.ix[pd_hours,:] = var_aux
    var_period.ix[pd_hours,:] = var_aux.reshape((25,nblat*nblon))
    #m = var_period.mean(axis=0).values.reshape(16,40) 
    #var_period.mean(axis=0).values.reshape(16,40)
    ncfile.close()

var_period.to_pickle('NE_2011_RF.PData')
sys.exit()

var_period = pd.read_pickle('NE_2012_RF.PData')

### Convert from Kg/m2/s to mm/month
var_period = var_period*3600
month = var_period.resample('M').sum()
### Convert from mm to inches
month = month/25.4


### Select specific pixels with Lincoln and Grand Island
linc_rain = month[153][:-1]
gi_rain = month[187][:-1]
### HPRCC monthly rainfall statistics for KLNK and KGRI
lin_2012 = [0.16,2.1,0.89,3.49,3,3.57,0.33,0.3,1.73,1.92,0.15,1.5]
lin_2009 = [0.38,0.64,0.18,1.52,1.17,6.18,1.84,3.2,1.25,4.24,0.06,2.42]
lin_2011 = [1.07,0.79,0.66,3.27,6.0,3.44,1.55,6.89,1.33,0.93,1.66,1.58]

gi_2012 = [0.16,1.05,0.83,1.41,2.29,1.27,0.16,0.94,0.47,0.78,0.53,1.66]
gi_2009 = [0.30,0.88,0.14,2.56,2.05,8.27,2.7,2.4,0.96,3.39,0.17,1.76]
gi_2011 = [1.49,0.28,1.02,2.93,8.7,1.72,4.37,2.3,0.83,2.2,0.21,1.11]

months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
x = np.arange(len(months))

width = 0.35   
#fig = plt.figure(figsize=(8,2.5),dpi=300,edgecolor='w')
fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
rects1 = ax.bar(x - width/2, linc_rain, width, label='ERA5',color='c')
rects2 = ax.bar(x + width/2, lin_2012, width, label='Station',color='b')

ax.set_title('Lincoln Precipitation : Station Obs vs ERA5 2012')
ax.set_ylabel('Precipitation [inches]')
#ax2.set_ylabel('Precipitation [inches]')
ax.set_xlabel('Months')
ax.set_xticks(x)
ax.set_xticklabels(months)
ax2.get_yaxis().set_visible(False)
fig.legend(loc=(0.85,0.85),fontsize=8)

plt.tight_layout()

plt.show()
