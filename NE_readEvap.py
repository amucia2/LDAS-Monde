import pandas as pd
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import sys

forcing_dir = '/cnrm/vegeo/bonanb/LDAS/ALEXI/ET_04/CONUS/'

### Select Years ###
begin_date = pd.to_datetime('20100101')
end_date = pd.to_datetime('20181231')

pd_period = pd.date_range(begin_date, end_date,freq='D')
pd_day = pd.date_range(begin_date,end_date,freq='D')

filename = forcing_dir+'CONUS_ET_'+begin_date.strftime("%Y%m%d")+'.nc'
ncfile = nc.Dataset(filename,'r')
nblat, nblon = ncfile.variables['ET'][:].shape
ncfile.close()

var_period = pd.DataFrame(index=pd_period,columns=range(0,4*12))

for dat in pd_day:
    print(dat)    
    
    filename = forcing_dir + 'CONUS_ET_' + dat.strftime("%Y%m%d") + '.nc'
    ncfile = nc.Dataset(filename,'r')
    
    lats = ncfile.dimension['lat'][:].data
    lons = ncfile.dimension['lon'][:].data
    
    ### Locations only within the Subregion of Nebraska
    if (lons[dat] >=-100.0 and lons[dat] <=-97.0) and (lats[dat] >=40.25 and lats[dat] <=41.25):
        var_aux = ncfile.variables['ET'][:].data
        var_aux[var_aux > 10000.0] = np.nan
    
        var_period.ix[:] = var_aux
        #m = var_period.mean(axis=0).values.reshape(16,40) 
        #var_period.mean(axis=0).values.reshape(16,40)
        ncfile.close()

var_period.to_pickle('EVAP_CONUS.PData')
sys.exit()

