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


def haversine(lon1,lat1,lon2,lat2,yy_lon,yy_lat):
    """
    calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal to radians
    lon1b, lat1b, lon2b, lat2b = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2b-lon1b
    dlat = lat2b-lat1b

    a = sin(dlat/2)**2 + cos(lat1b)*cos(lat2b)*sin(dlon/2)**2
    c = 2*asin(sqrt(a))
    r = 6371. # radius of Earth in km

    return c*r,lat2,lon2,yy_lat,yy_lon,lat1,lon1

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

def NS(s,o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
    s: simulated
    o: observed
    output:
    ns: Nash Sutcliffe efficient coefficient
    """
    s,o = filter_nan(s,o)
    return 1 - sum((s-o)**2)/sum((o-np.mean(o))**2)
    #return 1 - sum((np.log(s)-np.log(o))**2)/sum((np.log(o)-np.mean(np.log(o)))**2)

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


#fc = ('','fc1','fc2','fc3','fc4','fc5','fc6','fc7','fc8','fc9','fc10','fc11','fc12','fc13','fc14')
#fc = ('','fc1','fc2','fc3','fc4','fc5')
#fc = ('')
#for i in range(len(fc)):
    ## LOAD FILES
    #file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US15/results/ol{0}/WG4*.PData'.format(fc[i])),key=numericalSort)
    #df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
    #olwg4 = pd.concat(df)

    #file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US15/results/ekf{0}/WG4*.PData'.format(fc[i])),key=numericalSort)
    #df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
    #ekfwg4 = pd.concat(df)
    
    #file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US15/results/ol{0}/WG5*.PData'.format(fc[i])),key=numericalSort)
    #df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
    #olwg5 = pd.concat(df)

    #file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US15/results/ekf{0}/WG5*.PData'.format(fc[i])),key=numericalSort)
    #df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
    #ekfwg5 = pd.concat(df)
    
    #olwg_20 = ((olwg4*10)+(olwg5*20))/30
    #ekfwg_20 = ((ekfwg4*10)+(ekfwg5*20))/30
    
    #print("Saving 20cm layer PData {0}".format(fc[i]))
    #olwg_20.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US15/results/ol{0}/WG_20_2017_2018.PData'.format(fc[i]))
    #ekfwg_20.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US15/results/ekf{0}/WG_20_2017_2018.PData'.format(fc[i]))
    
    #print("OL WG4 mean : {0}".format(olwg4.mean(axis=1).mean()))
    #print("OL WG5 mean : {0}".format(olwg5.mean(axis=1).mean()))
    #print("OL 20cm mean: {0}".format(olwg_20.mean(axis=1).mean()))
    
### For non-FC Experiments
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/VOD/results/analysis_ssm_vod/WG4*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6])]
olwg4 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/VOD/results/analysis_vod/WG4*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6])]
ekfwg4 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/VOD/results/analysis_ssm_vod/WG5*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6])]
olwg5 = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/VOD/results/analysis_vod/WG5*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6])]
ekfwg5 = pd.concat(df)

olwg_20 = ((olwg4*10)+(olwg5*20))/30
ekfwg_20 = ((ekfwg4*10)+(ekfwg5*20))/30

print("Saving 20cm layer PData")
olwg_20.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/VOD/results/analysis_ssm_vod/WG_20_2017_2018.PData')
ekfwg_20.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/VOD/results/analysis_vod/WG_20_2017_2018.PData')    


print("Computation Complete")
sys.exit()
