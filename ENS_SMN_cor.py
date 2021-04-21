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



file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens*'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),\
    pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),pd.read_pickle(file[15]),\
    pd.read_pickle(file[16]),pd.read_pickle(file[17]),pd.read_pickle(file[18]),pd.read_pickle(file[19]),pd.read_pickle(file[20]),pd.read_pickle(file[21]),pd.read_pickle(file[22]),pd.read_pickle(file[23]),\
    pd.read_pickle(file[24]),pd.read_pickle(file[25]),pd.read_pickle(file[26]),pd.read_pickle(file[27]),pd.read_pickle(file[28]),pd.read_pickle(file[29]),pd.read_pickle(file[30]),pd.read_pickle(file[31]),\
    pd.read_pickle(file[32]),pd.read_pickle(file[33]),pd.read_pickle(file[34]),pd.read_pickle(file[35]),pd.read_pickle(file[36]),pd.read_pickle(file[37]),pd.read_pickle(file[38]),pd.read_pickle(file[39]),\
    pd.read_pickle(file[40]),pd.read_pickle(file[41]),pd.read_pickle(file[42]),pd.read_pickle(file[43]),pd.read_pickle(file[44]),pd.read_pickle(file[45]),pd.read_pickle(file[46]),pd.read_pickle(file[47]),\
    pd.read_pickle(file[48]),pd.read_pickle(file[49]),pd.read_pickle(file[50])]
R = pd.concat(df)

ens0 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens0.PData')
ens1 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens1.PData')
ens2 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens2.PData')
ens3 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens3.PData')
ens4 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens4.PData')
ens5 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens5.PData')
ens6 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens6.PData')
ens7 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens7.PData')
ens8 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens8.PData')
ens9 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens9.PData')
ens10 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens10.PData')

ens11 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens11.PData')
ens12 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens12.PData')
ens13 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens13.PData')
ens14 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens14.PData')
ens15 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens15.PData')
ens16 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens16.PData')
ens17 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens17.PData')
ens18 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens18.PData')
ens19 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens19.PData')
ens20 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens20.PData')

ens21 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens21.PData')
ens22 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens22.PData')
ens23 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens23.PData')
ens24 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens24.PData')
ens25 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens25.PData')
ens26 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens26.PData')
ens27 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens27.PData')
ens28 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens28.PData')
ens29 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens29.PData')
ens30 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens30.PData')

ens31 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens31.PData')
ens32 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens32.PData')
ens33 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens33.PData')
ens34 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens34.PData')
ens35 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens35.PData')
ens36 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens36.PData')
ens37 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens37.PData')
ens38 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens38.PData')
ens39 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens39.PData')
ens40 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens40.PData')

ens41 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens41.PData')
ens42 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens42.PData')
ens43 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens43.PData')
ens44 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens44.PData')
ens45 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens45.PData')
ens46 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens46.PData')
ens47 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens47.PData')
ens48 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens48.PData')
ens49 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens49.PData')
ens50 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens50.PData')

R_s = R[0]
Rd = R_s.values.reshape(51,27)
m_ens = list()
m_sta = list()

for i in range(len(Rd)):
    m_ens.append(Rd[i].mean())
    print 'ens'+ str(i) + ': ' + str(m_ens[i])
    
for i in range(len(Rd[0,:]-1)):
    m_sta.append(Rd[:,i].mean())
    print 'station'+ str(i) + ': ' + str(m_sta[i])



ticks = ['Ensemble Members','Stations']

f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
plt.title('Correlation Variability',fontsize=12)
ax1.boxplot(m_ens)
ax1.set_title('Ensemble Members')
ax2.boxplot(m_sta)
ax2.set_title('Stations')
plt.ylim(0.3, 0.9)
plt.tight_layout()
plt.show()

    




