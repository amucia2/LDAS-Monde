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
    s = s.flatten()
    o = o.flatten()
    s = np.transpose(s)
    o = np.transpose(o)
    s[np.invert(~np.isnan(o))] = np.nan
    o[np.invert(~np.isnan(s))] = np.nan
    s = s[~pd.isnull(s)]
    o = o[~pd.isnull(o)]
    return s,o
    

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

# LOAD EXTRA FILES
### CONUS Region

nc=Dataset('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/sfx-trip/pgd_prep/PREP.nc','r')
lat_1D=nc.variables['REG_LAT'][:]
lon_1D=nc.variables['REG_LON'][:]
#nc.close()


depth = ('0.05','0.20','0.50','1.00')

exp1_list = ('ol','ekf_ssm','ekf_lai_ssm')
exp2_list = ('ekf_lai','ekf_vod','ekf_vod_ssm')


lat = np.loadtxt('/cnrm/vegeo/muciaa/NO_SAVE/US00/USCRN_coords_wban.txt')[:,0]
lon = np.loadtxt('/cnrm/vegeo/muciaa/NO_SAVE/US00/USCRN_coords_wban.txt')[:,1]
wban = np.genfromtxt('/cnrm/vegeo/muciaa/NO_SAVE/US00/USCRN_coords_wban.txt',dtype='str')[:,2]

count = 1

for m in range(len(exp1_list)):
    exp1 = exp1_list[m]
    exp2 = exp2_list[m]
    
    for i in range(len(depth)):
        dep = depth[i]
        if dep == '0.05':
            wg ='WG3'
            smcol = 18
            stcol = 23
        elif dep == '0.20':
            wg = 'WG_20'
            smcol = 20
            stcol = 25
        elif dep == '0.50':
            wg = 'WG6'
            smcol = 21
            stcol = 26
        elif dep == '1.00':
            wg = 'WG8'
            smcol = 22
            stcol = 27
        
        print("Experiment 1 : {0}".format(exp1))
        print("Experiment 2 : {0}".format(exp2))
        print("Soil Layer : {0}".format(wg))
        print("Depth : {0}".format(dep))

        file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/{0}/{1}_201*.PData'.format(exp1,wg)),key=numericalSort)
        df = [pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8])]
        wg3_1D = pd.concat(df)

        file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/{0}/{1}_201*.PData'.format(exp2,wg)),key=numericalSort)
        df = [pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8])]
        wg3_ekf_1D = pd.concat(df)

        #file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/sfx-trip/SSM_COM*_201*.PData'),key=numericalSort)
        #df = [pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),pd.read_pickle(file[8])]
        #wg3_obs_1D = pd.concat(df)

        lat_2D = lat_1D
        lon_2D = lon_1D
        wg3_2D = wg3_1D
        wg3_ekf_2D = wg3_ekf_1D

        rng3 = pd.date_range('01/01/2011',periods=wg3_1D.shape[0],freq='D')
        rng4 = pd.date_range('01/01/2011',periods=wg3_ekf_1D.shape[0],freq='D')

        lat_keep = list()
        lon_keep = list()

        R_mod = list() ; bias_mod = list() ; ubrmse_mod = list() ; P_mod = list() ; E_mod = list()
        R_ekf = list() ; bias_ekf = list() ; ubrmse_ekf = list() ; P_ekf = list() ; E_ekf = list()
        R_wg3ol = list() ; R_wg6ol = list(); R_wg8ol = list()
        R_wg3ekf = list() ; R_wg6ekf = list(); R_wg8ekf = list()
        
        R_win_ol = list() ; R_spr_ol = list() ; R_sum_ol = list() ; R_aut_ol = list()
        R_win_ekf = list() ; R_spr_ekf = list() ; R_sum_ekf = list() ; R_aut_ekf = list()
        R_wg3ol = list() ; R_wg6ol = list(); R_wg8ol = list()
        R_wg3ekf = list() ; R_wg6ekf = list(); R_wg8ekf = list()
        R_grow_ol = list(); R_grow_ekf = list()

        

        for coord in range(len(lat)):
            print(coord)
            if (lat_1D[-1,0] >= lat[coord] >= lat_1D[0,0]) and (lon_1D[0,-1] >= lon[coord] >= lon_1D[0,0]):
                
                ### Load USCRN Soil Moisture (and Temperature) Observations   
                fn = glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/USCRN/concat/{0}*'.format(wban[coord]))
                obs_05 = np.genfromtxt(fn[0])[:,smcol]
                obs_05ts = np.genfromtxt(fn[0])[:,stcol]
                
                rng = pd.date_range('2011-01-01 00:00:00', periods=len(obs_05), freq='D')
                
                obs_05bis = pd.Series(obs_05, index = rng)
                obs_05tsbis = pd.Series(obs_05ts, index = rng)
                
                ### Filters Nan and frozen conditions from observations
                obs_05bis[obs_05bis<0.] = np.nan
                obs_05bis[obs_05tsbis<4.] = np.nan
                
                ### more than 100 days of obs check
                if (len(obs_05bis.resample('D').mean()[~np.isnan(obs_05bis.resample('D').mean())]) > 100):

                    dst = list()
                    for y_lat in range(len(lat_1D[:,0])):
                        for y_lon in range(len(lon_1D[0,:])):
                            dst.append(haversine(lon[coord],lat[coord],lon_1D[0,:][y_lon],lat_1D[:,0][y_lat],y_lon,y_lat))
                    dst = np.array(dst)
                    
                    print(dst[[np.argsort(dst[:,0])[0]],0],\
                            lat_1D[[int(dst[[np.argsort(dst[:,0])[0]],3][0])],0],\
                            lon_1D[0,[int(dst[[np.argsort(dst[:,0])[0]],4][0])]])

                    SSM = list() ; SSM2 = list()
                    ### Was reshape(52,80) for MIDW region, for US_FC 175,350, US_FC_IFS is 350,700, VOD 140,280
                    for l in range(wg3_2D.shape[0]):
                        SSM.append(wg3_2D.values[l].reshape(140,280)[int(dst[[np.argsort(dst[:,0])[0]],3][0]),int(dst[[np.argsort(dst[:,0])[0]],4][0])])
                    for l in range(wg3_ekf_2D.shape[0]):
                        SSM2.append(wg3_ekf_2D.values[l].reshape(140,280)[int(dst[[np.argsort(dst[:,0])[0]],3][0]),int(dst[[np.argsort(dst[:,0])[0]],4][0])])
                        
    
                    SSM_bis = pd.Series(np.array(SSM),index=rng3)
                    SSM2_bis = pd.Series(np.array(SSM2),index=rng4)
                    wg3_bis = SSM_bis
                    wg3_ekf_bis = SSM2_bis
                    
                    grow_ol = wg3_bis.loc[(wg3_bis.index.month==4)+(wg3_bis.index.month==5)+(wg3_bis.index.month==6)+(wg3_bis.index.month==7)+(wg3_bis.index.month==8)+(wg3_bis.index.month==9)]
                    grow_ekf = wg3_ekf_bis.loc[(wg3_ekf_bis.index.month==4)+(wg3_ekf_bis.index.month==5)+(wg3_ekf_bis.index.month==6)+(wg3_ekf_bis.index.month==7)+(wg3_ekf_bis.index.month==8)+(wg3_ekf_bis.index.month==9)]

                    win_ol = wg3_bis.loc[(wg3_bis.index.month==12)+(wg3_bis.index.month==1)+(wg3_bis.index.month==2)]
                    win_ekf = wg3_ekf_bis.loc[(wg3_ekf_bis.index.month==12)+(wg3_ekf_bis.index.month==1)+(wg3_ekf_bis.index.month==2)]
                    spr_ol = wg3_bis.loc[(wg3_bis.index.month==3)+(wg3_bis.index.month==4)+(wg3_bis.index.month==5)]
                    spr_ekf = wg3_ekf_bis.loc[(wg3_ekf_bis.index.month==3)+(wg3_ekf_bis.index.month==4)+(wg3_ekf_bis.index.month==5)]
                    sum_ol = wg3_bis.loc[(wg3_bis.index.month==6)+(wg3_bis.index.month==7)+(wg3_bis.index.month==8)]
                    sum_ekf = wg3_ekf_bis.loc[(wg3_ekf_bis.index.month==6)+(wg3_ekf_bis.index.month==7)+(wg3_ekf_bis.index.month==8)]
                    aut_ol = wg3_bis.loc[(wg3_bis.index.month==9)+(wg3_bis.index.month==10)+(wg3_bis.index.month==11)]
                    aut_ekf = wg3_ekf_bis.loc[(wg3_ekf_bis.index.month==9)+(wg3_ekf_bis.index.month==10)+(wg3_ekf_bis.index.month==11)]

                    
                    if len(wg3_ekf_bis[~np.isnan(wg3_ekf_bis)])>100.:
                        #obs_05bis_DD = obs_05bis.resample('D').mean()
                        obs_05bis_DD = obs_05bis["2011-01-01":"2018-12-31"] 
                        lat_keep.append(lat[coord]) ; lon_keep.append(lon[coord])
                
                        print("R ol : %.4f Bias ol : %.4f ubRMSE ol : %.4f"\
                            %(correlation((wg3_bis).values,\
                        obs_05bis_DD.values)[0],\
                        bias((wg3_bis).values,\
                        obs_05bis_DD.values),\
                        ubrmse((wg3_bis).values,\
                        obs_05bis_DD.values)))

                        print("R ekf : %.4f Bias ekf : %.4f ubRMSE ekf : %.4f"\
                            %(correlation((wg3_ekf_bis).values,\
                        obs_05bis_DD.values)[0],\
                        bias((wg3_ekf_bis).values,\
                        obs_05bis_DD.values),\
                        ubrmse((wg3_ekf_bis).values,\
                        obs_05bis_DD.values)))
                    ##else:
                        #obs_05bis_DD = np.nan
                    
                        R_mod.append(correlation((wg3_bis).values,obs_05bis_DD.values)[0])
                        P_mod.append(correlation((wg3_bis).values,obs_05bis_DD.values)[1])
                        bias_mod.append(bias((wg3_bis).values,obs_05bis_DD.values))
                        ubrmse_mod.append(ubrmse((wg3_bis).values,obs_05bis_DD.values))

                        R_ekf.append(correlation((wg3_ekf_bis).values,obs_05bis_DD.values)[0])
                        P_ekf.append(correlation((wg3_ekf_bis).values,obs_05bis_DD.values)[1])
                        bias_ekf.append(bias((wg3_ekf_bis).values,obs_05bis_DD.values))
                        ubrmse_ekf.append(ubrmse((wg3_ekf_bis).values,obs_05bis_DD.values))
                    else:
                        print("Not >100 LDAS Days")
                        R_mod.append(np.nan)
                        P_mod.append(np.nan)
                        bias_mod.append(np.nan)
                        ubrmse_mod.append(np.nan)
                        
                        R_ekf.append(np.nan)
                        P_ekf.append(np.nan)
                        bias_ekf.append(np.nan)
                        ubrmse_ekf.append(np.nan)
                else:
                    print("Not >100 Observations")
                    R_mod.append(np.nan)
                    P_mod.append(np.nan)
                    bias_mod.append(np.nan)
                    ubrmse_mod.append(np.nan)
                    
                    R_ekf.append(np.nan)
                    P_ekf.append(np.nan)
                    bias_ekf.append(np.nan)
                    ubrmse_ekf.append(np.nan)
            else:
                print("Lat/Lon outside of domain")
                R_mod.append(np.nan)
                P_mod.append(np.nan)
                bias_mod.append(np.nan)
                ubrmse_mod.append(np.nan)
                
                R_ekf.append(np.nan)
                P_ekf.append(np.nan)
                bias_ekf.append(np.nan)
                ubrmse_ekf.append(np.nan)

        R_mod = pd.DataFrame(R_mod)
        R_ekf = pd.DataFrame(R_ekf)
        P_mod = pd.DataFrame(P_mod)
        P_ekf = pd.DataFrame(P_ekf)
        E_mod = pd.DataFrame(ubrmse_mod)
        E_ekf = pd.DataFrame(ubrmse_ekf)

        #sys.exit()
        print('Saving Files in /cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/')
        R_mod.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/{0}_{1}_{2}_R.PData'.format(exp1,wg,dep))
        R_ekf.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/{0}_{1}_{2}_R.PData'.format(exp2,wg,dep))

        P_mod.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/{0}_{1}_{2}_P.PData'.format(exp1,wg,dep))
        P_ekf.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/{0}_{1}_{2}_P.PData'.format(exp2,wg,dep))

        E_mod.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/{0}_{1}_{2}_E.PData'.format(exp1,wg,dep))
        E_ekf.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/{0}_{1}_{2}_E.PData'.format(exp2,wg,dep))

print("All Correlations Computed and Saved")
sys.exit()

######################################
w = np.where(P_mod<0.05 and P_ekf<0.05)
fig, axes = plt.subplots(1,1)
axes.set_title("Evaluation vs USCRN")
map = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-55, urcrnrlon=180, urcrnrlat=90,resolution='c',ax = axes)
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines()
map.drawmapboundary(fill_color='lightblue')
map.drawcountries()
map.fillcontinents(color='beige',lake_color='lightblue')
x,y = map(np.array(lon_keep), np.array(lat_keep))
cs2 = map.scatter(x,y,marker='s',s=30,c=(np.array(R_ekf)-np.array(R_mod)),vmin=-0.1,vmax=0.1,cmap=cm,zorder=2)
cbar = map.colorbar(cs2,location='bottom',pad="5%")
plt.show()

from matplotlib import ticker
NIC = 100*(np.array(R_ekf)-np.array(R_mod))/(1-np.array(R_mod))

w1=np.where(NIC>3)[0]
w3=np.where((NIC>-3) & (NIC<3))[0]
w2=np.where(NIC<-3)[0]

plt.figure(figsize=(8,6),dpi=300,edgecolor='w')
plt.title("NIC R NoahMP vs LDAS-Monde")
map = Basemap(projection='aea',llcrnrlon=-122, llcrnrlat=22, urcrnrlon=-62, urcrnrlat=48,lat_0=40,lon_0=-98,resolution='l')
map.shadedrelief()
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.8)
map.drawcountries(1.8)
#map.drawmapscale(lon=-90,lat=27,lon0=-98,lat0=40,length=750,barstyle='fancy')
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
#x,y = map(np.array(lon_keep)[w3], np.array(lat_keep)[w3])
#cs2 = map.scatter(x,y,marker='D',s=10,c=NIC[w3],vmin=-20,vmax=20,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w3], np.array(lat_keep)[w3])
cs2 = map.scatter(x,y,marker='v',s=25,c=NIC[w3],vmin=-50,vmax=50,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w2], np.array(lat_keep)[w2])
cs2 = map.scatter(x,y,marker='o',s=45,c=NIC[w2],vmin=-50,vmax=50,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w1], np.array(lat_keep)[w1])
cs2 = map.scatter(x,y,marker='o',s=45,c=NIC[w1],vmin=-50,vmax=50,cmap=cm,zorder=2)
cbar = map.colorbar(cs2,location='bottom',pad="10%")
cbar.set_label("NIC R")
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
plt.rcParams.update({'font.size': 12})
#plt.savefig('images/SMN_NIC_NoahMP_monde.png',format='png',dpi=300)
plt.show()


plt.figure(figsize=(8,6),dpi=300,edgecolor='w')
plt.title("NIC R NoahMP vs LDAS-Monde")
map = Basemap(projection='aea',llcrnrlon=-122, llcrnrlat=22, urcrnrlon=-62, urcrnrlat=48,lat_0=40,lon_0=-98,resolution='l')
map.shadedrelief()
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.8)
map.drawcountries(1.8)
#map.drawmapscale(lon=-90,lat=27,lon0=-98,lat0=40,length=750,barstyle='fancy')
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
#x,y = map(np.array(lon_keep)[w3], np.array(lat_keep)[w3])
#cs2 = map.scatter(x,y,marker='D',s=10,c=NIC[w3],vmin=-20,vmax=20,cmap=cm,zorder=2)
x,y = map(np.array(lon)[w3], np.array(lat)[w3])
cs2 = map.scatter(x,y,marker='v',s=25,c=NIC[w3],vmin=-50,vmax=50,cmap=cm,zorder=2)
x,y = map(np.array(lon)[w2], np.array(lat)[w2])
cs2 = map.scatter(x,y,marker='o',s=45,c=NIC[w2],vmin=-50,vmax=50,cmap=cm,zorder=2)
x,y = map(np.array(lon)[w1], np.array(lat)[w1])
cs2 = map.scatter(x,y,marker='o',s=45,c=NIC[w1],vmin=-50,vmax=50,cmap=cm,zorder=2)
cbar = map.colorbar(cs2,location='bottom',pad="10%")
cbar.set_label("NIC R")
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
plt.rcParams.update({'font.size': 12})
#plt.savefig('images/SMN_NIC_NoahMP_monde.png',format='png',dpi=300)
plt.show()


