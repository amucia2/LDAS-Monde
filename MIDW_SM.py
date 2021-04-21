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

# LOAD EXTRA FILES
#nc=Dataset('mapping.nc', 'r')
#G2P=nc.variables['G2P'][:]
#nc.close()
nc=Dataset('/cnrm/vegeo/bonanb/LDAS/LDAS_curr/Copernicus_2019/MIDW_LDAS/pgd_prep/PREP.nc','r')
#nc=Dataset('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/sfx-trip/pgd_prep/PREP.nc','r')
#nc=Dataset('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/FC/sfx-trip/pgd_prep/PREP.nc','r')

lat_1D=nc.variables['REG_LAT'][:]
lon_1D=nc.variables['REG_LON'][:]
nc.close()

#file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ol/WG3*'),key=numericalSort)
#df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])] 
#wg3_1D = pd.concat(df)

#file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ekf/WG3*'),key=numericalSort)
#df = [pd.read_pickle(file[0]),pd.read_pickle(file[1])]
#wg3_ekf_1D = pd.concat(df)



#wg3_1D=pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ol/WG2_2018-01-01_2018-12-31.PData')
#wg3_ekf_1D=pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/results/ol/WG2_2018-01-01_2018-12-31.PData')
wg3_1D=pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens34/results/ol/WG3_2018-01-01_2018-12-31.PData')
wg3_ekf_1D=pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/ens35/results/ol/WG3_2018-01-01_2018-12-31.PData')

'''
lat_2D=np.full(G2P.shape, np.nan, dtype=np.float)
lat_2D[G2P>0]=lat_1D
lon_2D=np.full(G2P.shape, np.nan, dtype=np.float)
lon_2D[G2P>0]=lon_1D

wg3_2D=np.full((wg3_1D.shape[0],G2P.shape[0],G2P.shape[1]), np.nan, dtype=np.float)
wg3_2D[:,G2P>0]=wg3_1D

wg3_ekf_2D=np.full((wg3_ekf_1D.shape[0],G2P.shape[0],G2P.shape[1]), np.nan, dtype=np.float)
wg3_ekf_2D[:,G2P>0]=wg3_ekf_1D
'''

lat_2D = lat_1D
lon_2D = lon_1D
wg3_2D = wg3_1D
wg3_ekf_2D = wg3_ekf_1D

rng3 = pd.date_range('01/01/2018',periods=wg3_1D.shape[0],freq='D')
#rng3 = pd.date_range('01/01/2017',periods=wg3_1D.shape[0],freq='D')



lat_keep = list()
lon_keep = list()

lat = np.loadtxt('/cnrm/vegeo/albergelc/Project/SMOS/liste_2018.txt')[:,0]
lon = np.loadtxt('/cnrm/vegeo/albergelc/Project/SMOS/liste_2018.txt')[:,1]

R_mod = list() ; bias_mod = list() ; ubrmse_mod = list()
R_ekf = list() ; bias_ekf = list() ; ubrmse_ekf = list()

for coord in range(len(lat)):
    print coord
    if (lat_1D[-1,0] >= lat[coord] >= lat_1D[0,0]) and (lon_1D[0,-1] >= lon[coord] >= lon_1D[0,0]):
        obs_05 = np.loadtxt('/cnrm/vegeo/albergelc/ldas_chain_python/LDAS_Copernicus_2019/USCRN/2018/USCRN_0.05_{0:.4f}_{1:.4f}_obs_2018.DAT'.format(lat[coord],lon[coord]))[:,1]
        obs_05ts = np.loadtxt('/cnrm/vegeo/albergelc/ldas_chain_python/LDAS_Copernicus_2019/USCRN/2018/USCRN_0.05_{0:.4f}_{1:.4f}_obs_2018.DAT'.format(lat[coord],lon[coord]))[:,2]

        #obs_05 = np.loadtxt('/cnrm/vegeo/muciaa/NO_SAVE/US_FC15/USCRN/USCRN_0.05_{0:.4f}_{1:.4f}.DAT'.format(lat[coord],lon[coord]))[:,1]
        #obs_05ts = np.loadtxt('/cnrm/vegeo/muciaa/NO_SAVE/US_FC15/USCRN/USCRN_0.05_{0:.4f}_{1:.4f}.DAT'.format(lat[coord],lon[coord]))[:,2]
        rng = pd.date_range('2018-01-01 00:00:00', periods=len(obs_05), freq='H')


        #rng = pd.date_range('01/01/2018',periods=8760.,freq='H')
        #rng = pd.date_range('01/01/2017',periods=17520,freq='H')

        obs_05bis = pd.Series(obs_05, index = rng)
        obs_05tsbis = pd.Series(obs_05ts, index = rng)
        obs_05bis[obs_05bis<0.] = np.nan
        obs_05bis[obs_05tsbis<4.] = np.nan

        if (len(obs_05bis.resample('D').mean()[~np.isnan(obs_05bis.resample('D').mean())]) > 100):

            #w = np.where((lat_1D<lat[coord]+5) & (lat_1D>lat[coord]-5) & (lon_1D>lon[coord]-5) & (lon_1D<lon[coord]+5))
            dst = list()
            for y_lat in range(len(lat_1D[:,0])):
                for y_lon in range(len(lon_1D[0,:])):
                    dst.append(haversine(lon[coord],lat[coord],lon_1D[0,:][y_lon],lat_1D[:,0][y_lat],y_lon,y_lat))
            dst = np.array(dst)

            #print dst[[np.argsort(dst[:,0])[0]],0],\
            #        lat_1D[0,[dst[[np.argsort(dst[:,0])[0]],3][0]]],lon_1D[0,[dst[[np.argsort(dst[:,0])[0]],4][0]]]

            print dst[[np.argsort(dst[:,0])[0]],0],\
                    lat_1D[[dst[[np.argsort(dst[:,0])[0]],3][0]],0],lon_1D[0,[dst[[np.argsort(dst[:,0])[0]],4][0]]]

            #wlat = np.where(lat_2D == lat_1D[dst[[np.argsort(dst[:,0])[0]],3][0]])
            #wlon = np.where(lon_2D == lon_1D[dst[[np.argsort(dst[:,0])[0]],4][0]])

            #wg3_bis = wg3_2D[dst[[np.argsort(dst[:,0])[0]],3][0]*dst[[np.argsort(dst[:,0])[0]],4][0]]
            #wg3_ekf_bis= wg3_ekf_2D[dst[[np.argsort(dst[:,0])[0]],3][0]*dst[[np.argsort(dst[:,0])[0]],4][0]]


            SSM = list() ; SSM2 = list()
            ### Was reshape(52,80) for MIDW region, for US_FC 175,350
            for i in range(wg3_2D.shape[0]) : 
                SSM.append(wg3_2D.values[i].reshape(52,80)[dst[[np.argsort(dst[:,0])[0]],3][0], \
                        dst[[np.argsort(dst[:,0])[0]],4][0]])
                SSM2.append(wg3_ekf_2D.values[i].reshape(52,80)[dst[[np.argsort(dst[:,0])[0]],3][0], \
                        dst[[np.argsort(dst[:,0])[0]],4][0]])
            SSM_bis = pd.Series(np.array(SSM),index=rng3)
            SSM2_bis = pd.Series(np.array(SSM2),index=rng3)
            wg3_bis = SSM_bis
            wg3_ekf_bis = SSM2_bis

            if len(wg3_bis[~np.isnan(wg3_bis)])>100.:
                '''
                plt.subplot(312) ; plt.plot(obs_05bis.resample('D').mean(),label='in situ',color='k') ;
                plt.plot(wg3_bis,linewidth=1,color='b',label='LDAS-Monde ol')
                plt.plot(wg3_ekf_bis,linewidth=1,color='r',label='LDAS-Monde ekf')
                plt.xticks(rotation=45)
                plt.ylim(0,0.6)
                plt.show()
                #plt.legend(loc='best')
                '''
                obs_05bis_DD = obs_05bis.resample('D').mean()
                ### Anomaly Calculations
                #obs_05bis_DD = \
                #        (obs_05bis_DD-pd.rolling_mean(obs_05bis_DD,30,center=True,min_periods=5))/pd.rolling_std(obs_05bis_DD,30,center=True,min_periods=5)
                #wg3_bis = \
                #        (wg3_bis-pd.rolling_mean(wg3_bis,30,center=True,min_periods=5))/pd.rolling_std(wg3_bis,30,center=True,min_periods=5)
                #wg3_ekf_bis = \
                #        (wg3_ekf_bis-pd.rolling_mean(wg3_ekf_bis,30,center=True,min_periods=5))/pd.rolling_std(wg3_ekf_bis,30,center=True,min_periods=5)
                '''
                plt.subplot(313) ; plt.plot(obs_05bis_DD,label='in situ',color='k') 
                plt.plot(wg3_bis,linewidth=1,color='b',label='LDAS-Monde ol')
                plt.plot(wg3_ekf_bis,linewidth=1,color='r',label='LDAS-Monde ekf')
                plt.xticks(rotation=45)
                plt.show()
                '''
                lat_keep.append(lat[coord]) ; lon_keep.append(lon[coord])
                '''
                print "R ol : %.4f Bias ol : %.4f ubRMSE ol : %.4f"\
                    %(correlation((wg3_bis).values,\
                obs_05bis.resample('D').mean().values)[0],\
                bias((wg3_bis).values,\
                obs_05bis.resample('D').mean().values),\
                ubrmse((wg3_bis).values,\
                obs_05bis.resample('D').mean().values))
                print "R ekf : %.4f Bias ekf : %.4f ubRMSE ekf : %.4f"\
                    %(correlation((wg3_ekf_bis).values,\
                obs_05bis.resample('D').mean().values)[0],\
                bias((wg3_ekf_bis).values,\
                obs_05bis.resample('D').mean().values),\
                ubrmse((wg3_ekf_bis).values,\
                obs_05bis.resample('D').mean().values))
                '''
                print "R ol : %.4f Bias ol : %.4f ubRMSE ol : %.4f"\
                    %(correlation((wg3_bis).values,\
                obs_05bis_DD.values)[0],\
                bias((wg3_bis).values,\
                obs_05bis_DD.values),\
                ubrmse((wg3_bis).values,\
                obs_05bis_DD.values))
                
                print "R ekf : %.4f Bias ekf : %.4f ubRMSE ekf : %.4f"\
                    %(correlation((wg3_ekf_bis).values,\
                obs_05bis_DD.values)[0],\
                bias((wg3_ekf_bis).values,\
                obs_05bis_DD.values),\
                ubrmse((wg3_ekf_bis).values,\
                obs_05bis_DD.values))
                
                '''
                R_mod.append(correlation((wg3_bis).values,obs_05bis.resample('D').mean().values)[0])
                bias_mod.append(bias((wg3_bis).values,obs_05bis.resample('D').mean().values))
                ubrmse_mod.append(ubrmse((wg3_bis).values,obs_05bis.resample('D').mean().values))
                R_ekf.append(correlation((wg3_ekf_bis).values,obs_05bis.resample('D').mean().values)[0])
                bias_ekf.append(bias((wg3_ekf_bis).values,obs_05bis.resample('D').mean().values))
                ubrmse_ekf.append(ubrmse((wg3_ekf_bis).values,obs_05bis.resample('D').mean().values))
                '''
                R_mod.append(correlation((wg3_bis).values,obs_05bis_DD.values)[0])
                bias_mod.append(bias((wg3_bis).values,obs_05bis_DD.values))
                ubrmse_mod.append(ubrmse((wg3_bis).values,obs_05bis_DD.values))
                R_ekf.append(correlation((wg3_ekf_bis).values,obs_05bis_DD.values)[0])
                bias_ekf.append(bias((wg3_ekf_bis).values,obs_05bis_DD.values))
                ubrmse_ekf.append(ubrmse((wg3_ekf_bis).values,obs_05bis_DD.values))

### Save Correlations into files for later use
R_mod = pd.DataFrame(R_mod)
R_ekf = pd.DataFrame(R_ekf)
R_mod.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens34.PData')
R_ekf.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/ens_R/ens35.PData')


sys.exit()



plt.figure(1)
plt.subplot(3,3,1)
plt.scatter(R_ekf,R_mod,color='r',s=8)
plt.rcParams.update({'font.size':10})
plt.plot( [0,1],[0,1] )
plt.xlabel('LDAS-Monde EKF vs USCRN')
plt.ylabel('LDAS-Monde OL vs USCRN')# / e5_025-green')
plt.xlim(0.4,1)
plt.ylim(0.4,1)
plt.title('R : WG3')
plt.subplot(3,3,2)
plt.title('ubRMSE : WG3')
plt.scatter(ubrmse_ekf,ubrmse_mod,color='r',s=8)
plt.plot( [0,0.15],[0,0.15])
plt.xlim(0,0.15)
plt.ylim(0,0.15)
plt.subplot(3,3,3)
plt.title('abs(Bias) : LE')
plt.scatter(np.abs(bias_ekf),np.abs(bias_mod),color='r',s=8)
plt.xlim(0,0.30)
plt.ylim(0,0.30)
plt.plot( [0,0.30],[0,0.30] )
plt.show()


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


NIC = 100*(np.array(R_ekf)-np.array(R_mod))/(1-np.array(R_mod))

w1=np.where(NIC>3)
w3=np.where((NIC>-3) & (NIC<3))
w2=np.where(NIC<-3)

fig, axes = plt.subplots(1,1)
axes.set_title("Evaluation vs USCRN")
map = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-55, urcrnrlon=180, urcrnrlat=90,resolution='c',ax = axes)
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines()
map.drawmapboundary(fill_color='lightblue')
map.drawcountries()
map.fillcontinents(color='beige',lake_color='lightblue')
#x,y = map(np.array(lon_keep), np.array(lat_keep))
#cs2 = map.scatter(x,y,marker='s',s=30,c=(np.array(R_ekf)-np.array(R_mod)),vmin=-0.1,vmax=0.1,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w3], np.array(lat_keep)[w3])
cs2 = map.scatter(x,y,marker='D',s=10,c=NIC[w3],vmin=-20,vmax=20,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w2], np.array(lat_keep)[w2])
cs2 = map.scatter(x,y,marker='o',s=30,c=NIC[w2],vmin=-20,vmax=20,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w1], np.array(lat_keep)[w1])
cs2 = map.scatter(x,y,marker='o',s=30,c=NIC[w1],vmin=-20,vmax=20,cmap=cm,zorder=2)
cbar = map.colorbar(cs2,location='bottom',pad="5%")
plt.show()

