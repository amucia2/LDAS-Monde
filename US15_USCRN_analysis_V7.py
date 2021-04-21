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
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

# LOAD EXTRA FILES
### MIDW Region
#nc=Dataset('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/MIDW/FC/sfx-trip/pgd_prep/PREP.nc','r')

### CONUS Region
nc=Dataset('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15_SG/sfx-trip/pgd_prep/PREP.nc','r')

lat_1D=nc.variables['REG_LAT'][:]
lon_1D=nc.variables['REG_LON'][:]
nc.close()

lat_keep = list()
lon_keep = list()

lat = np.loadtxt('/cnrm/vegeo/albergelc/Project/SMOS/liste_2018.txt')[:,0]
lon = np.loadtxt('/cnrm/vegeo/albergelc/Project/SMOS/liste_2018.txt')[:,1]

layer = ('WG3','WG6','WG8')
ana = ('R')

### Negative List - that is which stations to exclude
lwg3 = [29]
lwg6 = [1,3,12,19,26,28,29,32,43,44,46,65,69,83,88,96,97,106,107]
lwg8 = [1,3,19,22,26,29,29,32,43,44,46,65,69,83,88,96,97,106,108] 
### Postive List - which stations to include for analysis of all three
lall = [6,7,9,11,13,23,25,27,30,31,34,35,36,37,38,39,40,41,47,48,49,50,51,53,54,55,56,57,59,61,66,67,68,70,71,72,74,77,78,87,89,91,92,93,98,99,101,103,104,107,109]

gldas = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R/gldas/GLDAS_0_WG3_0.05_R.PData')
gldas = gldas[gldas.index.isin(lall)]
gldas = [gldas.mean(),np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
nca = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R/nca/NoahMP_LAI_0_WG3_0.05_R.PData')
nca = nca[nca.index.isin(lall)]
nca = [nca.mean(),np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
era5 = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R/era5/ERA5_0_WG3_0.05_R.PData')
era5 = era5[era5.index.isin(lall)]
era5 = [era5.mean(),np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]

us6ol = [] ; us6ekf = [] ; us6ol8 = [] ; us6ekf8 = [] ; us6ol20 = [] ; us6ekf20 = [] ; usifsol = [] ; usifsol8 = []; usifsol2 = []
us15ol = [] ; us15ekf = [] ; us15ol8 = [] ; us15ekf8 = [] ; us15ol20 = [] ; us15ekf20 = [] ; usifsekf = [];usifsekf8 = [];usifsekf2 = []

### Import and Assign US_FC_IFS OL and EKF FC at diff Depths
usifsol.append(np.nan)
usifsekf.append(np.nan)
usifsol8.append(np.nan)
usifsekf8.append(np.nan)
usifsol2.append(np.nan)
usifsekf2.append(np.nan)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US_FC_IFS/ol*WG3*0.05*R*.PData'),key=natural_keys)
for k in range(len(file)):
    usifsol.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US_FC_IFS/ekf*WG3*0.05*R*.PData'),key=natural_keys)
for k in range(len(file)):
    usifsekf.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US_FC_IFS/ol*WG8*1.00*R*.PData'),key=natural_keys)
for k in range(len(file)):
    usifsol8.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US_FC_IFS/ekf*WG8*1.00*R*.PData'),key=natural_keys)
for k in range(len(file)):
    usifsekf8.append(pd.read_pickle(file[k]).mean())
    
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US_FC_IFS/ol*WG_20*0.20*R*.PData'),key=natural_keys)
for k in range(len(file)):
    usifsol2.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US_FC_IFS/ekf*WG_20*0.20*R*.PData'),key=natural_keys)
for k in range(len(file)):
    usifsekf2.append(pd.read_pickle(file[k]).mean())
    
### Import and Assign LDAS-Monde OL and EKF FC at diff Depths
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US6/ol*WG3*0.05*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us6ol.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US6/ekf*WG3*0.05*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us6ekf.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US6/ol*WG8*1.00*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us6ol8.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US6/ekf*WG8*1.00*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us6ekf8.append(pd.read_pickle(file[k]).mean())
    
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US6/ol*WG_20*0.20*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us6ol20.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US6/ekf*WG_20*0.20*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us6ekf20.append(pd.read_pickle(file[k]).mean())
###-------------------------------------------------------------------###    
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US15/ol*WG3*0.05*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us15ol.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US15/ekf*WG3*0.05*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us15ekf.append(pd.read_pickle(file[k]).mean())
    
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US15/ol*WG8*1.00*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us15ol8.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US15/ekf*WG8*1.00*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us15ekf8.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US15/ol*WG_20*0.20*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us15ol20.append(pd.read_pickle(file[k]).mean())

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US15/ekf*WG_20*0.20*R*.PData'),key=natural_keys)
for k in range(len(file)):
    us15ekf20.append(pd.read_pickle(file[k]).mean())
### ------------------PLOTTING---------------###


ll_5_ol = 0.027852, 0.027315, 0.028864, 0.030044, 0.0303129999999999, 0.031033, 0.033068, 0.034512, 0.035497, 0.036586, 0.037767, 0.039919, 0.040266, 0.039803 
ul_5_ol = 0.026944, 0.027208, 0.028818, 0.029709, 0.030926, 0.031538, 0.034183, 0.034679, 0.038021, 0.0367959999999999, 0.0402279999999999, 0.041157, 0.04137, 0.0414370000000001
ll_5_ekf =  0.0291710000000001, 0.0286150000000001, 0.029616, 0.030841, 0.031479, 0.031325, 0.034388, 0.0344749999999999, 0.0356620000000001, 0.036598, 0.0400700000000001, 0.0383039999999999, 0.040551, 0.040086
ul_5_ekf =  0.0280929999999999, 0.0268769999999999, 0.028387, 0.028269, 0.029847, 0.0326879999999999, 0.034265, 0.034408, 0.035496, 0.0368999999999999, 0.039362, 0.0407770000000001, 0.04236, 0.04295

ll_20_ol = 0.0440739999999999, 0.043661, 0.0443939999999999, 0.042863, 0.0447420000000001, 0.0426150000000001, 0.0444880000000001, 0.0424209999999999, 0.043498, 0.045072, 0.047054, 0.046868, 0.048196, 0.048075
ul_20_ol = 0.0424410000000001, 0.0421940000000001, 0.042142, 0.040633, 0.0418649999999999, 0.041655, 0.043195, 0.042863, 0.04453, 0.047011, 0.045716, 0.0445789999999999, 0.0507190000000001, 0.0493490000000001
ll_20_ekf = 0.049458, 0.0488729999999999, 0.049667, 0.048132, 0.0473479999999999, 0.0471, 0.047197, 0.047168, 0.04687, 0.047675, 0.0473479999999999, 0.05096, 0.049523, 0.05029
ul_20_ekf = 0.0449999999999999, 0.043328, 0.0442509999999999, 0.043385, 0.0448040000000001, 0.045913, 0.044975, 0.045234, 0.045231, 0.047647, 0.048078, 0.047099, 0.049648, 0.050917

us15ekf = us15ekf[1:]
us15ol = us15ol[1:]

us15ekf20 = us15ekf20[1:]
us15ol20 = us15ol20[1:]

ana = input("Combined (COM) or Individual (IND)? : ")

if ana == 'IND':
    fig, axes = plt.subplots(2,1)
    x = 0,1,2,3,4,5,6,7,8,9,10,11,12,13
    ticks = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14']
    axes[0].set_title('LDAS-Monde Forecast Soil Moisture vs In-Situ Observations (USCRN) : 5cm',fontsize=7)
    axes[1].set_title('LDAS-Monde Forecast Soil Moisture vs In-Situ Observations (USCRN) : 20cm',fontsize=7)
    axes[0].set_xticks(np.arange(14))
    axes[1].set_xticks(np.arange(14))
    axes[0].set_xticklabels(ticks)
    axes[1].set_xticklabels(ticks)
    axes[0].set_ylabel('Correlation')
    axes[1].set_ylabel('Correlation')
    #plt.ylabel('Correlation (R)')
    axes[1].set_xlabel('fc Lead Time (Days)')
    axes[0].set_xlim(-0.5,13.5)
    axes[1].set_xlim(-0.5,13.5)

    import matplotlib.transforms as transforms
    offset = transforms.ScaledTranslation(0.05, 0,fig.dpi_scale_trans)

    trans = axes[0].transData-offset
    trans2 = axes[0].transData+offset

    trans_ = axes[1].transData-offset
    trans2_ = axes[1].transData+offset

    axes[0].plot(us15ol,color='b',marker='h',linewidth=1,label='fc_init_ol',transform=trans,markeredgecolor='black')[0]
    axes[0].plot(us15ekf,color='r',marker='h',linewidth=1,label='fc_init_ekf',transform=trans2,markeredgecolor='black')[0]

    axes[1].plot(us15ol20,color='b',marker='h',linewidth=1,transform=trans_,markeredgecolor='black')[0]
    axes[1].plot(us15ekf20,color='r',marker='h',linewidth=1,transform=trans2_,markeredgecolor='black')[0]

    ### Error Bars
    '''
    axes[0].errorbar(x,us15ol,yerr=[ll_5_ol,ul_5_ol],transform=trans,color='b')
    axes[0].errorbar(x,us15ekf,yerr=[ll_5_ekf,ul_5_ekf],transform=trans2,color='r')

    axes[1].errorbar(x,us15ol20,yerr=[ll_20_ol,ul_20_ol],transform=trans_,color='b')
    axes[1].errorbar(x,us15ekf20,yerr=[ll_20_ekf,ul_20_ekf],transform=trans2_,color='r')

    axes[0].axvline(x=0.5,color='black',linestyle='--',linewidth=1)
    axes[1].axvline(x=0.5,color='black',linestyle='--',linewidth=1)
    '''
    axes[0].set_ylim(0.45,0.8)
    axes[1].set_ylim(0.45,0.8)

    axes[0].set_yticks(np.arange(0.50,0.80,step=0.1))
    axes[1].set_yticks(np.arange(0.50,0.80,step=0.1))
    axes[0].annotate("A", xy=(0.025, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')
    axes[1].annotate("B", xy=(0.025, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')

elif ana == 'COM':
    fig = plt.figure(figsize=(14,6))
    ax1 = fig.add_subplot(111)
    x = 0,1,2,3,4,5,6,7,8,9,10,11,12,13
    ticks = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14']
    ax1.set_title('LDAS-Monde Forecast Soil Moisture vs In-Situ Observations (USCRN)')
    ax1.set_xlabel('fc Lead Time (Days)')
    ax1.set_xlim(-0.5,13.5)
    ax1.set_ylabel('Correlation')
    ax1.set_ylim(0.45,0.8)
    ax1.set_yticks(np.arange(0.50,0.80,step=0.1))

    import matplotlib.transforms as transforms
    offset = transforms.ScaledTranslation(0.05, 0,fig.dpi_scale_trans)
    trans = ax1.transData-offset
    trans2 = ax1.transData+offset

    ax1.plot(us15ol,color='b',marker='.',linewidth=1,label='fc_init_ol 5cm (106 stations)',transform=trans,markeredgecolor='black',markersize=14)[0]
    ax1.plot(us15ekf,color='r',marker='.',linewidth=1,label='fc_init_sekf 5cm',transform=trans2,markeredgecolor='black',markersize=14)[0]
    ax1.plot(us15ol20,color='b',marker='s',linewidth=1,label='fc_init_ol 20cm (84 stations)',transform=trans,markeredgecolor='black',linestyle='--',markersize=8)[0]
    ax1.plot(us15ekf20,color='r',marker='s',linewidth=1,label='fc_init_sekf 20cm',transform=trans2,markeredgecolor='black',linestyle='--',markersize=8)[0]
    ax1.set_xticks(x)
    ax1.set_xticklabels(ticks)

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles,labels,numpoints=1,fontsize=12,bbox_to_anchor=(0.985,0.92))



#fig.legend(numpoints=1,fontsize=8,loc=(0.8,0.83))
#plt.rcParams.update({'font.size':8})
fig.tight_layout()
plt.show()

sys.exit()
















### NIC Maps ###

root_15 = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US15/'
root_ifs = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US_FC_IFS/'

ens_1 = 'ol'
ens_2 = 'ekf'

fc0 = ('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14')
wg = ('WG3','WG3','WG3','WG3','WG3','WG3','WG3','WG3','WG3','WG3','WG3','WG3','WG3','WG3','WG3')
#wg = ('WG8','WG8','WG8','WG8','WG8','WG8','WG8','WG8','WG8','WG8','WG8','WG8','WG8','WG8','WG8')
#wg = ('WG_20','WG_20','WG_20','WG_20','WG_20','WG_20','WG_20','WG_20','WG_20','WG_20','WG_20','WG_20','WG_20','WG_20','WG_20')
dep = ('0.05','0.05','0.05','0.05','0.05','0.05','0.05','0.05','0.05','0.05','0.05','0.05','0.05','0.05','0.05')
#dep = ('1.00','1.00','1.00','1.00','1.00','1.00','1.00','1.00','1.00','1.00','1.00','1.00','1.00','1.00','1.00')
#dep = ('0.20','0.20','0.20','0.20','0.20','0.20','0.20','0.20','0.20','0.20','0.20','0.20','0.20','0.20','0.20')

lat = np.loadtxt('/cnrm/vegeo/albergelc/Project/SMOS/liste_2018.txt')[:,0]
lon = np.loadtxt('/cnrm/vegeo/albergelc/Project/SMOS/liste_2018.txt')[:,1]
lat = pd.DataFrame(lat)
lon = pd.DataFrame(lon)



for i in range(len(fc0)):
    P_mod = pd.read_pickle(root_15+'ol_{0}_{1}_{2}_P.PData'.format(fc0[i],wg[i],dep[i]))
    R_mod = pd.read_pickle(root_15+'ol_{0}_{1}_{2}_R.PData'.format(fc0[i],wg[i],dep[i]))
    P_ekf = pd.read_pickle(root_15+'ekf_{0}_{1}_{2}_P.PData'.format(fc0[i],wg[i],dep[i]))
    R_ekf = pd.read_pickle(root_15+'ekf_{0}_{1}_{2}_R.PData'.format(fc0[i],wg[i],dep[i]))

    from matplotlib import ticker
    NIC = 100*(np.array(R_ekf)-np.array(R_mod))/(1-np.array(R_mod))
    lat = np.array(lat)[~np.isnan(np.array(NIC))]
    lon = np.array(lon)[~np.isnan(np.array(NIC))]
    NIC = np.array(NIC)[~np.isnan(np.array(NIC))]

    w1=np.where(NIC>3)[0]
    w3=np.where((NIC>-3) & (NIC<3))[0]
    w2=np.where(NIC<-3)[0]

#    print("FC Day : {0}".format(fc0[i]))
    print("Number of Stations with EKF Improvement : {0}".format(len(w1)))
    print("Number of Stations with EKF Degredation : {0}".format(len(w2)))
    print("Number of Stations with Neutral Impact : {0}".format(len(w3)))
    
    plt.figure(figsize=(8,5),dpi=250,edgecolor='w')
    plt.title("IFS NIC Correlation EKF vs OL - 20cm")
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
    x,y = map(lon[w3],lat[w3])
    cs2 = map.scatter(x,y,marker='v',s=20,c=NIC[w3],vmin=-15,vmax=15,cmap=cm,zorder=2,edgecolors='black')
    x,y = map(lon[w2],lat[w2])
    cs2 = map.scatter(x,y,marker='o',s=35,c=NIC[w2],vmin=-15,vmax=15,cmap=cm,zorder=2,edgecolors='black')
    x,y = map(lon[w1], lat[w1])
    cs2 = map.scatter(x,y,marker='o',s=35,c=NIC[w1],vmin=-15,vmax=15,cmap=cm,zorder=2,edgecolors='black')
    cbar = map.colorbar(cs2,location='bottom',pad="10%")
    cbar.set_label("NIC R")
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    plt.rcParams.update({'font.size': 10})
    #plt.tight_layout()
    plt.savefig('NIC_IFS_USCRN_20cm.png',format='png',dpi=300)
    plt.show()

sys.exit()
