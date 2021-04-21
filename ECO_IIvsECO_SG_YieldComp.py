import sys, glob, os, re, math, time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
import pylab as pl
import scipy.optimize as optim
import scipy as sp
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import zscore
from scipy import stats

###################################################################
### Useful Functions
###################################################################

def pearsonr_ci(x,y,alpha=0.05):
    ''' calculate Pearson correlation along with the confidence interval using scipy and numpy
    Parameters
    ----------
    x, y : iterable object such as a list or np.array
      Input for correlation calculation
    alpha : float
      Significance level. 0.05 by default
    Returns
    -------
    r : float
      Pearson's correlation coefficient
    pval : float
      The corresponding p value
    lo, hi : float
      The lower and upper bound of confidence intervals
    '''
    r, p = stats.pearsonr(x,y)
    r_z = np.arctanh(r)
    se = 1/np.sqrt(x.size-3)
    z = stats.norm.ppf(1-alpha/2)
    lo_z, hi_z = r_z-z*se, r_z+z*se
    lo, hi = np.tanh((lo_z, hi_z))
    return r, p, lo, hi

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

def ZscoreM(x):
    for ii in [1,2,3,4,5,6,7,8,9,10,11,12]:
        for name in x.columns:
            season = (x[name].index.month == ii)
            x[name].ix[season] = (x[name].ix[season]-x[name].ix[season].mean())/x[name].ix[season].std()
    return x

def ZscoreY(x):
    for ii in [1,2,3,4,5,6,7,8,9,10,11,12]:
        season = (x.index.month == ii)
        x.ix[season] = (x.ix[season]-x.ix[season].mean())/x.ix[season].std()
    return x

###################################################################
### Import Data  
###################################################################

### Import individual yearly files
lai_2= pd.read_pickle('Data/LAI_2000_2018_Nebraska.PData')
#lai_c4 = pd.read_pickle('../LAI_C4_2000_2018_Nebraska.PData')
#corn_yn = pd.read_csv('NE_Corn_Production_2000_2018.csv')
corn_y = pd.read_csv('Data/NE_Corn_Production_2000_2018_bu-acre.csv')
'''
f = open('CO2air.txt','r')
f = map(lambda s: s.strip(), f)

t = open('Tair.txt','r')
t = map(lambda s: s.strip(), t)

r = open('Rainf.txt', 'r')
r = map(lambda s: s.strip(), r)

d = open('DIR_SWdown.txt', 'r')
d = map(lambda s: s.strip(), d)

days = pd.date_range(start = '1/1/2000', end = '8/31/2018', freq='M')

CO2 = pd.DataFrame({'Date':days,'CO2':f})
CO2 = CO2.set_index('Date')
Tair = pd.DataFrame({'Date':days,'Tair':t})
Tair = Tair.set_index('Date')
Rainf = pd.DataFrame({'Date':days,'Rainf':r})
Rainf = Rainf.set_index('Date')
SWdown = pd.DataFrame({'Date':days, 'DIR_SWdown':d})
SWdown = SWdown.set_index('Date')
'''

### Import entire folders of yearly data
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    return parts
'''
### It is required to have as many 'read_pickle's as years in the folder 
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_LDAS/observations/sfx-trip/LAI_MC.t08*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15])]
lai_c4_obs = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_SGC/results/ol/LAI_2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
     pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
lai_ol_cco2 = pd.concat(df)
'''
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_ECOSG/results/ol/LAI_2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
lai_sg_ol = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_ECOSG/results/ekf/LAI_2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
lai_sg_ekf = pd.concat(df)


###################################################################
### Data Processing 
###################################################################
#lai_c4_obs_1 = lai_c4_obs.mean(axis=1)
#lai_ol_cco2 = lai_ol_cco2.mean(axis=1)

lai_xM = lai_2.resample('M').mean()
lai_xY = lai_2.resample('A').mean()

#lai_c4_xM = lai_c4.resample('M',label='left').mean()
#lai_c4_xY = lai_c4.resample('A').mean()
#lai_c4_obs_xY = lai_c4_obs_1.resample('A').mean()

#lai_ol_cco2_xY = lai_ol_cco2.resample('A',label='left').mean() 

#Tair = Tair.apply(pd.to_numeric)
#CO2 = CO2.apply(pd.to_numeric)
#Rainf = Rainf.apply(pd.to_numeric)
#SWdown = SWdown.apply(pd.to_numeric)

#CO2Y = CO2.resample('A',label='left').mean()
#TairY = Tair.resample('A',label='left').mean()
#RainfY = Rainf.resample('A',label='left').mean()
#SWdownY = SWdown.resample('A', label='left').mean()

lai_sg_olY = lai_sg_ol.resample('A',label='left').mean()
lai_sg_ekfY = lai_sg_ekf.resample('A',label='left').mean()

#sys.exit()
###################################################################
### Computer Z Scores
###################################################################
corn_yz = corn_y.apply(zscore)
#lai_yz = lai_2.apply(zscore)
#corn_ynz = corn_yn.apply(zscore)

#tairyz = TairY.apply(zscore)
#tairz = Tair.apply(zscore)
#co2yz = CO2Y.apply(zscore)
#co2z = CO2.apply(zscore)
#rainfyz = RainfY.apply(zscore)
#rainfz = Rainf.apply(zscore)
#swdownyz = SWdownY.apply(zscore)
#swdownz = SWdown.apply(zscore)

lai_M = ZscoreM(lai_xM)
lai_Y = ZscoreY(lai_xY)
#lai_c4_M = ZscoreM(lai_c4_xM)
#lai_c4_Y = ZscoreY(lai_c4_xY)
#lai_c4_obs_Y = ZscoreY(lai_c4_obs_xY)
#lai_ol_cco2_Y = ZscoreY(lai_ol_cco2_xY)

#lai_Y['Yield'] = corn_ynz['Value'].values
#lai_c4_Y['Yield'] = corn_ynz['Value'].values
lai_Y['Yield_b/a'] = corn_yz['Value'].values
#lai_c4_Y['Yield_b/a'] = corn_yz['Value'].values

lai_sg_olz = ZscoreY(lai_sg_olY)
lai_sg_ekfz = ZscoreY(lai_sg_ekfY)

#sys.exit()

#tt=list()
#for i in range(len(lai_c4_Y['Model'])):
#    tt.append(np.nan)

#tt[0:16] = lai_c4_obs_Y.values
#lai_c4_Y['Obs'] = tt

#corn_yz = corn_yz['Value'].values
#corn_ynz = corn_ynz['Value'].values

#lai_c4_Y.index = lai_c4_Y.index+pd.DateOffset(month=1,day=1)
lai_Y.index = lai_Y.index+pd.DateOffset(month=1,day=1)

###################################################################
### Computer Correlations
###################################################################

#cor_obsc4_yieldn = np.corrcoef(lai_c4_obs_Y.values,corn_ynz[0:16])
#cor_obs_yieldn = np.corrcoef(lai_Y['Obs'].values,corn_ynz) 
#cor_mod_yieldn = np.corrcoef(lai_Y['Model'].values,corn_ynz)
#cor_anl_yieldn = np.corrcoef(lai_Y['Analysis'].values,corn_ynz)

#cor_obsc4_yield = np.corrcoef(lai_c4_obs_Y.values,corn_yz[0:16])
cor_obs_yield = np.corrcoef(lai_Y['Obs'].values,lai_Y['Yield_b/a'].values)
cor_mod_yield = np.corrcoef(lai_Y['Model'].values,lai_Y['Yield_b/a'].values)
cor_anl_yield = np.corrcoef(lai_Y['Analysis'].values,lai_Y['Yield_b/a'].values)
#cor_c4mod_yield = np.corrcoef(lai_c4_Y['Model'].values,corn_yz)
#cor_c4anl_yield = np.corrcoef(lai_c4_Y['Analysis'].values,corn_yz)



cor_mod_ob = np.corrcoef(lai_Y['Model'].values,lai_Y['Obs'].values)
cor_anl_ob = np.corrcoef(lai_Y['Analysis'].values,lai_Y['Obs'].values)
cor_sg_yield = np.corrcoef(lai_sg_olz.mean(axis=1).values,lai_Y['Yield_b/a'].values)
cor_sg_ob = np.corrcoef(lai_sg_olz.mean(axis=1).values,lai_Y['Obs'].values)


cor_sg_ekf_yield = np.corrcoef(lai_sg_ekfz.mean(axis=1).values,lai_Y['Yield_b/a'].values)
cor_sg_ekf_ob = np.corrcoef(lai_sg_ekfz.mean(axis=1).values,lai_Y['Obs'].values)

correlations = pd.DataFrame(columns=['Exp','Cor','P-val','Lower bound','Upper bound'])

#c4Y_c4Obs = pearsonr_ci(lai_c4_obs_Y.values,corn_yz[0:16])
#c4Y_c4Mod = pearsonr_ci(lai_c4_Y['Model'].values,corn_yz)
#c4Y_c4Anl = pearsonr_ci(lai_c4_Y['Analysis'].values,corn_yz)
#c4Y_Obs = pearsonr_ci(lai_Y['Obs'].values,corn_yz)
#c4Y_Mod = pearsonr_ci(lai_Y['Model'].values,corn_yz)
#c4Y_Anl = pearsonr_ci(lai_Y['Analysis'].values,corn_yz)



#print("Yearly C4 LAI Observed Anomaly vs Corn Yield Anomaly: " + str(round(cor_obsc4_yield[0,1],4)))
#print("Yearly C4 LAI Open-Loop Anomaly vs Corn Yield Anomaly: " + str(round(cor_c4mod_yield[0,1],4)))
#print("Yearly C4 LAI EKF Anomaly vs Corn Yield Anomaly: " + str(round(cor_c4anl_yield[0,1],4)))
print("Yearly LAI Observed Anomaly vs Corn Yield Anomaly: " + str(round(cor_obs_yield[0,1],4)))


print("Yearly LAI Open-Loop (ECO_II) Anomaly vs Corn Yield Anomaly: " + str(round(cor_mod_yield[0,1],4)))
print("Yearly LAI Open-Loop (ECO_SG) Anomaly vs Corn Yield Anomaly: " + str(round(cor_sg_yield[0,1],4)))
print("Yearly LAI Open-Loop (ECO_II) Anomaly vs Observed LAI Anomaly: " +str(round(cor_mod_ob[0,1],4)))
print("Yearly LAI Open-Loop (ECO_SG) Anomaly vs Observed LAI Anomaly: " +str(round(cor_sg_ob[0,1],4)))


print("Yearly LAI EKF (ECO_II) Anomaly vs Corn Yield Anomaly: " + str(round(cor_anl_yield[0,1],4)))
print("Yearly LAI EKF (ECO_SG) Anomaly vs Corn Yield Anomaly: " + str(round(cor_sg_ekf_yield[0,1],4)))
print("Yearly LAI EKF (ECO_II) Anomaly vs Observed LAI Anomaly: " + str(round(cor_anl_ob[0,1],4)))
print("Yearly LAI EKF (ECO_SG) Anomaly vs Observed LAI Anomaly: " + str(round(cor_sg_ekf_ob[0,1],4)))


#sys.exit()
###################################################################
### Graphing Functions 
###################################################################

raw_input("Press Enter to continue ...")
### plt.subplots(num_rows,num_cols)

#plt.figure(figsize=(8,6),dpi=600,edgecolor='w')
fig, axs = plt.subplots(2,1)

axs[0].set_title('LAI and Corn Yield Interannual Variability over Nebraska')
#axs[0].plot(co2yz['CO2'],label='CO2',color='black',linewidth =1)
axs[0].plot(lai_Y['Yield_b/a'],label='Corn Yield - USDA',marker='o',linestyle='--',markersize=9,color='green')

axs[0].plot(lai_Y['Obs'],label='LAI Observations',color='black',marker='*',linestyle='--',markersize=9)
#axs[0].plot(lai_sg_olz.mean(axis=1),label='OL LAI',color='blue')
#axs[0].plot(lai_sg_ekfz.mean(axis=1),label='EKF LAI',color='red')
axs[0].plot(lai_Y['Analysis'],label='EKF LAI',color='red')
axs[0].plot(lai_Y['Model'],label='Open-Loop LAI',color='blue')


axs[0].legend(loc='upper left',fontsize=14)
axs[0].set_ylim([-2.5,2.5])
axs[0].axhline(color='black',linewidth=.5)
axs[0].set_ylabel('Standard Anomaly')
plt.rcParams.update({'font.size': 18})
fig.tight_layout()
fig.patch.set_facecolor('white')
#plt.savefig('images/NE_LAI_2000-2018.png',format='png',dpi=600)
plt.show()



#axs[1].set_title('LAI and Corn Yield Interannual Variability over Nebraska')
##axs[1].plot(co2yz['CO2'],label='CO2',color='black',linewidth =1)
#axs[1].plot(lai_Y['Yield_b/a'],label='Corn Yield - USDA',marker='o',linestyle='--',markersize=9,color='green')
#axs[1].plot(lai_Y['Obs'],label ='LAI Observations',marker='*',linewidth=0,markersize=9)
#axs[1].plot(lai_Y['Analysis'],label='EKF LAI',color='red',linewidth=2)
#axs[1].plot(lai_Y['Model'],label='Open-Loop LAI',color='blue',linewidth=2)
#axs[1].legend(loc='upper left')
#axs[1].set_ylim([-2.5,2.5])
#axs[1].axhline(color='black',linewidth=.5)
#axs[1].set_ylabel('Standard Anomaly')

#x= np.arange(2000,2019,1)

#y=tairyz['Tair'].values 
#ly=lai_Y['Model'].values
#oy=lai_Y['Yield_b/a'].values

#z = np.polyfit(x,y,1)
#p = np.poly1d(z)

#lz = np.polyfit(x,ly,1)
#lp = np.poly1d(lz)

#oz = np.polyfit(x,oy,1)
#op = np.poly1d(oz)

#lai_Y['t_trend'] = p(x)
#lai_Y['l_trend'] = lp(x)
#lai_Y['o_trend'] = op(x)
#diff = lai_Y['Model'].values-lai_ol_cco2_Y.values
'''
#plt.title('LAI and Corn Yield Interannual Variability over Nebraska - Variable and Const CO2')
#plt.title('Model LAI, Observed Corn Yield, and Rainfall Anomalies over Nebraska')

plt.title('LAI from ECO_II vs ECO_SG with Observed Corn Yield over Nebraska')
plt.plot(lai_Y['Yield_b/a'],label='Corn Yield - USDA',marker='o',linestyle='--',markersize=9,color='green')

plt.plot(lai_Y['Obs'],label='LAI Observations',color='black',marker='*',linestyle='--',markersize=9)
#plt.plot(co2yz['CO2'],label='CO2',color='black',linewidth =1)
#plt.plot(tairyz['Tair'],label='Tair',color='orange')
#plt.plot(rainfyz['Rainf'],label='Rainfall',color='red')
#plt.plot(swdownyz['DIR_SWdown'],label='SW Down',color='orange')
#plt.plot(lai_c4_Y['Obs'],label ='C4 LAI Observations',marker='*',linewidth=0,markersize=9)
#plt.plot(lai_c4_Y['Analysis'],label='C4 EKF LAI',color='red')
#plt.plot(lai_Y['Model'],label='Open-Loop LAI',color='blue')
#plt.plot(lai_ol_cco2_Y,label='Open-Loop LAI - Cons. CO2',color='red')
#plt.plot(diff,label='Open-Loop LAI - Cons. CO2',color='red')

#plt.plot(lai_Y['Model'],label='OL ECO_II LAI',color='blue')
#plt.plot(lai_Y['Analysis'],label='EKF ECO_II LAI',color='red')
plt.plot(lai_sg_olz.mean(axis=1),label='OL ECO_SG LAI',color='green')
plt.plot(lai_sg_ekfz.mean(axis=1),label='EKF ECO_SG LAI',color='orange')

plt.ylim([-2.5,2.5])

#plt.plot(lai_Y['t_trend'],label='Temperature Trend',linestyle="--",color='red')
#plt.plot(lai_Y['l_trend'],label='LAI Trend',linestyle="--",color='blue')
#plt.plot(lai_Y['o_trend'],label='Corn Yield Trend',linestyle="--",color='green')

plt.axhline(color='black',linewidth=.5)
plt.legend(loc='upper left')

plt.rcParams.update({'font.size': 10})
#fig.tight_layout()
#save('../vegeoFigures/Yield-RainF', ext='ps', close=True, verbose=True)
plt.show()
'''
