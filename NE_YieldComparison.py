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
lai_iir = pd.read_pickle('Data/Nebraska/LAI_2000_2018_Nebraska.PData')
lai_ii_c4r = pd.read_pickle('Data/Nebraska/LAI_C4_2000_2018_Nebraska.PData')
corn_y = pd.read_csv('Data/Nebraska/NE_Corn_Production_2000_2018_bu-acre.csv')
lai_sgr = pd.read_pickle('Data/old/ECOSG_LAI.PData')
swi_iir = pd.read_pickle('Data/old/SWI_ECO_II_2000-2018.PData')
swi_sgr = pd.read_pickle('Data/old/ECOSG_SWI.PData')

'''
f = open('Data/CO2air.txt','r')
f = map(lambda s: s.strip(), f)

t = open('Data/Tair.txt','r')
t = map(lambda s: s.strip(), t)

r = open('Data/Rainf.txt', 'r')
r = map(lambda s: s.strip(), r)

d = open('Data/DIR_SWdown.txt', 'r')
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

### It is required to have as many 'read_pickle's as years in the folder 
#file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_LDAS/observations/sfx-trip/LAI_MC.t08*.PData'),key=numericalSort)
#df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
#      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
#      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
#      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
#      pd.read_pickle(file[15])]
#lai_c4_obs = pd.concat(df)

#file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_SGC/results/ol/LAI_2*.PData'),key=numericalSort)
#df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
#      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
#      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
#      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
#      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
#      pd.read_pickle(file[18])]
#lai_ii_ol_cco2 = pd.concat(df)
'''
file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_LDAS/observations/sfx-trip/LAI_V2_2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
lai_obs = pd.concat(df)

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

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_LDAS/results/ol/LAI_2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
lai_2_ol = pd.concat(df)

file = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_LDAS/results/ekf/LAI_2*.PData'),key=numericalSort)
df = [pd.read_pickle(file[0]),pd.read_pickle(file[1]),pd.read_pickle(file[2]),pd.read_pickle(file[3]),
      pd.read_pickle(file[4]),pd.read_pickle(file[5]),pd.read_pickle(file[6]),pd.read_pickle(file[7]),
      pd.read_pickle(file[8]),pd.read_pickle(file[9]),pd.read_pickle(file[10]),pd.read_pickle(file[11]),
      pd.read_pickle(file[12]),pd.read_pickle(file[13]),pd.read_pickle(file[14]),
      pd.read_pickle(file[15]),pd.read_pickle(file[16]),pd.read_pickle(file[17]),
      pd.read_pickle(file[18])]
lai_2_ekf = pd.concat(df)
'''
#lai_sg = lai_sg_ol
#lai_sg_ol['Analysis'] = lai_sg_ekf.values

###################################################################
### Data Processing 
###################################################################

lai_ii = lai_iir.resample('A',label='left').mean()
lai_ii_c4 = lai_ii_c4r.resample('A',label='left').mean()
lai_sg = lai_sgr.resample('A',label='left').mean()
lai_sg_max = lai_sgr.resample('A',label='left').max()
lai_sg_med = lai_sgr.resample('A',label='left').median()

lai_iiM = lai_iir.resample('M',label='left').mean()
lai_sgM = lai_sgr.resample('M',label='left').mean()

swi_ii = swi_iir.resample('A',label='left').mean()
swi_sg = swi_sgr.resample('A',label='left').mean()

swi_iiM = swi_iir.resample('M',label='left').mean()
swi_sgM = swi_sgr.resample('M',label='left').mean()


### For some reason this .mean(axis=1) converts the df to a series - this is a problem
#lai_c4_obs = lai_c4_obs.mean(axis=1).resample('A',label='left').mean()
#lai_ii_ol_cco2 = lai_ii_ol_cco2.mean(axis=1).resample('A',label='left').mean() 

#Tair = Tair.apply(pd.to_numeric)
#CO2 = CO2.apply(pd.to_numeric)
#Rainf = Rainf.apply(pd.to_numeric)
#SWdown = SWdown.apply(pd.to_numeric)

#CO2 = CO2.resample('A',label='left').mean()
#Tair = Tair.resample('A',label='left').mean()
#Rainf = Rainf.resample('A',label='left').mean()
#SWdown = SWdown.resample('A', label='left').mean()

#sys.exit()
###################################################################
### Computer Z Scores
###################################################################
corn_yz = corn_y.apply(zscore)
lai_ii = lai_ii.apply(zscore)
#lai_ii_c4 = lai_ii_c4.apply(zscore)
lai_sg = lai_sg.apply(zscore)
lai_sg_max = lai_sg_max.apply(zscore)
lai_sg_med = lai_sg_med.apply(zscore)

lai_iiM = lai_iiM.apply(zscore)
lai_sgM = lai_sgM.apply(zscore)


swi_iiM = (swi_iiM - swi_iiM.mean())/swi_iiM.std(ddof=0)
swi_sgM = (swi_sgM - swi_sgM.mean())/swi_sgM.std(ddof=0)


### Having problems with applying zscore to single column dataframes
#lai_c4_obs = lai_c4_obs.apply(zscore)
#lai_ii_ol_cco2 = lai_ii_ol_cco2.apply(zscore)
#lai_sg_ol = lai_sg_ol.apply(zscore)
#lai_sg_ekf = lai_sg_ekf.apply(zscore)

#Tairyz = TairY.apply(zscore)
#tairz = Tair.apply(zscore)
#co2yz = CO2Y.apply(zscore)
#co2z = CO2.apply(zscore)
#rainfyz = RainfY.apply(zscore)
#rainfz = Rainf.apply(zscore)
#swdownyz = SWdownY.apply(zscore)
#swdownz = SWdown.apply(zscore)

lai_ii['Yield'] = corn_yz['Value'].values

#lai_ii_c4['Yield'] = corn_yz['Value'].values

#sys.exit()

### Add C4 (or other) observations into the lai_c4 data frame - because they were imported seperately
#tt=list()
#for i in range(len(lai_ii_c4['Model'])):
#    tt.append(np.nan)

#tt[0:16] = lai_c4_obs.values
#lai_ii_c4['Obs'] = tt

### Convert corn yield into raw numbers in order to make computing correlations easier
corn_yz = corn_yz['Value'].values
#corn_ynz = corn_ynz['Value'].values

###################################################################
### Compute Correlations
###################################################################

#cor_obsc4_yieldn = np.corrcoef(lai_c4_obs_Y.values,corn_ynz[0:16])
#cor_obs_yieldn = np.corrcoef(lai_Y['Obs'].values,corn_ynz) 
#cor_mod_yieldn = np.corrcoef(lai_Y['Model'].values,corn_ynz)
#cor_anl_yieldn = np.corrcoef(lai_Y['Analysis'].values,corn_ynz)


#cor_obsc4_yield = np.corrcoef(lai_c4_obs_Y.values,corn_yz[0:16])
cor_obs_yield = np.corrcoef(lai_ii['Obs'].values,corn_yz)
cor_obs_yield = np.corrcoef(obs.values,corn_yz[0].values)
cor_mod_yield = np.corrcoef(lai_ii['Model'].values,corn_yz)
cor_anl_yield = np.corrcoef(lai_ii['Analysis'].values,corn_yz)
#cor_c4mod_yield = np.corrcoef(lai_c4_Y['Model'].values,corn_yz)
#cor_c4anl_yield = np.corrcoef(lai_c4_Y['Analysis'].values,corn_yz)
cor_mod_ob = np.corrcoef(lai_ii['Model'].values,lai_ii['Obs'].values)
cor_anl_ob = np.corrcoef(lai_ii['Analysis'].values,lai_ii['Obs'].values)


cor_sg_mod_yield = np.corrcoef(lai_sg['Model'].values,corn_yz)
cor_sg_anl_yield = np.corrcoef(lai_sg['Analysis'].values,corn_yz)
cor_sg_mod_ob = np.corrcoef(lai_sg['Model'].values,lai_ii['Obs'].values)
cor_sg_anl_ob = np.corrcoef(lai_sg['Analysis'].values,lai_ii['Obs'].values)

#correlations = pd.DataFrame(columns=['Exp','Cor','P-val','Lower bound','Upper bound'])

#c4Y_c4Obs = pearsonr_ci(lai_c4_obs_Y.values,corn_yz[0:16])
#c4Y_c4Mod = pearsonr_ci(lai_c4_Y['Model'].values,corn_yz)
#c4Y_c4Anl = pearsonr_ci(lai_c4_Y['Analysis'].values,corn_yz)
#c4Y_Obs = pearsonr_ci(lai_Y['Obs'].values,corn_yz)
#c4Y_Mod = pearsonr_ci(lai_Y['Model'].values,corn_yz)
#c4Y_Anl = pearsonr_ci(lai_Y['Analysis'].values,corn_yz)
'''
print("Yearly Observed LAI Anomaly vs Corn Yield Anomaly: " + str(round(cor_obs_yield[0,1],4)))
print("------------------------------------------------------------------")
print("Yearly LAI Open-Loop (ECO_II) Anomaly vs Corn Yield Anomaly: " + str(round(cor_mod_yield[0,1],4)))
print("Yearly LAI Open-Loop (ECO_SG) Anomaly vs Corn Yield Anomaly: " + str(round(cor_sg_mod_yield[0,1],4)))
print("Yearly LAI Open-Loop (ECO_II) Anomaly vs Observed LAI Anomaly: " +str(round(cor_mod_ob[0,1],4)))
print("Yearly LAI Open-Loop (ECO_SG) Anomaly vs Observed LAI Anomaly: " +str(round(cor_sg_mod_ob[0,1],4)))
print("------------------------------------------------------------------")

print("Yearly LAI EKF (ECO_II) Anomaly vs Corn Yield Anomaly: " + str(round(cor_anl_yield[0,1],4)))
print("Yearly LAI EKF (ECO_SG) Anomaly vs Corn Yield Anomaly: " + str(round(cor_sg_anl_yield[0,1],4)))
print("Yearly LAI EKF (ECO_II) Anomaly vs Observed LAI Anomaly: " + str(round(cor_anl_ob[0,1],4)))
print("Yearly LAI EKF (ECO_SG) Anomaly vs Observed LAI Anomaly: " + str(round(cor_sg_anl_ob[0,1],4)))
'''

corswi_ii = swi_iiM.corr()
corswi_sg = swi_sgM.corr()
print("SWI with ECO II")
print(corswi_ii)
print("SWI with ECO SG")
print(corswi_sg)


#sys.exit()
###################################################################
### Graphing Functions 
###################################################################
### plt.subplots(num_rows,num_cols)
#fig, axs = plt.subplots(2,1)

input("Press Enter to continue ...")

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

#plt.title('LAI and Corn Yield Interannual Variability over Nebraska - Variable and Const CO2')
#plt.title('Model LAI, Observed Corn Yield, and Rainfall Anomalies over Nebraska')
#plt.title('LAI from ECO_II vs ECO_SG with Observed Corn Yield over Nebraska')

### Plot Yearly Results ###
fig = plt.figure(figsize=(8,2),dpi=300,edgecolor='w')
ax1 = fig.add_subplot(111)
plt.title('Interannual Variability of LAI and Corn Yield over Nebraska')
#'''
### ONLY YIELD AND OBS
ax1.plot(lai_ii['Yield'],label='Corn Yield - USDA',marker='o',linestyle='--',markersize=6,color='green')
ax1.plot(lai_sg['Obs'],label='LAI Observations',marker='*',linestyle='-',markersize=6,color='black')
'''
### YIELD, OBS, AND ECO_II OL, EKF
plt.plot(lai_ii['Yield'],label='Corn Yield - USDA',marker='o',linestyle='--',markersize=10,color='green')
plt.plot(lai_ii['Obs'],label='LAI Observations',marker='*',linestyle='--',markersize=10,color='black')
#plt.plot(lai_ii['Model'],label='ECO_II Model',color='blue',linewidth=1)
#plt.plot(lai_ii['Analysis'],label='ECO_II Analysis',color='red',linewidth=1)
plt.plot(lai_iiM['Model'],label='ECO_II Model',color='blue',linewidth=1)
#plt.plot(lai_iiM['Analysis'],label='ECO_II Analysis',color='red',linewidth=1)
'''
### YIELD, OBS, AND ECO_SG OL, EKF
#plt.plot(lai_ii['Yield'],label='Corn Yield - USDA',marker='o',linestyle='--',markersize=10,color='green')
#plt.plot(lai_ii['Obs'],label='LAI Observations',marker='*',linestyle='--',markersize=10,color='black')
#plt.plot(lai_sg['Model'],label='Model Mean',color='blue',linewidth=1,linestyle='-')
#plt.plot(lai_sg['Analysis'],label='Analysis Mean',color='red',linewidth=1,linestyle='-')
#plt.plot(lai_sg_med['Model'],label='Model Median',color='blue',linewidth=1,linestyle='--')
#plt.plot(lai_sg_med['Analysis'],label='Analysis Median',color='red',linewidth=1,linestyle='--')
#plt.plot(lai_sgM['Analysis'],label='ECO_SG Analysis',color='red',linewidth=1,linestyle='--')
#plt.plot(lai_sgM['Model'],label='ECO_SG Model',color='blue',linewidth=1,linestyle='--')
#'''



### SWI

#plt.plot(swi_iiM['Obs'],label='SWI Observations',marker='*',linestyle='--',markersize=10,color='black')
#plt.plot(lai_ii['Model'],label='ECO_II Model',color='blue',linewidth=1)
#plt.plot(lai_ii['Analysis'],label='ECO_II Analysis',color='red',linewidth=1)
#plt.plot(swi_iiM['Model'],label='ECO_II Model',color='blue',linewidth=1)
#plt.plot(swi_iiM['Analysis'],label='ECO_II Analysis',color='red',linewidth=1)

### YIELD, OBS, AND ECO_SG OL, EKF
#plt.plot(lai_ii['Yield'],label='Corn Yield - USDA',marker='o',linestyle='--',markersize=10,color='green')
#plt.plot(swi_iiM['Obs'],label='SWI Observations',marker='*',linestyle='--',markersize=10,color='black')
#plt.plot(swi_sgM['Analysis'],label='ECO_SG Analysis',color='red',linewidth=1,linestyle='--')
#plt.plot(swi_sgM['Model'],label='ECO_SG Model',color='blue',linewidth=1,linestyle='--')

#diffModel = swi_iiM['Model'] - swi_sgM['Model']
#diffAna = swi_iiM['Analysis'] - swi_sgM['Analysis']
#diffMLAI = lai_iiM['Model'] - lai_sgM['Model']
#diffALAI = lai_iiM['Analysis'] - lai_sgM['Analysis']

#plt.plot(diffModel,label='ECO II - ECO SG: Model SWI',color='b')
#plt.plot(diffAna,label='ECO II - ECO SG: Analysis SWI',color='r')
#plt.plot(diffMLAI,label='ECO II - ECO SG: Model LAI',color='b',linestyle ='--')
#plt.plot(diffALAI,label='ECO II - ECO SG: Analysis LAI',color='r',linestyle ='--')

### Plot Trendlines ###
#plt.plot(lai_Y['t_trend'],label='Temperature Trend',linestyle="--",color='red')
#plt.plot(lai_Y['l_trend'],label='OL LAI Trend',linestyle="--",color='blue')
#plt.plot(lai_Y['o_trend'],label='Corn Yield Trend',linestyle="--",color='green')

ax1.axhline(color='black',linewidth=.5)
plt.legend(loc='upper left',fontsize=6,numpoints=1)
plt.ylim([-2.5,2.5])
#plt.ylim([-1,1])
plt.ylabel('Anomaly')
#plt.xlim('2017-01-31','2017-12-31')

plt.rcParams.update({'font.size': 8})
#fig.tight_layout()
#save('../vegeoFigures/Yield-RainF', ext='ps', close=True, verbose=True)
plt.show()

###################################################################
### Graphing Functions 
###################################################################
'''
#v1 = lai_sg_ol.mean(axis=0).values.reshape(16,40)
#v2 = lai_obs.mean(axis=0).values.reshape(16,40)
#v3 = lai_sg_ekf.mean(axis=0).values.reshape(16,40)
#v4 = v3-v1

v1 = lai_2_ol.mean(axis=0).values.reshape(16,40)
v2 = lai_obs.mean(axis=0).values.reshape(16,40)
v3 = lai_2_ekf.mean(axis=0).values.reshape(16,40)
v4 = v3-v1

#v5 = c4_SG.reshape(16,40)

fig, axes = plt.subplots(1,4)

ax1 = plt.subplot2grid((1,4), (0,0))
ax2 = plt.subplot2grid((1,4), (0,1))
ax3 = plt.subplot2grid((1,4), (0,2))
ax4 = plt.subplot2grid((1,4), (0,3))

#ax5 = plt.subplot2grid((2,4), (1,0))
#ax6 = plt.subplot2grid((2,4), (1,1))
#ax7 = plt.subplot2grid((2,4), (1,2))
#ax8 = plt.subplot2grid((2,4), (1,3))

#axes[0].set_title("C4 Fraction - 25km")
ax1.set_title("OL")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax=ax1)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1),linewidth=.5)
map.drawparallels(np.array([40.5,42.5]),labels=(1,0,0,1),linewidth=.5)
map.drawrivers()
map.imshow(v1,interpolation='none',cmap='RdYlGn',vmin=0.,vmax=2.)

ax2.set_title("Obs")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax2)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.array([40.5,42.5]),labels=(1,0,0,1),linewidth=.5)
map.drawrivers()
map.imshow(v2,interpolation='none',cmap='RdYlGn',vmin=0.,vmax=2.)

#diff = v1-v2
ax3.set_title("EKF")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax3)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.array([40.5,42.5]),labels=(1,0,0,1),linewidth=.5)
map.drawrivers()
map.imshow(v3,interpolation='none',cmap='RdYlGn',vmin=0.,vmax=2.)
#map.imshow(v3,interpolation='none',cmap='seismic',vmin=-0.5,vmax=0.5)

ax4.set_title("EKF-OL")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax4)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.array([40.5,42.5]),labels=(1,0,0,1),linewidth=.5)
map.drawrivers()
im = map.imshow(v4,interpolation='none',cmap='RdBu_r',vmin=-0.75,vmax=0.75)


ax5.set_title("OL")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax5)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.array([40.5,42.5]),labels=(1,0,0,1),linewidth=.5)
map.drawrivers()
map.imshow(v5,interpolation='none',cmap='RdYlGn',vmin=0.,vmax=2.)

ax6.set_title("Obs")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax6)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.array([40.5,42.5]),labels=(1,0,0,1),linewidth=.5)
map.drawrivers()
map.imshow(v6,interpolation='none',cmap='RdYlGn',vmin=0.,vmax=2.)

ax7.set_title("OL")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax7)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.array([40.5,42.5]),labels=(1,0,0,1),linewidth=.5)
map.drawrivers()
map.imshow(v7,interpolation='none',cmap='RdYlGn',vmin=0.,vmax=2.)

ax8.set_title("EKF-OL")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = ax8)
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.array([40.5,42.5]),labels=(1,0,0,1),linewidth=.5)
map.drawrivers()
map.imshow(v8,interpolation='none',cmap='RdBu_r',vmin=-0.75,vmax=0.75)


### for 'add_axes' it is [left,bottom,witdth,height], plt.figure(#) is seperate image, or on the same plot
#fig = plt.figure(0)
#cbaxes = fig.add_axes([0.08,0.85,0.7,0.025])
#fig.colorbar(im,ax=axes.ravel().tolist())
#cbaxes = matplotlib.colorbar.make_axes(location='bottom')
#cbar = fig.colorbar(im,cax=cbaxes,orientation='horizontal',ticks=[-.75,0,.75])
#cbar.ax.set_xticklabels(['-0.75','0','0.75'])


#save('../vegeoFigures/Yield-RainF', ext='ps', close=True, verbose=True)
plt.show()
'''
