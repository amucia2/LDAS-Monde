import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#from mpl_toolkits.basemap import Basemap
import pandas as pd
#from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr
import pylab as pl
import sys, glob, os, re
import scipy.io as spio 
os.chdir('C:/Users/antho/Documents/SXVGO1/ldas_curr')
from ldasBoot import *
os.chdir('C:/Users/antho/Documents/SXVGO1/CustomPostProcessing/')

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts


corfiles = sorted(glob.glob('C:/Users/antho/Documents/SXVGO1/CustomPostProcessing/Data/US00/USCRN/*R.PData'),key=numericalSort)
pfiles = sorted(glob.glob('C:/Users/antho/Documents/SXVGO1/CustomPostProcessing/Data/US00/USCRN/*P.PData'),key=numericalSort)

alpha = float(input("Alpha? : "))

names = [];cor = [];pv=[];ci =[]
WG3 = [] ; WG_20 = [] ; WG6 = [] ; WG8 = []
cWG3 = [] ; cWG_20 = [] ; cWG6 = [] ; cWG8 = []
sWG3 = [] ; sWG_20 = [] ; sWG6 = [] ; sWG8 = []
pWG3 = [] ; pWG_20 = [] ; pWG6 = [] ; pWG8 = []
uWG3 = [] ; uWG_20 = [] ; uWG6 = [] ; uWG8 = []
lWG3 = [] ; lWG_20 = [] ; lWG6 = [] ; lWG8 = []
mWG3 = [] ; mWG_20 = [] ; mWG6 = [] ; mWG8 = []
stats = []
pmin=[]

for i in range(len(corfiles)):
    #exp = corfiles[i][84:-13]
    exp = corfiles[i][69:-13]
    names.append(exp)
    f = pd.read_pickle(corfiles[i])
    p = pd.read_pickle(pfiles[i])
    cor.append(f.mean()[0])
    pv.append(p.mean()[0])
    #pmin.append(min(p)[0])
    stats.append(111-f.isna().sum()[0])
    R = f.values
    upp = bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000,alpha=alpha)[1]
    low = bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000,alpha=alpha)[0]
    ci.append((abs(f.mean()[0]-low)+abs(f.mean()[0]-upp))/2)
    if 'WG3' in exp:
        WG3.append(exp)
        cWG3.append(f.mean()[0])
        sWG3.append(111-f.isna().sum()[0])
        pWG3.append(p.mean()[0])
        upp = bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000,alpha=alpha)[1]
        low = bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000,alpha=alpha)[0]
        mWG3.append((abs(f.mean()[0]-low)+abs(f.mean()[0]-upp))/2)
    elif 'WG_20' in exp:
        WG_20.append(exp)
        cWG_20.append(f.mean()[0])
        sWG_20.append(111-f.isna().sum()[0])
        pWG_20.append(p.mean()[0])
        upp = bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000,alpha=alpha)[1]
        low = bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000,alpha=alpha)[0]
        mWG_20.append((abs(f.mean()[0]-low)+abs(f.mean()[0]-upp))/2)
    elif 'WG6' in exp:
        WG6.append(exp)
        cWG6.append(f.mean()[0])
        sWG6.append(111-f.isna().sum()[0])
        pWG6.append(p.mean()[0])
        upp = bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000,alpha=alpha)[1]
        low = bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000,alpha=alpha)[0]
        mWG6.append((abs(f.mean()[0]-low)+abs(f.mean()[0]-upp))/2)
    elif 'WG8' in exp:
        WG8.append(exp)
        cWG8.append(f.mean()[0])  
        sWG8.append(111-f.isna().sum()[0])
        pWG8.append(p.mean()[0])
        upp = bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000,alpha=alpha)[1]
        low = bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000,alpha=alpha)[0]
        mWG8.append((abs(f.mean()[0]-low)+abs(f.mean()[0]-upp))/2)
    print(exp)

#order = [6,2,0,1,4,3,5]
order = [6,0,4,3,1,2,5]
order = [6,0,3,4,1,2,5]
l = list(zip(names,cor,ci,stats))
df = pd.DataFrame(l,columns=['Exp','Correlation','{0} CI +/-'.format(1-alpha),'# of Stations'])

l3 = list(zip(WG3,cWG3,mWG3,sWG3))
l3_r = [l3[i] for i in order]
wg3 = pd.DataFrame(l3_r,columns=['Exp','Correlation','{0} CI +/-'.format(1-alpha),'# of Stations'])

l20 = list(zip(WG_20,cWG_20,mWG_20,sWG_20))
l20_r = [l20[i] for i in order]
wg20 = pd.DataFrame(l20_r,columns=['Exp','Correlation','{0} CI +/-'.format(1-alpha),'# of Stations'])

l6 = list(zip(WG6,cWG6,mWG6,sWG6))
l6_r = [l6[i] for i in order]
wg6 = pd.DataFrame(l6_r,columns=['Exp','Correlation','{0} CI +/-'.format(1-alpha),'# of Stations'])

l8 = list(zip(WG8,cWG8,mWG8,sWG8))
l8_r = [l8[i] for i in order]
wg8 = pd.DataFrame(l8_r,columns=['Exp','Correlation','{0} CI +/-'.format(1-alpha),'# of Stations'])

wg3 = wg3.round(3)
wg20 = wg20.round(3)
wg6 = wg6.round(3)
wg8 = wg8.round(3)

### For only VOD output
'''
wg3 = wg3.drop([4,5,6])
wg20 = wg20.drop([4,5,6])
wg6 = wg6.drop([4,5,6])
wg8 = wg8.drop([4,5,6])

'''
ekf_lai = [wg3["Correlation"][0], wg20["Correlation"][0],wg6["Correlation"][0],wg8["Correlation"][0]]
ekf_lai_ssm = [wg3["Correlation"][1], wg20["Correlation"][1],wg6["Correlation"][1],wg8["Correlation"][1]]
ekf_ssm = [wg3["Correlation"][2], wg20["Correlation"][2],wg6["Correlation"][2],wg8["Correlation"][2]]
ekf_vod = [wg3["Correlation"][3], wg20["Correlation"][3],wg6["Correlation"][3],wg8["Correlation"][3]]
ekf_vod_ssm = [wg3["Correlation"][4], wg20["Correlation"][4],wg6["Correlation"][4],wg8["Correlation"][4]]
ol = [wg3["Correlation"][5], wg20["Correlation"][5],wg6["Correlation"][5],wg8["Correlation"][5]]

print("Soil Layer WG3")
print(wg3)
print("Soil Layer WG4-WG5 (20cm)")
print(wg20)
print("Soil Layer WG6")
print(wg6)
print("Soil Layer WG8")
print(wg8)

print("Computation Complete")   
sys.exit()

#wg3 = wg3.reindex([6,2,0,1,4,3,5])
#wg20 = wg20.reindex([6,2,0,1,4,3,5])
#wg6 = wg6.reindex([6,2,0,1,4,3,5])
#wg8 = wg8.reindex([6,2,0,1,4,3,5])

#wg3 = wg3.loc[[6,2,0,1,4,3,5]]
#wg20 = wg20.loc[[6,2,0,1,4,3,5]]
#wg6 = wg6.loc[[6,2,0,1,4,3,5]]
#wg8 = wg8.loc[[6,2,0,1,4,3,5]]

### Graphing/Plotting 

#fig, axes = plt.subplots(4,1)
#fig = plt.figure(figsize=(3,2),dpi=150,edgecolor='w')
fig = plt.figure()
axes = fig.add_subplot(111)

plt.title('USCRN In Situ SM vs Assimilation Scenarios correlations')

#x = 0,1,2,3,4,5,6
#x = 6,2,0,1,4,3,5
#ticks = ['EKF LAI','EKF LAI+SSM','EKF SSM','EKF VOD10','EKF VOD','EKF VOD SSM','OL']
ticks = ['OL','EKF SSM','EKF LAI','EKF LAI+SSM','EKF VOD','EKF VOD10','EKF VOD+SSM']
ticks = ['OL','EKF LAI','EKF VOD','EKF VOD10']
ticks = ['OL','EKF LAI','EKF VOD10','EKF VOD']
plt.ylabel('Correlation')

#plt.xticks(np.arange(4))
#plt.xticklabels(ticks)
#plt.xlim(-0.5,3.5)

axes.set_xticks(np.arange(4))
axes.set_xticklabels(ticks)
axes.set_xlim(-0.5,3.5)

axes.plot(wg3["Correlation"],color='b',marker='h',linewidth=1,label='5cm',markeredgecolor='black')
axes.plot(wg20["Correlation"],color='r',marker='h',linewidth=1,label='20cm',markeredgecolor='black',linestyle='dashed')
axes.plot(wg6["Correlation"],color='brown',marker='h',linewidth=1,label='50cm',markeredgecolor='black',linestyle='dotted')
axes.plot(wg8["Correlation"],color='green',marker='h',linewidth=1,label='100cm',markeredgecolor='black',linestyle='dashdot')

fig.tight_layout()
#plt.legend(framealpha=1)
plt.grid(True)
plt.show()

'''
fig, axes = plt.subplots(2,2)

plt.xlabel('Correlation')

axes[0][0].plot(wg3["Correlation"],color='b',marker='h',linewidth=1,label='5cm',markeredgecolor='black')
axes[0][1].plot(wg20["Correlation"],color='r',marker='h',linewidth=1,label='20cm',markeredgecolor='black')
axes[1][0].plot(wg6["Correlation"],color='r',marker='h',linewidth=1,label='50cm',markeredgecolor='black')
axes[1][1].plot(wg8["Correlation"],color='r',marker='h',linewidth=1,label='100cm',markeredgecolor='black')

fig.tight_layout()
plt.show()
'''

### Mapping taken from comp_SMN python files


R_mod = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/ol_WG3_0.05_R.PData')
R_ekf = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/ekf_vod10_WG3_0.05_R.PData')

R_mod = pd.read_pickle('C:/Users/antho/Documents/SXVGO1/CustomPostProcessing/Data/US00/USCRN/ol_WG8_1.00_R.PData')
R_ekf = pd.read_pickle('C:/Users/antho/Documents/SXVGO1/CustomPostProcessing/Data/US00/USCRN/ekf_vod10_WG8_1.00_R.PData')

NIC = 100*(np.array(R_ekf)-np.array(R_mod))/(1-np.array(R_mod))

w1=np.where(NIC>3)[0]
w3=np.where((NIC>-3) & (NIC<3))[0]
w2=np.where(NIC<-3)[0]

print("# of stations with Improved correlation : {0}".format(len(w1)))
print("# of stations with Degraded correlation : {0}".format(len(w2)))
print("# of stations with No Improvement or Degredation: {0}".format(len(w3)))

#P_mod = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/ol_WG3_0.05_P.PData')
#P_ekf = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/ekf_lai_WG3_0.05_P.PData')

#lat = np.loadtxt('/cnrm/vegeo/muciaa/NO_SAVE/US00/USCRN_coords_wban.txt')[:,0]
#lon = np.loadtxt('/cnrm/vegeo/muciaa/NO_SAVE/US00/USCRN_coords_wban.txt')[:,1]

w = np.where(P_mod<0.05 and P_ekf<0.05)
fig, axes = plt.subplots(1,1)
axes.set_title("Evaluation vs USCRN")
map = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-55, urcrnrlon=180, urcrnrlat=90,resolution='c',ax = axes)
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines()
map.drawmapboundary(fill_color='lightblue')
map.drawcountries()
map.fillcontinents(color='beige',lake_color='lightblue')
x,y = map(np.array(lon), np.array(lat))
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
x,y = map(np.array(lon)[w3], np.array(lat)[w3])
cs2 = map.scatter(x,y,marker='v',s=25,c=np.array(NIC[w3]),vmin=-50,vmax=50,cmap=cm,zorder=2)
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
