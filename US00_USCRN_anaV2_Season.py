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


corfiles = sorted(glob.glob('C:/Users/antho/Documents/SXVGO1/CustomPostProcessing/Data/US00/USCRN/Seasonal/*R.PData'),key=numericalSort)
pfiles = sorted(glob.glob('C:/Users/antho/Documents/SXVGO1/CustomPostProcessing/Data/US00/USCRN/Seasonal/*P.PData'),key=numericalSort)

alpha = float(input("Alpha? : "))

names = [];cor = [];pv=[];ci =[]
WG3 = [] ; WG_20 = [] ; WG6 = [] ; WG8 = []
cWG3 = [] ; cWG_20 = [] ; cWG6 = [] ; cWG8 = []
sWG3 = [] ; sWG_20 = [] ; sWG6 = [] ; sWG8 = []
pWG3 = [] ; pWG_20 = [] ; pWG6 = [] ; pWG8 = []
uWG3 = [] ; uWG_20 = [] ; uWG6 = [] ; uWG8 = []
lWG3 = [] ; lWG_20 = [] ; lWG6 = [] ; lWG8 = []
mWG3 = [] ; mWG_20 = [] ; mWG6 = [] ; mWG8 = []
sprWG3 = [];sprWG_20 = [];sprWG6 = [];sprWG8 = [];
sumWG3 = [];sumWG_20 = [];sumWG6 = [];sumWG8 = [];
autWG3 = [];autWG_20 = [];autWG6 = [];autWG8 = [];
winWG3 = [];winWG_20 = [];winWG6 = [];winWG8 = [];
growWG3 = [];growWG_20 = [];growWG6 = [];growWG8 = [];
sprcWG3 = [];sprcWG_20 = [];sprcWG6 = [];sprcWG8 = [];
sumcWG3 = [];sumcWG_20 = [];sumcWG6 = [];sumcWG8 = [];
autcWG3 = [];autcWG_20 = [];autcWG6 = [];autcWG8 = [];
wincWG3 = [];wincWG_20 = [];wincWG6 = [];wincWG8 = [];
growcWG3 = [];growcWG_20 = [];growcWG6 = [];growcWG8 = [];

stats = []
pmin=[]

for i in range(len(corfiles)):
    #exp = corfiles[i][84:-13]
    #exp = corfiles[i][69:-13]
    exp = corfiles[i][78:-8]
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
        if 'spr' in exp:
            sprWG3.append(exp)
            sprcWG3.append(f.mean()[0])
        if 'sum' in exp:
            sumWG3.append(exp)
            sumcWG3.append(f.mean()[0])
        if 'aut' in exp:
            autWG3.append(exp)
            autcWG3.append(f.mean()[0])
        if 'win' in exp:
            winWG3.append(exp)
            wincWG3.append(f.mean()[0])
        if 'grow' in exp:
            growWG3.append(exp)
            growcWG3.append(f.mean()[0])

    elif 'WG_20' in exp:
        if 'spr' in exp:
            sprWG_20.append(exp)
            sprcWG_20.append(f.mean()[0])
        if 'sum' in exp:
            sumWG_20.append(exp)
            sumcWG_20.append(f.mean()[0])
        if 'aut' in exp:
            autWG_20.append(exp)
            autcWG_20.append(f.mean()[0])
        if 'win' in exp:
            winWG_20.append(exp)
            wincWG_20.append(f.mean()[0])
        if 'grow' in exp:
            growWG_20.append(exp)
            growcWG_20.append(f.mean()[0])

    elif 'WG6' in exp:
        if 'spr' in exp:
            sprWG6.append(exp)
            sprcWG6.append(f.mean()[0])
        if 'sum' in exp:
            sumWG6.append(exp)
            sumcWG6.append(f.mean()[0])
        if 'aut' in exp:
            autWG6.append(exp)
            autcWG6.append(f.mean()[0])
        if 'win' in exp:
            winWG6.append(exp)
            wincWG6.append(f.mean()[0])
        if 'grow' in exp:
            growWG6.append(exp)
            growcWG6.append(f.mean()[0])

    elif 'WG8' in exp:
        if 'spr' in exp:
            sprWG8.append(exp)
            sprcWG8.append(f.mean()[0])
        if 'sum' in exp:
            sumWG8.append(exp)
            sumcWG8.append(f.mean()[0])
        if 'aut' in exp:
            autWG8.append(exp)
            autcWG8.append(f.mean()[0])
        if 'win' in exp:
            winWG8.append(exp)
            wincWG8.append(f.mean()[0])
        if 'grow' in exp:
            growWG8.append(exp)
            growcWG8.append(f.mean()[0])

    print(exp)

sprWG3[:] = (elem[:-9] for elem in sprWG3)
sprWG_20[:] = (elem[:-9] for elem in sprWG_20)
sprWG6[:] = (elem[:-9] for elem in sprWG6)
sprWG8[:] = (elem[:-9] for elem in sprWG8)

order = [6,2,0,1,4,3,5]
#l3 = list(zip(WG3,cWG3))
#l3_r = [l3[i] for i in order]
#wg3 = pd.DataFrame(l3_r,columns=['Exp','Correlation'])

sprl3 = list(zip(sprWG3,sprcWG3))
sprl3_r = [sprl3[i] for i in order]
sprwg3 = pd.DataFrame(sprl3_r,columns=['Exp','Correlation'])

suml3 = list(zip(sumWG3,sumcWG3))
suml3_r = [suml3[i] for i in order]
sumwg3 = pd.DataFrame(suml3_r,columns=['Exp','Correlation'])

autl3 = list(zip(autWG3,autcWG3))
autl3_r = [autl3[i] for i in order]
autwg3 = pd.DataFrame(autl3_r,columns=['Exp','Correlation'])

winl3 = list(zip(winWG3,wincWG3))
winl3_r = [winl3[i] for i in order]
winwg3 = pd.DataFrame(winl3_r,columns=['Exp','Correlation'])

growl3 = list(zip(growWG3,growcWG3))
growl3_r = [growl3[i] for i in order]
growwg3 = pd.DataFrame(growl3_r,columns=['Exp','Correlation'])

wg3 = pd.DataFrame({'Exp':sprwg3['Exp'],
                    'Spr':sprwg3['Correlation'].values,
                    'Sum':sumwg3['Correlation'].values,
                    'Aut':autwg3['Correlation'].values,
                    'Win':winwg3['Correlation'].values,
                    'Grow':growwg3['Correlation'].values})

sprl20 = list(zip(sprWG_20,sprcWG_20))
sprl20_r = [sprl20[i] for i in order]
sprwg20 = pd.DataFrame(sprl20_r,columns=['Exp','Correlation'])

suml20 = list(zip(sumWG_20,sumcWG_20))
suml20_r = [suml20[i] for i in order]
sumwg20 = pd.DataFrame(suml20_r,columns=['Exp','Correlation'])

autl20 = list(zip(autWG_20,autcWG_20))
autl20_r = [autl20[i] for i in order]
autwg20 = pd.DataFrame(autl20_r,columns=['Exp','Correlation'])

winl20 = list(zip(winWG_20,wincWG_20))
winl20_r = [winl20[i] for i in order]
winwg20 = pd.DataFrame(winl20_r,columns=['Exp','Correlation'])

growl20 = list(zip(growWG_20,growcWG_20))
growl20_r = [growl20[i] for i in order]
growwg20 = pd.DataFrame(growl20_r,columns=['Exp','Correlation'])

wg20 = pd.DataFrame({'Exp':sprwg20['Exp'],
                    'Spr':sprwg20['Correlation'].values,
                    'Sum':sumwg20['Correlation'].values,
                    'Aut':autwg20['Correlation'].values,
                    'Win':winwg20['Correlation'].values,
                    'Grow':growwg20['Correlation'].values})

sprl6 = list(zip(sprWG6,sprcWG6))
sprl6_r = [sprl6[i] for i in order]
sprwg6 = pd.DataFrame(sprl6_r,columns=['Exp','Correlation'])

suml6 = list(zip(sumWG6,sumcWG6))
suml6_r = [suml6[i] for i in order]
sumwg6 = pd.DataFrame(suml6_r,columns=['Exp','Correlation'])

autl6 = list(zip(autWG6,autcWG6))
autl6_r = [autl6[i] for i in order]
autwg6 = pd.DataFrame(autl6_r,columns=['Exp','Correlation'])

winl6 = list(zip(winWG6,wincWG6))
winl6_r = [winl6[i] for i in order]
winwg6 = pd.DataFrame(winl6_r,columns=['Exp','Correlation'])

growl6 = list(zip(growWG6,growcWG6))
growl6_r = [growl6[i] for i in order]
growwg6 = pd.DataFrame(growl6_r,columns=['Exp','Correlation'])

wg6 = pd.DataFrame({'Exp':sprwg6['Exp'],
                    'Spr':sprwg6['Correlation'].values,
                    'Sum':sumwg6['Correlation'].values,
                    'Aut':autwg6['Correlation'].values,
                    'Win':winwg6['Correlation'].values,
                    'Grow':growwg6['Correlation'].values})

sprl8 = list(zip(sprWG8,sprcWG8))
sprl8_r = [sprl8[i] for i in order]
sprwg8 = pd.DataFrame(sprl8_r,columns=['Exp','Correlation'])

suml8 = list(zip(sumWG8,sumcWG8))
suml8_r = [suml8[i] for i in order]
sumwg8 = pd.DataFrame(suml8_r,columns=['Exp','Correlation'])

autl8 = list(zip(autWG8,autcWG8))
autl8_r = [autl8[i] for i in order]
autwg8 = pd.DataFrame(autl8_r,columns=['Exp','Correlation'])

winl8 = list(zip(winWG8,wincWG8))
winl8_r = [winl8[i] for i in order]
winwg8 = pd.DataFrame(winl8_r,columns=['Exp','Correlation'])

growl8 = list(zip(growWG8,growcWG8))
growl8_r = [growl8[i] for i in order]
growwg8 = pd.DataFrame(growl8_r,columns=['Exp','Correlation'])

wg8 = pd.DataFrame({'Exp':sprwg8['Exp'],
                    'Spr':sprwg8['Correlation'].values,
                    'Sum':sumwg8['Correlation'].values,
                    'Aut':autwg8['Correlation'].values,
                    'Win':winwg8['Correlation'].values,
                    'Grow':growwg8['Correlation'].values})

wg3 = wg3.round(3)
wg20 = wg20.round(3)
wg6 = wg6.round(3)
wg8 = wg8.round(3)
'''
ekf_lai = [wg3["Correlation"][0], wg20["Correlation"][0],wg6["Correlation"][0],wg8["Correlation"][0]]
ekf_lai_ssm = [wg3["Correlation"][1], wg20["Correlation"][1],wg6["Correlation"][1],wg8["Correlation"][1]]
ekf_ssm = [wg3["Correlation"][2], wg20["Correlation"][2],wg6["Correlation"][2],wg8["Correlation"][2]]
ekf_vod = [wg3["Correlation"][3], wg20["Correlation"][3],wg6["Correlation"][3],wg8["Correlation"][3]]
ekf_vod_ssm = [wg3["Correlation"][4], wg20["Correlation"][4],wg6["Correlation"][4],wg8["Correlation"][4]]
ol = [wg3["Correlation"][5], wg20["Correlation"][5],wg6["Correlation"][5],wg8["Correlation"][5]]
'''
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

### Graphing/Plotting 

#fig, axes = plt.subplots(4,1)
#fig = plt.figure(figsize=(3,2),dpi=150,edgecolor='w')
fig = plt.figure()
axes = fig.add_subplot(111)

plt.title('USCRN In Situ SM vs Assimilation Scenarios correlations : Growing Season')

#x = 0,1,2,3,4,5,6
#x = 6,2,0,1,4,3,5
#ticks = ['EKF LAI','EKF LAI+SSM','EKF SSM','EKF VOD10','EKF VOD','EKF VOD SSM','OL']
ticks = ['OL','EKF SSM','EKF LAI','EKF LAI+SSM','EKF VOD','EKF VOD10','EKF VOD+SSM']
plt.ylabel('Correlation')

#plt.xticks(np.arange(4))
#plt.xticklabels(ticks)
#plt.xlim(-0.5,3.5)

axes.set_xticks(np.arange(7))
axes.set_xticklabels(ticks)
axes.set_xlim(-0.5,6.5)
axes.set_ylim(0.3,0.8)

axes.plot(wg3["Grow"],color='b',marker='h',linewidth=1,label='5cm',markeredgecolor='black')
axes.plot(wg20["Grow"],color='r',marker='h',linewidth=1,label='20cm',markeredgecolor='black',linestyle='dashed')
axes.plot(wg6["Grow"],color='brown',marker='h',linewidth=1,label='50cm',markeredgecolor='black',linestyle='dotted')
axes.plot(wg8["Grow"],color='green',marker='h',linewidth=1,label='100cm',markeredgecolor='black',linestyle='dashdot')

fig.tight_layout()
#plt.legend(framealpha=1)
plt.grid(True)
plt.show()

### Seasonal Graphs 2x2
fig, axes = plt.subplots(2,2)
#plt.title('USCRN In Situ SM vs Assimilation Scenarios correlations - Seasonal')
ticks = ['OL','EKF SSM','EKF LAI','EKF LAI+SSM','EKF VOD','EKF VOD10','EKF VOD+SSM']

### Spring
axes[0][0].set_title("Spring")
axes[0][0].set_ylabel('Correlation')
axes[0][0].set_xticks(np.arange(7))
axes[0][0].set_xticklabels(ticks)
axes[0][0].set_xlim(-0.5,6.5)
axes[0][0].set_ylim(0.3,0.8)
axes[0][0].grid(True)
axes[0][0].plot(wg3["Spr"],color='b',marker='h',linewidth=1,label='5cm',markeredgecolor='black')
axes[0][0].plot(wg20["Spr"],color='r',marker='h',linewidth=1,label='20cm',markeredgecolor='black',linestyle='dashed')
axes[0][0].plot(wg6["Spr"],color='brown',marker='h',linewidth=1,label='50cm',markeredgecolor='black',linestyle='dotted')
axes[0][0].plot(wg8["Spr"],color='green',marker='h',linewidth=1,label='100cm',markeredgecolor='black',linestyle='dashdot')

### Summer
axes[0][1].set_title("Summer")
axes[0][1].set_ylabel('Correlation')
axes[0][1].set_xticks(np.arange(7))
axes[0][1].set_xticklabels(ticks)
axes[0][1].set_xlim(-0.5,6.5)
axes[0][1].set_ylim(0.3,0.8)
axes[0][1].grid(True)
axes[0][1].plot(wg3["Sum"],color='b',marker='h',linewidth=1,label='5cm',markeredgecolor='black')
axes[0][1].plot(wg20["Sum"],color='r',marker='h',linewidth=1,label='20cm',markeredgecolor='black',linestyle='dashed')
axes[0][1].plot(wg6["Sum"],color='brown',marker='h',linewidth=1,label='50cm',markeredgecolor='black',linestyle='dotted')
axes[0][1].plot(wg8["Sum"],color='green',marker='h',linewidth=1,label='100cm',markeredgecolor='black',linestyle='dashdot')

### Fall
axes[1][0].set_title("Autumn")
axes[1][0].set_ylabel('Correlation')
axes[1][0].set_xticks(np.arange(7))
axes[1][0].set_xticklabels(ticks)
axes[1][0].set_xlim(-0.5,6.5)
axes[1][0].set_ylim(0.3,0.8)
axes[1][0].grid(True)
axes[1][0].plot(wg3["Aut"],color='b',marker='h',linewidth=1,label='5cm',markeredgecolor='black')
axes[1][0].plot(wg20["Aut"],color='r',marker='h',linewidth=1,label='20cm',markeredgecolor='black',linestyle='dashed')
axes[1][0].plot(wg6["Aut"],color='brown',marker='h',linewidth=1,label='50cm',markeredgecolor='black',linestyle='dotted')
axes[1][0].plot(wg8["Aut"],color='green',marker='h',linewidth=1,label='100cm',markeredgecolor='black',linestyle='dashdot')

### Winter
axes[1][1].set_title("Winter")
axes[1][1].set_ylabel('Correlation')
axes[1][1].set_xticks(np.arange(7))
axes[1][1].set_xticklabels(ticks)
axes[1][1].set_xlim(-0.5,6.5)
axes[1][1].set_ylim(0.3,0.8)
axes[1][1].grid(True)
axes[1][1].plot(wg3["Win"],color='b',marker='h',linewidth=1,label='5cm',markeredgecolor='black')
axes[1][1].plot(wg20["Win"],color='r',marker='h',linewidth=1,label='20cm',markeredgecolor='black',linestyle='dashed')
axes[1][1].plot(wg6["Win"],color='brown',marker='h',linewidth=1,label='50cm',markeredgecolor='black',linestyle='dotted')
axes[1][1].plot(wg8["Win"],color='green',marker='h',linewidth=1,label='100cm',markeredgecolor='black',linestyle='dashdot')

plt.legend(framealpha=1,bbox_to_anchor = (-0.04,1))
#plt.figlegend()
#fig.legend()
#plt.legend()
plt.show()


sys.exit()
### NIC Analysis
R_mod = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/ol_WG_20_0.20_R.PData')
R_ekf = pd.read_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/USCRN/ekf_ssm_WG_20_0.20_R.PData')

NIC = 100*(np.array(R_ekf)-np.array(R_mod))/(1-np.array(R_mod))

w1=np.where(NIC>3)[0]
w3=np.where((NIC>-3) & (NIC<3))[0]
w2=np.where(NIC<-3)[0]

print("# of stations with Improved correlation : {0}".format(len(w1)))
print("# of stations with Degraded correlation : {0}".format(len(w2)))
print("# of stations with No Improvement or Degredation: {0}".format(len(w3)))

