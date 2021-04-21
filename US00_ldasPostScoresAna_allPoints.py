import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.mlab as mla
#from mpl_toolkits.basemap import Basemap
import pandas as pd
#from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr
import pylab as pl
import sys, glob, os, re
import scipy.io as spio 
import pickle
import xarray as xr
#os.chdir('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr')
#from ldasBoot import *
#os.chdir('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing')

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts



files = sorted(glob.glob('C:/Users/antho/Documents/SXVGO1/CustomPostProcessing/Data/US00/ldasPostScores/V4/SPlains*scores_for_all_points*'),key=numericalSort)
#files = sorted(glob.glob('C:/Users/antho/Documents/SXVGO1/CustomPostProcessing/Data/US00/ldasPostScores/V4/California*per_point*'),key=numericalSort)

for i in range(len(files)):
    #[93:-3] when using "scores_for_all_points" on sxvgo1
    #[78:-3] when on home pc
    exp = files[i][78:-3]
    print(exp)
    
    if 'LAI_V2' in exp:
        df = xr.open_dataarray(files[i])
        ### Opens all stats in the "Model" place
        if '_LAI_VOD_' in exp:
            lai_lai_R = df[1][0].values
        if '_LAISSM_' in exp:
            lai_laissm_R = df[1][0].values
        if '_VOD10_' in exp:
            lai_vod10_R = df[1][0].values
        if '_OL_' in exp:
            lai_ol_R = df[1][0].values
        ### Opens all stats in the "Analysis" place
        if '_VOD_' in exp:
            lai_vod_R = df[1][1].values
        if '_VODSSM_' in exp:
            lai_vodssm_R = df[1][1].values
        if '_SSM_' in exp:
            lai_ssm_R = df[1][1].values
    ''' 
    if 'LAI_V2' in exp:
        df = xr.open_dataarray(files[i])
        ### Opens all stats in the "Model" place
        if '_LAI_' in exp:
            lai_lai_R = np.nanmean(df[1][0].values,axis=1)
        if '_LAISSM_' in exp:
            lai_laissm_R = np.nanmean(df[1][0].values,axis=1)
        if '_VOD10_' in exp:
            lai_vod10_R = np.nanmean(df[1][0].values,axis=1)
        if '_OL_' in exp:
            lai_ol_R = np.nanmean(df[1][0].values,axis=1)
        ### Opens all stats in the "Analysis" place
        if '_VOD_' in exp:
            lai_vod_R = np.nanmean(df[1][1].values,axis=1)
        if '_VODSSM_' in exp:
            lai_vodssm_R = np.nanmean(df[1][1].values,axis=1)
        if '_SSM_' in exp:
            lai_ssm_R = np.nanmean(df[1][1].values,axis=1)
    '''
    if 'FLCGPP' in exp:
        df = xr.open_dataarray(files[i])
        if '_LAI_' in exp:
            gpp_lai_R = df[1][0].values
        if '_LAISSM_' in exp:
            gpp_laissm_R = df[1][0].values
        if '_VOD10_' in exp:
            gpp_vod10_R = df[1][0].values
        if '_OL_' in exp:
            gpp_ol_R = df[1][0].values
        if '_VOD_' in exp:
            gpp_vod_R = df[1][1].values
        if '_VODSSM_' in exp:
            gpp_vodssm_R = df[1][1].values
        if '_SSM_' in exp:
            gpp_ssm_R = df[1][1].values

    if 'ET_04' in exp:
        df = xr.open_dataarray(files[i])
        if '_LAI_' in exp:
            et_lai_R = df[1][0].values
        if '_LAISSM_' in exp:
            et_laissm_R = df[1][0].values
        if '_VOD10_' in exp:
            et_vod10_R = df[1][0].values
        if '_OL_' in exp:
            et_ol_R = df[1][0].values
        if '_VOD_' in exp:
            et_vod_R = df[1][1].values
        if '_VODSSM_' in exp:
            et_vodssm_R = df[1][1].values
        if '_SSM_' in exp:
            et_ssm_R = df[1][1].values

    if 'SSM_COMBINED' in exp:
        df = xr.open_dataarray(files[i])
        if '_LAI_' in exp:
            ssm_lai_R = df[1][0].values
        if '_LAISSM_' in exp:
            ssm_laissm_R = df[1][0].values
        if '_VOD10_' in exp:
            ssm_vod10_R = df[1][0].values
        if '_OL_' in exp:
            ssm_ol_R = df[1][0].values
        if '_VOD_' in exp:
            ssm_vod_R = df[1][1].values
        if '_VODSSM_' in exp:
            ssm_vodssm_R = df[1][1].values
        if '_SSM_' in exp:
            ssm_ssm_R = df[1][1].values

### Graphing/Plotting 

fig, axes = plt.subplots(2,2)
#fig = plt.figure(figsize=(3,2),dpi=150,edgecolor='w')
#fig = plt.figure()
#axes = fig.add_subplot(111)

#plt.title('Monthly Correlations per Assimilation scheme')

x = 0,1,2,3,4,5,6,7,8,9,10,11
ticks = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

### LAI
axes[0][0].set_title('LAI vs CGLS LAI')
axes[0][0].set_xticks(np.arange(12))
axes[0][0].set_xticklabels(ticks)
axes[0][0].set_xlim(-0.5,11.5)
axes[0][0].grid(True)
axes[0][0].set_ylabel('Correlation')

axes[0][0].plot(lai_ol_R,color='b',linewidth=1,label='OL')

axes[0][0].plot(lai_lai_R,color='g',linewidth=1,label='EKF LAI')
axes[0][0].plot(lai_vod_R,color='r',linewidth=1,label='EKF VOD')
axes[0][0].plot(lai_vod10_R,color='maroon',linewidth=1,label='EKF VOD10')

#axes[0][0].plot(lai_ssm_R,color='c',linewidth=1,label='EKF SSM')
#axes[0][0].plot(lai_laissm_R,color='g',linewidth=1,label='EKF LAI+SSM')
#axes[0][0].plot(lai_vodssm_R,color='r',linewidth=1,label='EKF VOD+SSM')

### GPP
axes[0][1].set_title('GPP vs FLUXCOM GPP')
axes[0][1].set_xticks(np.arange(12))
axes[0][1].set_xticklabels(ticks)
axes[0][1].set_xlim(-0.5,11.5)
axes[0][1].grid(True)
axes[0][1].set_ylabel('Correlation')

axes[0][1].plot(gpp_ol_R,color='b',linewidth=1,label='OL')

axes[0][1].plot(gpp_lai_R,color='g',linewidth=1,label='EKF LAI')
axes[0][1].plot(gpp_vod_R,color='r',linewidth=1,label='EKF VOD')
axes[0][1].plot(gpp_vod10_R,color='maroon',linewidth=1,label='EKF VOD10')

#axes[0][1].plot(gpp_ssm_R,color='c',linewidth=1,label='EKF SSM')
#axes[0][1].plot(gpp_laissm_R,color='g',linewidth=1,label='EKF LAI+SSM')
#axes[0][1].plot(gpp_vodssm_R,color='r',linewidth=1,label='EKF VOD+SSM')

### ET
axes[1][0].set_title('ET vs ALEXI ET')
axes[1][0].set_xticks(np.arange(12))
axes[1][0].set_xticklabels(ticks)
axes[1][0].set_xlim(-0.5,11.5)
axes[1][0].grid(True)
axes[1][0].set_ylabel('Correlation')

axes[1][0].plot(et_ol_R,color='b',linewidth=1,label='OL')

axes[1][0].plot(et_lai_R,color='g',linewidth=1,label='EKF LAI')
axes[1][0].plot(et_vod_R,color='r',linewidth=1,label='EKF VOD')
axes[1][0].plot(et_vod10_R,color='maroon',linewidth=1,label='EKF VOD10')

#axes[1][0].plot(et_ssm_R,color='c',linewidth=1,label='EKF SSM')
#axes[1][0].plot(et_laissm_R,color='g',linewidth=1,label='EKF LAI+SSM')
#axes[1][0].plot(et_vodssm_R,color='r',linewidth=1,label='EKF VOD+SSM')

### SSM
axes[1][1].set_title('SSM vs ESA SSM')
axes[1][1].set_xticks(np.arange(12))
axes[1][1].set_xticklabels(ticks)
axes[1][1].set_xlim(-0.5,11.5)
axes[1][1].grid(True)
axes[1][1].set_ylabel('Correlation')

axes[1][1].plot(ssm_ol_R,color='b',linewidth=1,label='OL')

axes[1][1].plot(ssm_lai_R,color='g',linewidth=1,label='EKF LAI')
axes[1][1].plot(ssm_vod_R,color='r',linewidth=1,label='EKF VOD')
axes[1][1].plot(ssm_vod10_R,color='maroon',linewidth=1,label='EKF VOD10')

#axes[1][1].plot(ssm_ssm_R,color='c',linewidth=1,label='EKF SSM')
#axes[1][1].plot(ssm_laissm_R,color='g',linewidth=1,label='EKF LAI+SSM')
#axes[1][1].plot(ssm_vodssm_R,color='r',linewidth=1,label='EKF VOD+SSM')

fig.tight_layout()
plt.legend(framealpha=1,bbox_to_anchor = (-0.10,1))
#plt.figlegend()
#fig.legend()
plt.show()
