import sys, glob, os, re, math, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
import pylab as pl
import scipy.optimize as optim
from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr
#from mpl_toolkits.axes_grid1 import make_axes_locatable

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts

doms=[]
directs=[]
dobs=[]
###################################################################
### For use on Home computer [66:0]
###################################################################
#rootdir = r'C:/Users/antho/Documents/SXVGO1/ldas_curr/US00/observations/Zones'
#imdir = r'C:/Users/antho/Documents/SXVGO1/Images/US00/VOD_Comp_V3/'

###################################################################
### For use on Work computer [71:]
### For V2 use [74:] and [:71]
###################################################################
rootdir = r'/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/Zones/V2'
imdir = r'/cnrm/vegeo/muciaa/Images/US00/VOD_Comparisons/V4/'

###################################################################
### Reading in all Data / Looping over
###################################################################

for root, dirs, files in os.walk(rootdir,topdown=False):
    for name in dirs:
        ### Change the numbers in brackets depending on which computer
        doms.append(os.path.join(root,name)[74:])
        directs.append(os.path.join(root,name))
        dobs.append(os.path.join(root,name)[:71])
        print(os.path.join(root,name)[74:])

for i in range(len(directs)):
    print("Working on Domain {0}".format(doms[i]))
    lai_isba = pd.read_pickle(directs[i] + '/LAI_ISBA.Pdata')
    vodx = pd.read_pickle(directs[i] + '/VODX.Pdata')
    vodc = pd.read_pickle(directs[i] + '/VODC.Pdata')
    lai_cgls = pd.read_pickle(directs[i] + '/LAI_CGLS.Pdata')
    
    ### Take rolling mean
    #vodx = vodx.rolling(window=30,min_periods=5,center=True).mean()
    #vodc = vodc.rolling(window=30,min_periods=5,center=True).mean()

    lai_cgls_12 = lai_cgls['2012-01-01':'2012-12-31']
    vodx_12 = vodx['2012-01-01':'2012-12-31']
    vodc_12 = vodc['2012-01-01':'2012-12-31']

    lai_cgls_10 = lai_cgls['2010-01-01':'2010-12-31']
    vodx_10 = vodx['2010-01-01':'2010-12-31']
    vodc_10 = vodc['2010-01-01':'2010-12-31']

    lai_cgls_09 = lai_cgls['2009-01-01':'2009-12-31']
    vodx_09 = vodx['2009-01-01':'2009-12-31']
    vodc_09 = vodc['2009-01-01':'2009-12-31']

###################################################################
### All Years Time Series
###################################################################
    plt.rcParams['figure.figsize'] = (15.5,5.5)
    #plt.rcParams['figure.figsize'] = (13.5,4.5)
    '''
    fig, ax1 = plt.subplots()

    plt.title('LAI vs VOD over {0} Region'.format(doms[i]))
    color='tab:green'
    ax1.set_ylabel('LAI [m^2/m^-2]',color=color)
    ax1.set_xlabel('Date',color='black')
    ax1.plot(lai_cgls.mean(axis=1),color='black',label='LAI CGLS')

    ax2 = ax1.twinx()
    color='tab:red'
    ax2.set_ylabel('VOD',color=color)
    ax2.plot(vodx.mean(axis=1),color='green',label='VODX')
    ax2.plot(vodc.mean(axis=1),color='red',label='VODC')

    fig.legend(bbox_to_anchor = (0.955,0.93))
    fig.tight_layout()
    #plt.show()
    plt.savefig(imdir+'{0}/TS_{1}_LAI_VOD_Obs.png'.format(doms[i],doms[i]),format='png',dpi=300)

###################################################################
### One year Time Series
###################################################################
    fig, ax1 = plt.subplots()

    plt.title('LAI vs VOD over {0} Region - 2012'.format(doms[i]))
    color='tab:green'
    ax1.set_ylabel('LAI [m^2/m^-2]',color=color)
    ax1.set_xlabel('Date',color='black')
    ax1.plot(lai_cgls_12.mean(axis=1),color='black',label='LAI CGLS')

    ax2 = ax1.twinx()
    color='tab:red'
    ax2.set_ylabel('VOD',color=color)
    ax2.plot(vodx_12.mean(axis=1),color='green',label='VODX')
    ax2.plot(vodc_12.mean(axis=1),color='red',label='VODC')
    fig.legend(bbox_to_anchor = (0.955,0.93))
    fig.tight_layout()
    plt.savefig(imdir+'{0}/TS_{1}_LAI_VOD_Obs_2012.png'.format(doms[i],doms[i]),format='png',dpi=300)

    fig, ax1 = plt.subplots()
    plt.title('LAI vs VOD over {0} Region - 2009'.format(doms[i]))
    color='tab:green'
    ax1.set_ylabel('LAI [m^2/m^-2]',color=color)
    ax1.set_xlabel('Date',color='black')
    ax1.plot(lai_cgls_09.mean(axis=1),color='black',label='LAI CGLS')

    ax2 = ax1.twinx()
    color='tab:red'
    ax2.set_ylabel('VOD',color=color)
    ax2.plot(vodx_09.mean(axis=1),color='green',label='VODX')
    ax2.plot(vodc_09.mean(axis=1),color='red',label='VODC')
    fig.legend(bbox_to_anchor = (0.955,0.93))
    fig.tight_layout()
    plt.savefig(imdir+'{0}/TS_{1}_LAI_VOD_Obs_2009.png'.format(doms[i],doms[i]),format='png',dpi=300)

    '''
###################################################################
### Scatter Plots
###################################################################
    
    idx = pd.DatetimeIndex(vodx.index)
    lai_cgls = lai_cgls.reindex(idx, fill_value=np.nan)
    lai_isba = lai_isba.reindex(idx, fill_value=np.nan)
    
    fig, axs = plt.subplots(1,2)
    fig.suptitle('{0}'.format(doms[i]),fontsize=16)
    axs[0].set_title('LAI CGLS vs VODX')
    axs[0].scatter(lai_cgls.mean(axis=1),vodx.mean(axis=1),c='black')
    axs[0].set_xlabel('LAI [m^2/m^-2]',color='green')
    axs[0].set_ylabel('VODX',color='red')

    axs[1].set_title('LAI ISBA vs VODX')
    axs[1].scatter(lai_isba.mean(axis=1),vodx.mean(axis=1),c='black')
    axs[1].set_xlabel('LAI [m^2/m^-2]',color='green')
    axs[1].set_ylabel('VODX',color='red')
    #plt.show()
    plt.savefig(imdir+'{0}/VODX_Scatter_{1}.png'.format(doms[i],doms[i]),format='png',dpi=300)


    idx = pd.DatetimeIndex(lai_cgls.index)
    vodc = vodc.reindex(idx, fill_value=np.nan)

    fig, axs = plt.subplots(1,2)
    fig.suptitle('{0}'.format(doms[i]),fontsize=16)
    axs[0].set_title('LAI CGLS vs VODC')
    axs[0].scatter(lai_cgls.mean(axis=1),vodc.mean(axis=1),c='black')
    axs[0].set_xlabel('LAI [m^2/m^-2]',color='green')
    axs[0].set_ylabel('VODC',color='red')

    axs[1].set_title('LAI ISBA vs VODC')
    axs[1].scatter(lai_isba.mean(axis=1),vodc.mean(axis=1),c='black')
    axs[1].set_xlabel('LAI [m^2/m^-2]',color='green')
    axs[1].set_ylabel('VODC',color='red')
    #plt.show()
    plt.savefig(imdir+'{0}/VODC_Scatter_{1}.png'.format(doms[i],doms[i]),format='png',dpi=300)
    
###################################################################
### Scatter Plots, with seasons
###################################################################
    '''
    idx = pd.DatetimeIndex(vodx.index)
    lai_cgls = lai_cgls.reindex(idx, fill_value=np.nan)
    lai_isba = lai_isba.reindex(idx, fill_value=np.nan)
    
    lai_cgls_grow = lai_cgls[~lai_cgls.index.month.isin([1,2,3,10,11,12])]
    lai_isba_grow = lai_isba[~lai_isba.index.month.isin([1,2,3,10,11,12])]
    vodx_grow = vodx[~vodx.index.month.isin([1,2,3,10,11,12])]
    vodc_grow = vodc[~vodc.index.month.isin([1,2,3,10,11,12])]

    fig, axs = plt.subplots(1,2)
    fig.suptitle('{0}'.format(doms[i]),fontsize=16)
    axs[0].set_title('LAI CGLS vs VODX')
    axs[0].scatter(lai_cgls.mean(axis=1),vodx.mean(axis=1),c='black')
    axs[0].scatter(lai_cgls_grow.mean(axis=1),vodx_grow.mean(axis=1),c='red')
    axs[0].set_xlabel('LAI [m^2/m^-2]',color='green')
    axs[0].set_ylabel('VODX',color='red')

    axs[1].set_title('LAI ISBA vs VODX')
    axs[1].scatter(lai_isba.mean(axis=1),vodx.mean(axis=1),c='black')
    axs[1].scatter(lai_isba_grow.mean(axis=1),vodx_grow.mean(axis=1),c='red')
    axs[1].set_xlabel('LAI [m^2/m^-2]',color='green')
    axs[1].set_ylabel('VODX',color='red')
    #plt.show()
    plt.savefig(imdir+'{0}/VODX_Scatter_Grow_{1}.png'.format(doms[i],doms[i]),format='png',dpi=300)
    
    idx = pd.DatetimeIndex(lai_cgls.index)
    vodc = vodc.reindex(idx, fill_value=np.nan)

    fig, axs = plt.subplots(1,2)
    fig.suptitle('{0}'.format(doms[i]),fontsize=16)
    axs[0].set_title('LAI CGLS vs VODC')
    axs[0].scatter(lai_cgls.mean(axis=1),vodc.mean(axis=1),c='black')
    axs[0].scatter(lai_cgls_grow.mean(axis=1),vodc_grow.mean(axis=1),c='red')
    axs[0].set_xlabel('LAI [m^2/m^-2]',color='green')
    axs[0].set_ylabel('VODC',color='red')

    axs[1].set_title('LAI ISBA vs VODC')
    axs[1].scatter(lai_isba.mean(axis=1),vodc.mean(axis=1),c='black')
    axs[1].scatter(lai_isba_grow.mean(axis=1),vodc_grow.mean(axis=1),c='red')
    axs[1].set_xlabel('LAI [m^2/m^-2]',color='green')
    axs[1].set_ylabel('VODC',color='red')
    #plt.show()
    plt.savefig(imdir+'{0}/VODC_Scatter_Grow_{1}.png'.format(doms[i],doms[i]),format='png',dpi=300)
    '''
