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

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

### Define Statistics
### FC15 only, every day PER POINT and Filter_NaN is ON
R_lai_ol15 = 0.554, 0.556, 0.55, 0.551, 0.553, 0.553, 0.556, 0.569, 0.569, 0.561, 0.558, 0.561, 0.556, 0.553 
R_lai_ekf15 = 0.689, 0.692, 0.687, 0.689, 0.692, 0.694, 0.696, 0.709, 0.712, 0.709, 0.658, 0.649, 0.644, 0.642
RMSD_lai_ol15 = 1.020,1.019,1.018,1.016,1.015,1.013,1.010,1.007,1.010,1.011,1.013,1.011,1.012,1.015
RMSD_lai_ekf15 =0.736,0.735,0.735,0.734,0.733,0.731,0.729,0.727,0.727,0.729,0.802,0.816,0.819,0.821

R_ssm_ol15 = 0.633, 0.62, 0.601, 0.579, 0.552, 0.522, 0.494, 0.463, 0.432, 0.408, 0.386, 0.364, 0.353, 0.346
R_ssm_ekf15 = 0.659, 0.639, 0.614, 0.588, 0.558, 0.526, 0.497, 0.465, 0.433, 0.408, 0.386, 0.363, 0.352, 0.345
RMSD_ssm_ol15 = 0.043,0.044,0.045,0.046,0.048,0.050,0.052,0.053,0.055,0.056,0.058,0.059,0.059,0.060
RMSD_ssm_ekf15 = 0.041,0.042,0.044,0.046,0.048,0.049,0.051,0.053,0.055,0.056,0.057,0.059,0.059,0.060

R_eva_ol15 = 0.734, 0.735, 0.733, 0.73, 0.724, 0.717, 0.713, 0.705, 0.7, 0.694, 0.687, 0.68, 0.673, 0.672
R_eva_ekf15 = 0.74, 0.741, 0.739, 0.735, 0.729, 0.723, 0.718, 0.711, 0.705, 0.699, 0.691, 0.685, 0.678, 0.676
RMSD_eva_ol15 = 1.040, 1.040, 1.043,1.048,1.055,1.062,1.067,1.076,1.082,1.087,1.099,1.104,1.113,1.114
RMSD_eva_ekf15 = 1.025, 1.025, 1.029,1.035,1.043,1.050,1.055,1.065,1.071,1.077,1.089,1.095,1.103,1.104

R_et_ol15 = 0.569, 0.572, 0.573, 0.571, 0.57, 0.566, 0.564, 0.558, 0.555, 0.551, 0.549, 0.545, 0.539, 0.539
R_et_ekf15 = 0.583, 0.585, 0.585, 0.583, 0.581, 0.577, 0.574, 0.569, 0.565, 0.561, 0.558, 0.554, 0.548, 0.548
RMSD_et_ol15 = 1.374, 1.368, 1.368,1.369,1.372,1.377,1.383,1.390,1.392,1.397,1.399,1.405,1.418,1.420
RMSD_et_ekf15 = 1.358, 1.353, 1.353,1.355,1.359,1.364,1.371,1.378,1.381,1.387,1.389,1.395,1.408,1.410


### ErrorBars 99% CI
### LAI
elai_ol = 0.004, 0.004, 0.004, 0.004, 0.004, 0.0035, 0.004, 0.004, 0.004, 0.004, 0.004, 0.0035, 0.0035, 0.0035
elai_ekf = 0.0035, 0.0035, 0.0035, 0.004, 0.0035, 0.0035, 0.0035, 0.0035, 0.003, 0.0035, 0.0035, 0.003, 0.003, 0.0035

### SSM
essm_ol = 0.003, 0.0025, 0.003, 0.00299999999999995, 0.0025, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.00299999999999997, 0.003
essm_ekf = 0.0025, 0.003, 0.003, 0.0025, 0.003, 0.0015, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.00299999999999997, 0.003

### EVAP
eeva_ol = 0.0015, 0.0015, 0.0015, 0.0015, 0.0015, 0.0015, 0.0015, 0.0015, 0.002, 0.0015, 0.002, 0.002, 0.002, 0.0015
eeva_ekf = 0.0015, 0.001, 0.0015, 0.002, 0.002, 0.0015, 0.002, 0.0015, 0.00199999999999995, 0.002, 0.0015, 0.002, 0.002, 0.002

eet_ol = 0.0025, 0.003, 0.0025, 0.003, 0.0025, 0.003, 0.0025, 0.003, 0.003, 0.003, 0.003, 0.0025, 0.0025, 0.0025
eet_ekf = 0.003, 0.0025, 0.003, 0.003, 0.003, 0.00299999999999995, 0.00249999999999995, 0.0025, 0.0025, 0.0025, 0.0025, 0.003, 0.0025, 0.0025


#'''

### Create Graphs ###
def indplot():
    fig = plt.figure(figsize=(8,2.5),dpi=300,edgecolor='w')
    ax1 = fig.add_subplot(111)
    x = 0,1,2,3,4,5,6,7,8,9,10,11,12,13
    ticks = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14']
    plt.xticks(np.arange(14),ticks)
    ax1.set_ylabel('R')
    ax1.set_xlabel('fc Lead Time (Days)')
    ax1.set_xlim(-0.5,13.5)
    import matplotlib.transforms as transforms
    offset = transforms.ScaledTranslation(0.025, 0,fig.dpi_scale_trans)
    trans = ax1.transData-offset
    trans2 = ax1.transData+offset
    ax2 = ax1.twinx()
    ax2.set_xlim(-0.5,13.5)
    trans_ = ax2.transData-offset
    trans2_ = ax2.transData+offset
    #plt.axvline(x=0.5,color='black',linestyle='--')
    #plt.axvline(x=5.5,color='black',linestyle='--')

ana = input("Which variable to plot (SSM, LAI, EVAP, ALL)? : ")

if ana == "EVAP":
    indplot()
    fig = plt.figure(figsize=(8,2.5),dpi=300,edgecolor='w')
    ax1 = fig.add_subplot(111)
    x = 0,1,2,3,4,5,6,7,8,9,10,11,12,13
    ticks = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14']
    plt.xticks(np.arange(14),ticks)
    ax1.set_ylabel('R')
    ax1.set_xlabel('fc Lead Time (Days)')
    ax1.set_xlim(-0.5,13.5)
    import matplotlib.transforms as transforms
    offset = transforms.ScaledTranslation(0.025, 0,fig.dpi_scale_trans)
    trans = ax1.transData-offset
    trans2 = ax1.transData+offset
    ax2 = ax1.twinx()
    ax2.set_xlim(-0.5,13.5)
    trans_ = ax2.transData-offset
    trans2_ = ax2.transData+offset
    #plt.axvline(x=0.5,color='black',linestyle='--')
    #plt.axvline(x=5.5,color='black',linestyle='--')
    plt.title('ISBA Forecast Evapotranspiration vs ALEXI Evapotranspiration',fontsize=7)
    ### R
    #plt.yticks(np.arange(0.66,0.751,step=0.03))
    ### GLEAM
    #l1 = ax1.plot(R_eva_ol15,color='b',marker='.',linewidth=0.1,label='R fc_init_ol',transform=trans,markeredgecolor='black',markersize=5)[0]
    #l2 = ax1.plot(R_eva_ekf15,color='r',marker='.',linewidth=0.1,label='R fc_init_ekf',transform=trans2,markeredgecolor='black',markersize=5)[0]
    #ax1.errorbar(x,R_eva_ol15,yerr=eeva_ol,xerr=None,transform=trans,ecolor='black',linewidth=1,color='b')
    #ax1.errorbar(x,R_eva_ekf15,yerr=eeva_ekf,xerr=None,transform=trans2,ecolor='black',linewidth=1,color='r')
    ### ALEXI
    ax1.set_ylim(0.3,0.75)
    ax1.set_yticks(np.arange(0.30,0.75,step=0.2))
    l3 = ax1.plot(R_et_ol15,color='b',marker='.',linewidth=1,label='R fc_init_ol',transform=trans,markeredgecolor='black',markersize=5)[0]
    l4 = ax1.plot(R_et_ekf15,color='r',marker='.',linewidth=1,label='R fc_init_ekf',transform=trans2,markeredgecolor='black',markersize=5)[0]
    #ax1.errorbar(x,R_et_ol15,yerr=eet_ol,xerr=None,transform=trans,ecolor='black',linewidth=1,color='b',linestyle='--')
    #ax1.errorbar(x,R_et_ekf15,yerr=eet_ekf,xerr=None,transform=trans2,ecolor='black',linewidth=1,color='r',linestyle='--')
    ### RMSD
    ax2.set_ylabel('RMSD[mm/day]')
    ax2.set_ylim(1.34,1.43)
    #plt.yticks(np.arange(1.3,1.46,step=0.05))
    l5 = ax2.plot(RMSD_et_ol15,color='b',marker='s',linewidth=1,label='RMSD fc_init_ol',transform=trans_,markeredgecolor='black',markersize=4,linestyle='--')[0]
    l6 = ax2.plot(RMSD_et_ekf15,color='r',marker='s',linewidth=1,label='RMSD fc_init_ekf',transform=trans2_,markeredgecolor='black',markersize=4,linestyle='--')[0]
    handles, labels = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    plt.legend(handles+handles2,labels+labels2,numpoints=1,fontsize=6,loc='lower right')

elif ana == "SSM":
    indplot()
    fig = plt.figure(figsize=(8,2.5),dpi=300,edgecolor='w')
    ax1 = fig.add_subplot(111)
    x = 0,1,2,3,4,5,6,7,8,9,10,11,12,13
    ticks = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14']
    plt.xticks(np.arange(14),ticks)
    ax1.set_ylabel('R')
    ax1.set_xlabel('fc Lead Time (Days)')
    ax1.set_xlim(-0.5,13.5)
    import matplotlib.transforms as transforms
    offset = transforms.ScaledTranslation(0.025, 0,fig.dpi_scale_trans)
    trans = ax1.transData-offset
    trans2 = ax1.transData+offset
    ax2 = ax1.twinx()
    ax2.set_xlim(-0.5,13.5)
    trans_ = ax2.transData-offset
    trans2_ = ax2.transData+offset
    #plt.axvline(x=0.5,color='black',linestyle='--')
    #plt.axvline(x=5.5,color='black',linestyle='--')
    plt.title('ISBA Forecast SSM vs CGLS SSM',fontsize=7)
    ### R
    ax1.set_ylim(0.30,0.75)
    ax1.set_yticks(np.arange(0.30,0.75,step=0.2))
    l1 = ax1.plot(R_ssm_ol15,color='b',marker='.',linewidth=1,label='R fc_init_ol',transform=trans,markeredgecolor='black',markersize=6)[0]
    l2 = ax1.plot(R_ssm_ekf15,color='r',marker='.',linewidth=1,label='R fc_init_ekf',transform=trans2,markeredgecolor='black',markersize=6)[0]
    #ax1.errorbar(x,R_ssm_ol15,yerr=essm_ol,xerr=None,transform=trans,ecolor='black',linewidth=1,color='b')
    #ax1.errorbar(x,R_ssm_ekf15,yerr=essm_ekf,xerr=None,transform=trans2,ecolor='black',linewidth=1,color='r')
    ### RMSD
    ax2.set_ylabel('RMSD[$m^3$/$m^3$]')
    ax2.set_ylim(0.0375,0.062)
    plt.yticks(np.arange(0.040,0.061,step=0.01))
    l5 = ax2.plot(RMSD_ssm_ol15,color='b',marker='s',linewidth=1,label='RMSD fc_init_ol',transform=trans_,markeredgecolor='black',markersize=4,linestyle='--')[0]
    l6 = ax2.plot(RMSD_ssm_ekf15,color='r',marker='s',linewidth=1,label='RMSD fc_init_ekf',transform=trans2_,markeredgecolor='black',markersize=4,linestyle='--')[0]
    handles, labels = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles+handles2,labels+labels2,numpoints=1,fontsize=6,bbox_to_anchor=(0.985,0.75))

elif ana == "LAI": 
    indplot()
    fig = plt.figure(figsize=(8,2.5),dpi=300,edgecolor='w')
    ax1 = fig.add_subplot(111)
    x = 0,1,2,3,4,5,6,7,8,9,10,11,12,13
    ticks = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14']
    plt.xticks(np.arange(14),ticks)
    ax1.set_ylabel('R')
    ax1.set_xlabel('fc Lead Time (Days)')
    ax1.set_xlim(-0.5,13.5)
    import matplotlib.transforms as transforms
    offset = transforms.ScaledTranslation(0.025, 0,fig.dpi_scale_trans)
    trans = ax1.transData-offset
    trans2 = ax1.transData+offset
    ax2 = ax1.twinx()
    ax2.set_xlim(-0.5,13.5)
    trans_ = ax2.transData-offset
    trans2_ = ax2.transData+offset
    #plt.axvline(x=0.5,color='black',linestyle='--')
    #plt.axvline(x=5.5,color='black',linestyle='--')
    plt.title('ISBA Forecast LAI vs CGLS LAI',fontsize=7)
    ### R
    ax1.set_ylim(0.30,0.75)
    ax1.set_yticks(np.arange(0.30,0.75,step=0.2))
    l1 = ax1.plot(R_lai_ol15,color='b',marker='.',linewidth=1,label='R fc_init_ol',transform=trans,markeredgecolor='black',markersize=6)[0]
    l2 = ax1.plot(R_lai_ekf15,color='r',marker='.',linewidth=1,label='R fc_init_ekf',transform=trans2,markeredgecolor='black',markersize=6)[0] 
    #ax1.errorbar(x,R_lai_ol15,yerr=elai_ol,xerr=None,transform=trans,ecolor='black',linewidth=1,color='b')
    #ax1.errorbar(x,R_lai_ekf15,yerr=elai_ekf,xerr=None,transform=trans2,ecolor='black',linewidth=1,color='r') 
    ### RMSD
    ax2.set_ylabel('RMSD[$m^2$/$m^2$]')
    ax2.set_ylim(0.5,1.20)
    plt.yticks(np.arange(0.6,1.11,step=0.2))
    l5 = ax2.plot(RMSD_lai_ol15,color='b',marker='s',linewidth=1,label='RMSD fc_init_ol',transform=trans_,markeredgecolor='black',markersize=4,linestyle='--')[0]
    l6 = ax2.plot(RMSD_lai_ekf15,color='r',marker='s',linewidth=1,label='RMSD fc_init_ekf',transform=trans2_,markeredgecolor='black',markersize=4,linestyle='--')[0]

    handles, labels = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles+handles2,labels+labels2,numpoints=1,fontsize=6,loc='lower right')

elif ana == 'ALL': 
    x = 0,1,2,3,4,5,6,7,8,9,10,11,12,13
    #x= 1,2,3,4,5,6,7,8,9,10,11,12,13,14
    ticks = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14']
    fig, axes = plt.subplots(3,1,sharex = True)
    import matplotlib.transforms as transforms
    offset = transforms.ScaledTranslation(0.025, 0,fig.dpi_scale_trans)
    trans = axes[0].transData-offset
    trans2 = axes[0].transData+offset
    ax1 = axes[0].twinx()
    trans_ = ax1.transData-offset
    trans2_ = ax1.transData+offset
    ### SSM  
    axes[0].set_ylabel('Correlation')
    axes[0].set_title('ISBA Forecast SSM vs CGLS SSM')
    axes[0].set_ylim(0.30,0.75)
    axes[0].set_yticks(np.arange(0.30,0.75,step=0.2))
    plt.yticks(np.arange(0.50,0.75,step=0.1))
    l1 = axes[0].plot(x,R_ssm_ol15,color='b',marker='.',linewidth=1,label='R fc_init_ol',transform=trans,markeredgecolor='black',markersize=10)[0]
    l2 = axes[0].plot(x,R_ssm_ekf15,color='r',marker='.',linewidth=1,label='R fc_init_sekf',transform=trans2,markeredgecolor='black',markersize=10)[0]
    ax1.set_ylabel('RMSD[$m^3$/$m^3$]')
    ax1.set_ylim(0.0375,0.062)
    ax1.set_yticks(np.arange(0.040,0.061,step=0.01))
    l5 = ax1.plot(x,RMSD_ssm_ol15,color='b',marker='s',linewidth=1,label='RMSD fc_init_ol',transform=trans_,markeredgecolor='black',markersize=6,linestyle='--')[0]
    l6 = ax1.plot(x,RMSD_ssm_ekf15,color='r',marker='s',linewidth=1,label='RMSD fc_init_sekf',transform=trans2_,markeredgecolor='black',markersize=6,linestyle='--')[0]
    ### LAI
    axes[1].set_ylabel('Correlation')
    axes[1].set_title('ISBA Forecast LAI vs CGLS LAI')
    ax2 = axes[1].twinx()
    trans_2 = axes[1].transData-offset
    trans2_2 = axes[1].transData+offset
    trans__ = ax2.transData-offset
    trans2__ = ax2.transData+offset
    axes[1].set_ylim(0.30,0.75)
    axes[1].set_yticks(np.arange(0.30,0.75,step=0.2))
    m1 = axes[1].plot(x,R_lai_ol15,color='b',marker='.',linewidth=1,label='R fc_init_ol',transform=trans_2,markeredgecolor='black',markersize=10)[0]
    m2 = axes[1].plot(x,R_lai_ekf15,color='r',marker='.',linewidth=1,label='R fc_init_ekf',transform=trans2_2,markeredgecolor='black',markersize=10)[0]
    ax2.set_ylabel('RMSD[$m^2$/$m^2$]')
    ax2.set_ylim(0.55,1.20)
    ax2.set_yticks(np.arange(0.6,1.21,step=0.2))
    m5 = ax2.plot(x,RMSD_lai_ol15,color='b',marker='s',linewidth=1,label='RMSD fc_init_ol',transform=trans__,markeredgecolor='black',markersize=6,linestyle='--')[0]
    m6 = ax2.plot(x,RMSD_lai_ekf15,color='r',marker='s',linewidth=1,label='RMSD fc_init_ekf',transform=trans2__,markeredgecolor='black',markersize=6,linestyle='--')[0]
    ### EVAP
    axes[2].set_ylabel('Correlation')
    axes[2].set_title('ISBA Forecast Evapotranspiration vs ALEXI Evapotranspiration')
    ax3 = axes[2].twinx()
    trans_3 = axes[2].transData-offset
    trans2_3 = axes[2].transData+offset
    trans___ = ax3.transData-offset
    trans2___ = ax3.transData+offset
    axes[2].set_ylim(0.3,0.75)
    axes[2].set_yticks(np.arange(0.3,0.75,step=0.2))
    n1 = axes[2].plot(x,R_et_ol15,color='b',marker='.',linewidth=1,label='R fc_init_ol',transform=trans_3,markeredgecolor='black',markersize=10)[0]
    n2 = axes[2].plot(x,R_et_ekf15,color='r',marker='.',linewidth=1,label='R fc_init_ekf',transform=trans2_3,markeredgecolor='black',markersize=10)[0]
    ax3.set_ylabel('RMSD[mm/day]')
    ax3.set_ylim(1.34,1.43)
    ax3.set_yticks(np.arange(1.34,1.43,step=0.04))
    n5 = ax3.plot(x,RMSD_et_ol15,color='b',marker='s',linewidth=1,label='RMSD fc_init_ol',transform=trans___,markeredgecolor='black',markersize=6,linestyle='--')[0]
    n6 = ax3.plot(x,RMSD_et_ekf15,color='r',marker='s',linewidth=1,label='RMSD fc_init_ekf',transform=trans2___,markeredgecolor='black',markersize=6,linestyle='--')[0]
    
    #axes[2].set_xlim(0.5,14.5)
    axes[2].set_xlim(-0.5,13.5)
    axes[2].set_xlabel('fc Lead Time (Days)')
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(ticks)
    handles, labels = axes[0].get_legend_handles_labels()
    handles2, labels2 = ax1.get_legend_handles_labels()
    ax1.legend(handles+handles2,labels+labels2,numpoints=1,fontsize=10,bbox_to_anchor=(0.985,0.75))
    axes[0].annotate("A", xy=(0.015, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')
    axes[1].annotate("B", xy=(0.015, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')
    axes[2].annotate("C", xy=(0.015, 0.025), xycoords="axes fraction",fontsize=14,weight='bold')


plt.tight_layout()
#fig.legend(fontsize=5,bbox_to_anchor=(0.9,0.65))
#fig.legend(numpoints=1,fontsize=5,bbox_to_anchor=(0.9,0.65))
plt.show()
#sys.exit()



