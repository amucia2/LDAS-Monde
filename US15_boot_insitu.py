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

os.chdir('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr')
from ldasBoot import *
os.chdir('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing')
from glob import glob

def numericalSort(value):
    parts = numbers.split(value)
    return parts

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

files = glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/conus_R_US15/*0.20_R*')

files.sort(key=natural_keys)
f = []
m = []
llim = []
ulim = []

### Go through saved PData to calculate indivdual CI
for n in range(len(files)):
    print(files[n])
    dat = pd.read_pickle(files[n])
    R = dat.values
    f.append(files[n][86:-17])
    m.append(dat.mean()[0])
    llim.append(bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000)[0])
    ulim.append(bootci(R[~np.isnan(R)],stat=np.mean,nboot=10000)[1])
    print("File     : " + str(files[n]))
    print("Mean R   : %.3f"%m[n])
    print("Lower CI : %.3f"%llim[n])
    print("Upper CI : %.3f"%ulim[n])

m = np.array(m)
llim = np.array(llim)
ulim = np.array(ulim)

df = pd.DataFrame(list(zip(f,m,llim,ulim)),columns=['File','Mean','LowerLim','UpperLim'])
df.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/Stats/US15_WG20_CI.PData')

