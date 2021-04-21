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

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts

#exp = ('ol','ekf_lai','ekf_ssm','ekf_vod','ekf_lai_ssm','ekf_vod_ssm')
exp = ('ekf_vod10')

for i in range(len(exp)):
    
    l4 = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/{0}/WG4*'.format(exp[i])),key=numericalSort)
    l5 = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/{0}/WG5*'.format(exp[i])),key=numericalSort)
    
    for n in range(len(l4)):
        wg4 = pd.read_pickle(l4[n])
        wg5 = pd.read_pickle(l5[n])
        wg_20 = ((wg4*10)+(wg5*20))/30
        print("Saving 20cm layer PData")
        fname = l4[n][-28:]
        wg_20.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/{0}/WG_20{1}'.format(exp[i],fname))


print("Computation Complete")
sys.exit()


l4 = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/{0}/WG4*'.format(exp)),key=numericalSort)
l5 = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/{0}/WG5*'.format(exp)),key=numericalSort)

for n in range(len(l4)):
    wg4 = pd.read_pickle(l4[n])
    wg5 = pd.read_pickle(l5[n])
    wg_20 = ((wg4*10)+(wg5*20))/30
    print("Saving 20cm layer PData")
    fname = l4[n][-28:]
    wg_20.to_pickle('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/results/{0}/WG_20{1}'.format(exp,fname))
