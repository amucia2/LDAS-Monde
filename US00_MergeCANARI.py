import sys, glob, os, re, math, time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
import pylab as pl
import scipy.optimize as optim
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from math import radians, cos, sin, asin, sqrt
from scipy.stats.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    return parts


rep = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/CANARI/ssmv2/'
#files = sorted(glob.glob(rep+'*.DAT'),key=numericalSort)
for i in sorted(glob.glob(rep+'*.DAT'),key=numericalSort):
    print(i)
    name=i[-27:]
    
    #ssm = np.genfromtxt('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/CANARI/ssmv2/'+name,dtype='str')
    #vod = np.genfromtxt('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/CANARI/vodx/'+name,dtype='str')
    ssm = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/CANARI/ssmv2/'+name
    vod = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/CANARI/vodx/'+name
    output = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US00/observations/CANARI/vodx_ssm/'+name
    
    #os.system('paste -d" "' + ssm + ' ' + vod + ' > ' + output)
    
    with open(output, 'w') as file3:
        with open(ssm, 'r') as file1:
            with open(vod, 'r') as file2:
                for line1, line2 in zip(file1, file2):
                    print(line1.strip(), line2.strip(), file=file3)
                    
