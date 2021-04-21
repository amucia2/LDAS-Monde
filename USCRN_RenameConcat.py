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
import csv

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts

files = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/USCRN/original/2011/*'),key=numericalSort)
fold = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/USCRN/original/20*'),key=numericalSort)

filenames = []

for i in range(len(files)): 
    ### import wban and isolate file name
    wban = int(np.loadtxt(files[i],delimiter=' ',usecols=[0])[0])
    fname = files[i][81:]
    
    ### Search for the same station in other years
    ### Concatenating all Years of the same station and output a single file
    for y in range(len(fold)):
        nlist = sorted(glob.glob(fold[y]+'/*'),key=numericalSort)
        for f in range(len(nlist)):
            if fname in nlist[f]:
                print("Found!")
                filenames.append(nlist[f])
    with open('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/USCRN/concat/'+str(wban)+'_'+fname, 'w') as outfile:
        for fn in filenames:
            with open(fn) as infile:
                outfile.write(infile.read())
    print('File Saved')
    filenames = []
    

### Match Coords and wban
coords = np.loadtxt('/cnrm/vegeo/muciaa/NO_SAVE/US00/USCRN_coords.txt') 
files = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/USCRN/concat/*'),key=numericalSort)
a = []
b = []
c = []
count = 0
### import wban and isolate file name
for i in range(len(coords)):
    lat = coords[i][0]
    lon = coords[i][1]
    for n in range(len(files)): 
        xlat = np.loadtxt(files[n],usecols=[4])[0]
        xlon = np.loadtxt(files[n],usecols=[3])[0]
        if lat == xlat and lon == xlon:
            wban = int(np.loadtxt(files[n],usecols=[0])[0])
            a.append(xlat)
            b.append(xlon)
            c.append(wban)
            print("Match found")
            print(files[n])
            print("Lat in list: " + str(lat))
            print("Lat in file: " + str(xlat))
            count += 1
            print("Match Count = " + str(count))
            print("i = " + str(i))
            
            

with open('USCRN_coords_wban.txt', 'w') as f:
    writer = csv.writer(f, delimiter=' ')
    writer.writerows(zip(a,b,c))


