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
from nco import Nco

nco = Nco()

files = glob.glob('/cnrm/vegeo/muciaa/NO_SAVE/US00/ekf_lai_swi/*ANA*.nc')

for i in range(len(files)):
    H1 = nco.ncks(input=files[i], options='-v HO1_1_1')[:20]
    H2 = nco.ncks(input=files[i], options='-v HO2_1_1')
    if H1.find('for observation LAI and control variable LAI') < 0 && H2.find('for observation LAI and control variable LAI') < 0:
        print('LAI Obs NOT Assimilated - ' + str(files[i]))
    if H2.find('for observation WG2 and control variable WG2') < 0 && H2.find('for observation WG2 and control variable WG2') < 0:
        print('SSM Obs NOT Assimilated - ' + str(files[i]))
