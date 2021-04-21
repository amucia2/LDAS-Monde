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
import pickle
os.chdir('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr')
from ldasBoot import *
os.chdir('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing')

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)

    return parts


files = sorted(glob.glob('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/CustomPostProcessing/Data/US00/ldasPostScores/*scores_for_all_points*'),key=numericalSort)

for i in range(len(files)):
    
