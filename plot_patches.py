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
import math
import scipy.optimize as opt

FILL = 999.

def readFromFile(file_name,var_names,row_slice=None,col_slice=None,row_invert=False):
    '''
    Read values from netcdf or hdf5 file
    Arguments:
        file_name (str): full path to the file to read
        var_names (list): names of variables to read
        row_slice/col_slice (slices): indices to slice the variables to the output region
        row_invert (bool): it True, row order is inverted
    Returns:
        obs_values (list): values read from file
    '''
    # determine format
    if file_name[-2:] not in ('h5','nc'):
        raise Exception('filename extension must be h5 of nc')
    file_format = file_name[-2:]
    if isinstance(var_names,basestring): var_names = [var_names]
    obs_values = []
    if file_format == 'nc':
        fid = Dataset(file_name)
        for i in range(len(var_names)):
            v = fid.variables[var_names[i]]
            tmp = v[:].squeeze()
            if row_invert: tmp = tmp[::-1,:]
            if row_slice is not None and col_slice is not None: tmp = tmp[row_slice,col_slice]
            tmp = tmp.astype(float)
            if not type(tmp) == pl.ma.core.MaskedArray: tmp = pl.ma.masked_array(tmp)
            for a in v.ncattrs():
                if a.lower() in ['missing_value','_fillvalue']: tmp.mask = pl.logical_or(tmp.mask,(tmp==float(v.getncattr(a))))
            tmp = tmp.filled(FILL)
            obs_values.append(tmp)
        fid.close()
    return obs_values

def getModelPatchFractions(mod_pgd_path):
    '''
    Read the patch fractions of model grid points from file.
    Arguments:
        mod_pgd_path (str): directory of PREP.nc
    Returns:
        mod_patch_frac (numpy.ndarray): patch fractions of each model grid point; shape = n_points,n_patches
    '''
    prep_file = mod_pgd_path+'PREP.nc'
    if not os.path.isfile(prep_file):
        raise Exception(prep_file+' does not exist.')
    frac = readFromFile(prep_file,['PATCH'])[0]
    return frac.reshape((frac.shape[0],-1))

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

patch_frac_25 = getModelPatchFractions('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_ECOSG/sfx-trip/pgd_prep/')
patch_frac_10 = getModelPatchFractions('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_HR/sfx-trip/pgd_prep/SG/')
patch_frac_1 = getModelPatchFractions('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_VHR/sfx-trip/pgd_prep/')

patch_frac_II = getModelPatchFractions('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_LDAS/sfx-trip/pgd_prep/')
patch_frac_SG = getModelPatchFractions('/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/Nebraska_ECOSG/sfx-trip/pgd_prep/')

patch_frac_25[patch_frac_25==999.0]=0.
patch_frac_10[patch_frac_10==999.0]=0.
patch_frac_1[patch_frac_1==999.0]=0.
patch_frac_II[patch_frac_II==999.0]=0.
patch_frac_SG[patch_frac_SG==999.0]=0.

c4_25 = patch_frac_25[7]
c4_10 = patch_frac_10[7]
c4_1 = patch_frac_1[7]
c4_II = patch_frac_II[7]
c4_SG = patch_frac_SG[7]

v1 = c4_25.reshape((16,40))
v2 = c4_10.reshape(40,100)
v3 = c4_1.reshape(400,1000)
v4 = c4_II.reshape(16,40)
v5 = c4_SG.reshape(16,40)

#fig, axes = plt.subplots(3,1)

#axes[0].set_title("C4 Fraction - 25km")
plt.title("C4 Fraction - 10km")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='i')
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1),linewidth=.5)
map.drawparallels(np.array([40.5,42.5]),labels=(1,0,0,1),linewidth=.5)
map.drawrivers()
im = map.imshow(v2,interpolation='none',cmap='YlGn',vmin=0.,vmax=1.)
'''
axes[1].set_title("C4 Fraction - 10km")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = axes[1])
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map.drawrivers()
map.imshow(v2,interpolation='none',cmap='YlGn',vmin=0.,vmax=1.)

#diff = v1-v2

axes[2].set_title("C4 Fraction - 1km")
map = Basemap(projection='cyl', llcrnrlon=-105, llcrnrlat=39.5, urcrnrlon=-95, urcrnrlat=43.5,resolution='h',ax = axes[2])
map.drawstates(linewidth=1.8)
map.drawcoastlines(1.8)
map.drawmapboundary()
map.drawmeridians(np.arange(-105,-95,2),labels=(1,0,0,1))
map.drawparallels(np.arange(39.5,43.5,2),labels=(1,0,0,1))
map.drawrivers()
map.imshow(v3,interpolation='none',cmap='YlGn',vmin=0.,vmax=1.)
#map.imshow(v3,interpolation='none',cmap='seismic',vmin=-0.5,vmax=0.5)
'''


### for 'add_axes' it is [left,bottom,witdth,height], plt.figure(#) is seperate image, or on the same plot
#fig = plt.figure(0)
#cbaxes = fig.add_axes([0.92,0.15,0.025,0.7])
#fig.colorbar(im,ax=axes.ravel().tolist())
#fig.colorbar(im,cax=cbaxes,label='C4 Fraction')

#save('vegeoFigures/C4Frac10km_SG', ext='ps', close=True, verbose=True)
plt.show()


