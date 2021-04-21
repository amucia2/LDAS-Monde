# Generate figures related to patch fractions
# author: B. Bonan
# Last modified: 15 July 2020

import matplotlib as mpl
#mpl.use('Agg') # case when no graphic interface is available

from ldasPost import *
from ldasMapSet import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Load options files
my_options_file = 'my_options_Afrique_Ouest.py'
initPost(my_options_file)

##################################################################################

# Reading patch fractions
PREPdir = '/cnrm/vegeo/bonanb/LDAS/LDAS_curr/Afrique_Ouest/pgd_prep/'
file_prep = 'PREP.nc'
PREP_set = Dataset(PREPdir+file_prep,'r')
Auxi = PREP_set.variables['FRAC_NATURE'][:]
FRAC_NATURE = Auxi
Auxi = PREP_set.variables['PATCH'][:]
patch_frac = Auxi.data
patch_frac[patch_frac == Auxi.fill_value] = np.nan

PREP_set.close()

dim_domain = Auxi.shape[1]*Auxi.shape[2]

# patch_domin 
patch_test = patch_frac.reshape((12,dim_domain))
patch_test2 = patch_frac.reshape((12,dim_domain))

patch_test2[np.isnan(patch_test2)] = -1

patch_max = np.argmax(patch_test2,axis=0)*1.0
patch_max2 = np.max(patch_test2,axis=0)
patch_max[patch_max2 == -1.0] = np.nan

patch_domin = pd.DataFrame(index=range(0,1),columns=range(0,dim_domain))
patch_domin.values[0,:] = patch_max

# for ip in range(0,12):

    # patch_test = patch_frac[ip,:,:].reshape((dim_domain)) 
    # patch_domin.values[:,patch_test > 0.5] = ip

grid_spec, fig, title_rel_h, cbar_aspect = map_set.figureLayout(title=None,n_rows=1,n_cols=1)

ax = fig.add_subplot(1,1,1)

color_ticks = np.linspace(0,11,12)
colorlist = [(0.7,0.4,0.0),(0.5,0.5,0.5),(0.8,0.8,0.8),(116.0/256.0,195.0/256.0,101.0/256.0),(79.0/256.0,121.0/256.0,66.0/256.0),(0.0,1.0,0.0),(1.0,0.8,0.4),(1.0,1.0,0.0),(1.0,0.3,0.3),(144.0/256.0,252.0/256.0,144.0/256.0),(1.0,0.6,0.6),(0.6,0.6,1.0)]
customCmap = mpl.colors.ListedColormap(colorlist)

pc = map_set.drawMap(patch_domin.values,color_ticks,customCmap,ax,to_draw='countries',line_widths=[1.5,1.5])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.15)

nn = 12
clb_var = fig.colorbar(pc,cax=cax,orientation='vertical', ticks=np.linspace((nn-1)/(2.0*nn),nn-1/2.0+1/(2.0*nn)-1,nn), format='%.2g')
clb_var.ax.set_yticklabels(['Bare soil','Bare rock (City)','Snow','Deciduous','Coniferous','Evergreen','C3 crops','C4 crops','C4 irrigated crops','Grasslands','Tropical','Wetlands'],rotation=0)

fig_name_tmp = 'domain_patch_Afrique_Ouest.png'
fig.savefig(fig_name_tmp, dpi=DEF_FIG_DPI, bbox_inches='tight')

pl.close()

##################################################################################

## Reading patch fractions
#PREPdir = '/cnrm/vegeo/muciaa/ldas_chain_python/ldas_curr/US_FC15/sfx-trip/pgd_prep/'
#file_prep = 'PREP.nc'
#PREP_set = Dataset(PREPdir+file_prep,'r')
#Auxi = PREP_set.variables['FRAC_NATURE'][:]
#FRAC_NATURE = Auxi
#Auxi = PREP_set.variables['PATCH'][:]
#patch_frac = Auxi.data
#patch_frac[patch_frac == Auxi.fill_value] = np.nan

#PREP_set.close()

#dim_domain = Auxi.shape[1]*Auxi.shape[2]

#patch_prairie = pd.DataFrame(index=range(0,1),columns=range(0,dim_domain))
#patch_prairie.ix[0,:] = patch_frac[0:1,:,:].reshape((dim_domain))

##patch_prairie_aux = pd.DataFrame(index=range(0,3),columns=range(0,dim_domain))
##patch_prairie_aux.ix[0:3,:] = patch_frac[0:3,:,:].reshape((3,dim_domain))
##patch_prairie = patch_prairie_aux.sum(axis=0)
##patch_prairie[patch_prairie.values == 0.0] = np.nan

#color_ticks, color_scale = discreteColorscale(patch_prairie, 'Paired', None, 13, None, False, 0.0, 0.2)

#grid_spec, fig, title_rel_h, cbar_aspect = map_set.figureLayout(title=None,n_rows=1,n_cols=2)

#ax = fig.add_subplot(1,2,1)

#pc = map_set.drawMap(patch_prairie.ix[0,:].values,color_ticks,color_scale,ax,to_draw='countries',line_widths=[.4,.4])
##pc = map_set.drawMap(patch_prairie.values,color_ticks,color_scale,ax,to_draw='countries',line_widths=[.4,.4])
#ax.set_title('Fraction patch bare soil',fontsize=16)

#cax = fig.add_subplot(grid_spec[1,0:1], aspect=cbar_aspect)
#clb_var = fig.colorbar(pc,cax=cax,orientation='horizontal', ticks=color_ticks, format='%.2g')
#clb_var.ax.tick_params(labelsize=14)
#clb_var.set_label('[-], if white pixel no bare soil',fontsize=16)

#pl.tight_layout()

#fig_name_tmp = 'fraction_patch_bare_soil_CONUS__0_02_ECOII.png'
#fig.savefig(fig_name_tmp, dpi=DEF_FIG_DPI, bbox_inches='tight')

#pl.close()
