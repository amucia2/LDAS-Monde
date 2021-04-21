fig, axes = plt.subplots(1,1)
map = Basemap(projection='aea',llcrnrlon=-122, llcrnrlat=22, urcrnrlon=-62, urcrnrlat=48,lat_0=40,lon_0=-98,resolution='l')
map.drawcoastlines()
map.drawcountries()
m1 = map.imshow(v2,interpolation='none',cmap='coolwarm_r',vmin=0,vmax=1)
plt.show()

fig, axes = plt.subplots(1,1)
map = Basemap(llcrnrlon=-130, llcrnrlat=20, urcrnrlon=-60, urcrnrlat=55,lat_0=40,lon_0=-98,resolution='l')
map.drawcoastlines(1.5)
map.drawcountries(1.5)
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
m1 = map.imshow(v2,origin='lower',vmin=-0.5,vmax=1.)
plt.title('R : VODCA-X vs. LAI GEOV2') 
cbar = map.colorbar(m1,location='bottom')
cbar.set_label("R")
plt.show()








plt.imshow(v2,origin='lower',vmin=-0.5,vmax=1.)
plt.title('R : VODCA-X vs. LAI GEOV2')
cbar = plt.colorbar(orientation='horizontal')

