NIC = 100*(np.array(R_ekf)-np.array(R_mod))/(1-np.array(R_mod))



#NIC = 100*(np.array(analysis_vod2)-np.array(model2))/(1-np.array(model2))
w1=np.where(NIC>3)[0]
w3=np.where((NIC>-3) & (NIC<3))[0]
w2=np.where(NIC<-3)[0]

w1=np.where(NIC2>3)[0]
w3=np.where((NIC2>-3) & (NIC2<3))[0]
w2=np.where(NIC2<-3)[0]

from matplotlib import ticker
plt.figure(figsize=(8,6),dpi=300,edgecolor='w')
plt.title("NIC R ERA5 vs R ERALand")
map = Basemap(projection='aea',llcrnrlon=-122, llcrnrlat=22, urcrnrlon=-62, urcrnrlat=48,lat_0=40,lon_0=-98,resolution='l')
map.shadedrelief()
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.8)
map.drawcountries(1.8)
#map.drawmapscale(lon=-90,lat=27,lon0=-98,lat0=40,length=750,barstyle='fancy')
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
#x,y = map(np.array(lon_keep)[w3], np.array(lat_keep)[w3])
#cs2 = map.scatter(x,y,marker='D',s=10,c=NIC[w3],vmin=-20,vmax=20,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w3], np.array(lat_keep)[w3])
cs2 = map.scatter(x,y,marker='v',s=25,c=NIC[w3],vmin=-15,vmax=15,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w2], np.array(lat_keep)[w2])
cs2 = map.scatter(x,y,marker='o',s=45,c=NIC[w2],vmin=-15,vmax=15,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w1], np.array(lat_keep)[w1])
cs2 = map.scatter(x,y,marker='o',s=45,c=NIC[w1],vmin=-15,vmax=15,cmap=cm,zorder=2)
cbar = map.colorbar(cs2,location='bottom',pad="10%")
cbar.set_label("NIC R")
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
plt.rcParams.update({'font.size': 12})
#plt.savefig('images/SMN_NIC_NoahMP_monde.png',format='png',dpi=300)
plt.show()

from matplotlib import ticker
plt.figure(figsize=(8,6),dpi=300,edgecolor='w')
plt.title("NIC R ERA5 vs R ERALand")
map = Basemap(projection='aea',llcrnrlon=-122, llcrnrlat=22, urcrnrlon=-62, urcrnrlat=48,lat_0=40,lon_0=-98,resolution='l')
map.shadedrelief()
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.8)
map.drawcountries(1.8)
#map.drawmapscale(lon=-90,lat=27,lon0=-98,lat0=40,length=750,barstyle='fancy')
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
#x,y = map(np.array(lon_keep)[w3], np.array(lat_keep)[w3])
#cs2 = map.scatter(x,y,marker='D',s=10,c=NIC[w3],vmin=-20,vmax=20,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w3], np.array(lat_keep)[w3])
cs2 = map.scatter(x,y,marker='v',s=25,c=NIC2[w3],vmin=-15,vmax=15,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w2], np.array(lat_keep)[w2])
cs2 = map.scatter(x,y,marker='o',s=45,c=NIC2[w2],vmin=-15,vmax=15,cmap=cm,zorder=2)
x,y = map(np.array(lon_keep)[w1], np.array(lat_keep)[w1])
cs2 = map.scatter(x,y,marker='o',s=45,c=NIC2[w1],vmin=-15,vmax=15,cmap=cm,zorder=2)
cbar = map.colorbar(cs2,location='bottom',pad="10%")
cbar.set_label("NIC R")
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
plt.rcParams.update({'font.size': 12})
plt.savefig('images/SMN_NIC_ERA5-ERALand.png',format='png',dpi=300)
plt.show()




from matplotlib import ticker
plt.figure(figsize=(8,6),dpi=300,edgecolor='w')
plt.title("R ASCAT vs USCRN")
map = Basemap(projection='aea',llcrnrlon=-122, llcrnrlat=22, urcrnrlon=-62, urcrnrlat=48,lat_0=40,lon_0=-98,resolution='l')
map.shadedrelief()
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.8)
map.drawcountries(1.8)
#map.drawmapscale(lon=-90,lat=27,lon0=-98,lat0=40,length=750,barstyle='fancy')
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
x,y = map(np.array(lon_keep), np.array(lat_keep))
#cs2 = map.scatter(x,y,marker='D',s=10,c=np.array(R_mod)[0],vmin=-1,vmax=1,zorder=2)
cs2 = map.scatter(x,y,marker='o',s=45,c=np.array(R_mod),vmin=-0.5,vmax=1,cmap=cm)
cbar = map.colorbar(cs2,location='bottom',pad="10%")
cbar.set_label("R")
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
plt.rcParams.update({'font.size': 12})
plt.savefig('images/SMN_ASCAT_R.png',format='png',dpi=300)
plt.show()



from matplotlib import ticker
plt.figure(figsize=(8,6),dpi=300,edgecolor='w')
plt.title("NIC R ERALand vs LDAS-Monde EKF ")
map = Basemap(projection='aea',llcrnrlon=-122, llcrnrlat=22, urcrnrlon=-62, urcrnrlat=48,lat_0=40,lon_0=-98,resolution='l')
map.shadedrelief()
cm = plt.cm.get_cmap('bwr_r')
map.drawcoastlines(1.8)
map.drawcountries(1.8)
#map.drawmapscale(lon=-90,lat=27,lon0=-98,lat0=40,length=750,barstyle='fancy')
map.drawmeridians(np.arange(-120,-60,10),labels=[0,0,0,1],linewidth=.5)
map.drawparallels(np.arange(25,55,5),labels=[1,0,0,0],linewidth=.5)
x,y = map(np.array(lon)[w3], np.array(lat)[w3])
cs2 = map.scatter(x,y,marker='v',s=25,c=NIC[w3],vmin=-15,vmax=15,cmap=cm,zorder=2)
x,y = map(np.array(lon)[w2], np.array(lat)[w2])
cs2 = map.scatter(x,y,marker='o',s=45,c=NIC[w2],vmin=-15,vmax=15,cmap=cm,zorder=2)
x,y = map(np.array(lon)[w1], np.array(lat)[w1])
cs2 = map.scatter(x,y,marker='o',s=45,c=NIC[w1],vmin=-15,vmax=15,cmap=cm,zorder=2)
cbar = map.colorbar(cs2,location='bottom',pad="10%")
cbar.set_label("NIC R")
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
plt.rcParams.update({'font.size': 12})
plt.savefig('images/SMN_NIC_ERALand-LDASMondeEKF.png',format='png',dpi=300)
plt.show()
