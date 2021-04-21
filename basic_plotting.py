
f = et['2017-12-17']
extent = [-125,24.8,-66,50]
f = f.values.reshape((175,350))

#map1 = Basemap(projection='cyl', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3],resolution='i')
map1 = Basemap(projection='cyl',resolution='i')
#map1.drawstates(linewidth=1.8)
#map1.drawcoastlines(1.8)
#map1.drawmapboundary()
map1.imshow(f,interpolation='none',vmin=0.,vmax=5.)
plt.show()
