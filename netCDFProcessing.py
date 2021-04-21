nc = Dataset('FORCING_2017010109.nc')
nc2 = Dataset('FCFORCING_2017010109.nc')

s = nc.variables['Tair'][:]
s2 = nc2.variables['Tair'][:]


nofc0 = s[0]
fc0 = s2[0].reshape(175,350)

nofc3 = s[3]
fc3 = s2[1].reshape(175,350)

nofc6 = s[6]
fc6 = s2[2].reshape(175,350)
