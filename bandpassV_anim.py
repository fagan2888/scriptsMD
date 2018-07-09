#!/usr/bin/env python3

#shows the bandpass V data from a dataset, centered on India-ish
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.basemap import Basemap

#data collection
rootgroup = nc.Dataset("850mb_oct04-oct16_2014.nc", "r", format="NETCDF4")

lat = rootgroup['g4_lat_1']
ymin = 0
ymax = len(lat) -  1
while lat[ymin] > 40:
    ymin += 1
while lat[ymax] < -20:
    ymax -= 1
lat = lat[ymin:ymax]
lat.sort()

lon = rootgroup['g4_lon_2']
xmin = 0
xmax = len(lon) -  1
while lon[xmin] < 60:
    xmin += 1
while lon[xmax] > 120:
    xmax -= 1
lon = lon[xmin:xmax]

time = rootgroup['time']

v = rootgroup['v']
v = v[:, ymin:ymax, xmin:xmax]
v = np.flip(v, 1)

t = 0
#plot creation
f, ax = plt.subplots()
#draw map
m = Basemap(llcrnrlon=lon[0],llcrnrlat=lat[0],urcrnrlon=lon[-1],urcrnrlat=lat[-1],
             resolution='c', projection='merc', lat_0 = lat[int(len(lat)/2)], lon_0 = lon[int(len(lon)/2)])
#m = Basemap(resolution='c')
m.drawcoastlines()
parallels = np.arange(-90,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
meridians = np.arange(0.,360.,10.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
im = m.imshow(v[t, :, :])
#ny = v[t, :, :].shape[0]; nx = v[t, :, :].shape[1]
#lons, lats = m.makegrid(nx, ny)
#x, y = m(lons, lats)
#levels = np.linspace(np.amin(v[t, :, :]), np.amax(v[t, :, :]), 10)
#cf = m.contourf(x, y, v[t, :, :], levels)

def update(i):
    """ for animation """
    global t
    t += 1 
    t = t % len(time)
    #cf = m.contourf(x, y, v[t, :, :], levels)
    im = m.imshow(v[t, :, :])
    return im,

dt = 100
ani = FuncAnimation(plt.gcf(), update, interval=dt)
plt.show()
