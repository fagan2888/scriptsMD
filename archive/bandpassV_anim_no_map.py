#!/usr/bin/env python3

#shows the bandpass V data from a dataset, centered on India-ish
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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

t = 0

#plot creation
f, ax = plt.subplots()
im = ax.imshow(v[t, :, :])

xticks = range(0, len(lon), 10)
yticks = range(0, len(lat), 10)
ax.set_xticks(xticks)
ax.set_yticks(yticks)
xlabels = [int(lon[i]) for i in xticks]
ylabels = [int(lat[i]) for i in yticks]
ax.set_xticklabels(xlabels)
ax.set_yticklabels(ylabels)
ax
ax.set_xlabel("deg lon")
ax.set_ylabel("deg lat")

def update(i):
    """ for animation """
    global t
    t += 1 
    t = t % len(time)
    im.set_array(v[t, :, :])
    return im,

dt = 50
ani = FuncAnimation(f, update, interval=dt, blit=True)
plt.show()
