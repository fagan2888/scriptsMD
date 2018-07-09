#!/usr/bin/env python3

################################################################################
# This script goes through the ERA-Interim data in Boos' scratch and
#   combines data from a single variable and level into one file
#   for an entire year.
################################################################################

import pygrib
import pandas as pd
import xarray as xr
import numpy as np

idir = '/global/scratch/users/williamb/data/ERA_Interim/'
odir = '/global/scratch/users/hpeter/data/ERA_Interim/'
level = 850
name = 'Geopotential'
for year in range(1979, 2017): # 2017 only has half a year, stop at 2016
    times = pd.date_range(start=str(year)+'-1-1', end=str(year)+'-12-31', freq='6H')
    data = np.zeros((len(times), 256, 512))
    for i in range(len(times)):
        date = times[i].strftime('%Y%m%d%H')
        grbs = pygrib.open('{}{}/ei.oper.an.pl.regn128sc.{}'.format(idir, year, date))
        grb = grbs.select(name=name, level=level)[0]
        data[i, :, :] = grb.values
    lats, lons = grb.latlons()
    lats = lats[:, 0]
    lons = lons[0, :]
    da = xr.DataArray(data, coords=[times, lats, lons], dims=['time', 'lat', 'lon'])
    da.to_netcdf('{}{}/ei.oper.an.pl.regn128sc.{}_{}mb_{}.nc'.format(odir, year, name, level, year))
