#!/usr/bin/env python3

################################################################################
# Gather ERA data following a chosen Monsoon Depression 
#   in the Bay of Bengal
################################################################################

import numpy as np
import pandas as pd
import xarray as xr
import pygrib
import datetime as dt
import os

print('\nCollecting Track ERA Data\n')

scratch_dir = '/global/scratch/hpeter/'
data_dir = scratch_dir + 'data/ERA_Interim/to_netCDF/'
harm_dir = scratch_dir + 'data/ERA_Interim/all_years/'
follow_dir = scratch_dir + 'track_following/'
track_dir = scratch_dir + 'public_trackdata/'

'''
ERA-int Variables:
uv:
    U_GDS4_ISBL
    V_GDS4_ISBL
sc:
    PV_GDS4_ISBL
    Z_GDS4_ISBL
    T_GDS4_ISBL
    Q_GDS4_ISBL
    W_GDS4_ISBL
    VO_GDS4_ISBL
    D_GDS4_ISBL
    R_GDS4_ISBL
    O3_GDS4_ISBL
    CLWC_GDS4_ISBL
    CIWC_GDS4_ISBL
    CC_GDS4_ISBL
'''
# get variable as input
variable_name = 'UV'
if variable_name == 'UV':
    grb_names = ['U component of wind', 'V component of wind']
    data_type = 'uv'
else:
    os.sys.exit('Error: Variable name not supported. Exiting now.')

print('Variable name: {}\nTranslated to: \n\t{}.\n'.format(variable_name, grb_names))

# Load event data
print('Loading event data...')
mode = 'biweekly'
phases = '1234'
date = '19940706'
trx_data_fname = '{}BoB_followed_track_{}_{}{}.npz'.format(track_dir, date, mode, phases)
trx_data = np.load(trx_data_fname)
trx_data = trx_data['arr_0']
print('Data loaded from {}'.format(trx_data_fname))

# Get time coordinate
lon, lat, year, month, day, hour, intensity = trx_data[0, :]
first_date = dt.datetime(int(year), int(month), int(day), int(hour))
lon, lat, year, month, day, hour, intensity = trx_data[-1, :]
last_date = dt.datetime(int(year), int(month), int(day), int(hour))
track_dates = pd.date_range(start=first_date, end=last_date, freq='6H')

L = len(trx_data)

# Pressure levels
levels = [300, 500, 850]

# Counting
nLevels = len(levels)
nNames = len(grb_names)

# Put all the data in a giant array for processing
data_total = np.zeros( (L, nNames, nLevels, 256, 512) )
print('data_total size: {}'.format(data_total.shape))

# Get Harmonics data to use for calculating anomaly
print('Loading seasonal harmonics data...')
harms = {}
for grb_name in grb_names:
    harms[grb_name] = {}
    for l in levels:
        harms[grb_name][str(l)] = np.load('{}ei.oper.an.pl.regn128uv.{}_{}mb_all_years_harmonics.npz'.format(harm_dir, grb_name[0], l))['arr_0']

# Begin looping through events
print('Looping through {} events.'.format(L))
for I in range(L):
    # Read genesis point data
    event_day = track_dates[I]
    event_day_str = event_day.strftime('%Y%m%d%H')
    # Get number of data points for the year
    dates_b4_event = pd.date_range(start=event_day.strftime('%Y')+"-1-1", end=event_day, freq='6H')
    nDatesB4Event = len(dates_b4_event)

    # Find the corresponding ERA data file
    idir = '/global/scratch/williamb/data/ERA_Interim/{}/'.format(event_day.year)
    ifile_name = 'ei.oper.an.pl.regn128{}.{}'.format(data_type, event_day_str)
    ifile = idir + ifile_name

    # Convert to netCDF if haven't before
    print('Reading {}'.format(ifile))

    # Read the ERA data and pick the variable we want
    grbs = pygrib.open(ifile)
    grb = grbs.select(name=grb_names, level=levels)

    # Transfer data to main array
    for K in range(nNames * nLevels):
        i = K % nNames
        j = K % nLevels
        data_total[I, i, j, :, :] = grb[K].values

    # Find mean from this day and subtract it
    time_index = nDatesB4Event - 1
    for i in range(nNames):
        for j in range(nLevels):
            data_total[I, i, j, :, :] -= harms[grb_names[i]][str(levels[j])][time_index, :, :]

    print('{} events loaded.\n'.format(I+1))

print('\nDone. Saving data.')

# Get Lats/Lons from last grib
lats, lons = grb[0].latlons()
lats = lats[:, 0]
lons = lons[0, :]

u_da = xr.DataArray(data_total[:, 0, :, :, :], 
                    coords=[track_dates, levels, lats, lons],
                    dims=['time', 'lvl', 'lat', 'lon'])
v_da = xr.DataArray(data_total[:, 1, :, :, :], 
                    coords=[track_dates, levels, lats, lons],
                    dims=['time', 'lvl', 'lat', 'lon'])
comp_ds = xr.Dataset({'u': u_da, 'v': v_da},
                    coords={'time': track_dates, 
                            'lvl': levels,
                            'lat': lats,
                            'lon': lons})
# Save for plotting
fname = '{}track_{}_{}_{}{}.nc'.format(follow_dir, variable_name, track_dates[0].strftime('%Y%m%d%H'), mode, phases)
comp_ds.to_netcdf(fname)
print('Data saved in {}.'.format(fname))

print('\nJob complete, terminating program.')
