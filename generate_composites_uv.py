#!/usr/bin/env python3

################################################################################
# Generate composites for Monsoon Depression genesis points
# in the Bay of Bengal (data from Boos)
################################################################################

import numpy as np
import pandas as pd
import xarray as xr
import pygrib
from scipy.stats import ttest_1samp
import datetime as dt
import os

# Simple tools
def fix_num(n):
    ''' Make sure each month/day/hour has a zero on the front when needed '''
    n = int(n)
    if n < 10:
        return '0' + str(n)
    else:
        return str(n)

scratch_dir = '/global/scratch/hpeter/'
data_dir = scratch_dir + 'data/ERA_Interim/to_netCDF/'
harm_dir = scratch_dir + 'data/ERA_Interim/all_years/'
composite_dir = scratch_dir + 'composites/'
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
    variables = ['U_GDS4_ISBL', 'V_GDS4_ISBL']
    grb_names = ['U component of wind', 'V component of wind']
    data_type = 'uv'
else:
    os.sys.exit('Error: Variable name not supported. Exiting now.')

print('Variable name: {}\nTranslated to: \n{},\n{}.\n'.format(variable_name, variables, grb_names))

print('Creating composites: printing to latest out file in {}.\n'.format(os.getcwd()))

# Load event data
print('Loading event data...')
#strat = 'weekly1238'
#trx_data = np.load(track_dir + 'BoB_genesis_pts_{}.npz'.format(strat))
#strat = 'biweekly1234'
strat = 'weekly1278'
trx_data_fname = track_dir + 'BoB_track_pts_{}.npz'.format(strat)
trx_data = np.load(trx_data_fname)
trx_data = trx_data['arr_0']
print('Data loaded from {}'.format(trx_data_fname))
# Center point of data
centery = (21 + 16) / 2
centerx = (93 + 83) / 2

#256 lat, 512 lon
dx = 180 / 256 # dy = 360 / 512 = dx
dy = dx

# Pick N evenly spaced genesis pts
L = len(trx_data)
#N = 10
N = L
n = 0
I = 0
# List of lag days to gather
lags = [-2, -1, 0, 1, 2]

# Pressure levels
levels = [300, 500, 850]

# Counting
nLevels = len(levels)
nNames = len(grb_names)

# Put all the data in a giant array for processing
data_total = np.zeros( (len(lags), N, nNames, nLevels, 256, 512) )
print('data_total size: {}'.format(data_total.shape))

#print('Testing Memory...')
#data = np.random.random( data_total.shape )
#t, prob = ttest_1samp(data, popmean=0.0, axis=1)
#print('ttest succeeded')
#data_mean = np.mean(data, axis=1)
#print('mean succeeded')
#os.sys.exit()

# Get Harmonics data to use for calculating anomaly
print('Loading seasonal harmonics data...')
harms = {}
for grb_name in grb_names:
    harms[grb_name] = {}
    for l in levels:
        harms[grb_name][str(l)] = np.load('{}ei.oper.an.pl.regn128uv.{}_{}mb_all_years_harmonics.npz'.format(harm_dir, grb_name[0], l))['arr_0']

# Begin looping through events
print('Looping through {} events.'.format(N))
skip = int(L/N)
if skip == 0:
    skip = 1
while I < L and n < N:
    # Read genesis point data
    lon, lat, year, month, day, hour, intensity = trx_data[I, :]
    hour_s = fix_num(hour)
    gen_day = dt.datetime(int(year), int(month), int(day), int(hour))
    # Get number of data points for the year
    dates_b4_gen = pd.date_range(start=fix_num(year)+"-1-1", end=gen_day, freq='6H')
    nDatesB4Gen = len(dates_b4_gen)

    # Loop lags
    for J in range(len(lags)):
        lag = lags[J]
        lag_day = gen_day + dt.timedelta(days=lag)
        year_s = fix_num(lag_day.year)
        month_s = fix_num(lag_day.month)
        day_s = fix_num(lag_day.day)
        # Find the corresponding ERA data file
        idir = '/global/scratch/williamb/data/ERA_Interim/{}/'.format(year_s)
        ifile_name = 'ei.oper.an.pl.regn128{}.{}{}{}{}'.format(data_type, year_s, month_s, day_s, hour_s)
        ifile = idir + ifile_name

        # Convert to netCDF if haven't before
        print('Reading {}'.format(ifile))

        # Read the ERA data and pick the variable we want
        grbs = pygrib.open(ifile)
        grb = grbs.select(name=grb_names, level=levels)

        # Transfer data to an array
        genpt_data = np.zeros( (nNames, nLevels, 256, 512) )
        for K in range(nNames * nLevels):
            i = K % nNames
            j = K % nLevels
            genpt_data[i, j, :, :] = grb[K].values

        # Find mean from this day and subtract it
        time_index = nDatesB4Gen + 4 * lags[J] - 1
        for i in range(nNames):
            for j in range(nLevels):
                genpt_data[i, j, :, :] -= harms[grb_names[i]][str(levels[j])][time_index, :, :]

        # Shift to have genesis (lon, lat) at (180, 0)
        genpt_data = np.roll(genpt_data, - int((90-lat)/dx) + int((90-centery)/dy), axis=1)      # move lat to 90, then move to centery 
        genpt_data = np.roll(genpt_data, - int(lon/dx)      + int(centerx/dx), axis=2)      # move lon to 0, then move to centerx

        # Add ERA data at this genesis point to full data array
        data_total[J, n, :, :, :, :] = genpt_data[:, :, :, :]

        print('Lag {} complete.'.format(lag))

    # Next
    I += skip
    n += 1
    print('{} events composited.\n'.format(n))
    if I == L and n < N:
        print('Ran out of genesis points before getting {} data points!'.format(N))
        data_total = data_total[:, :n, :, :, :, :]
        print('New data_total size: {}'.format(data_total.shape))

print('\nDone. Processing data...')

# Get Lats/Lons from last grib
lats, lons = grb[0].latlons()
lats = lats[:, 0]
lons = lons[0, :]

# t test -- null hypothesis that (data - daily mean) not different from 0.0
# going along axis 1 (the n gen pts) thus outputs prob array of size (lags, names, z, y, x)
t, prob = ttest_1samp(data_total, popmean=0.0, axis=1)

# First save raw data
data_total_composite = np.mean(data_total, axis=1)
for i in range(len(lags)):
    comp_ds = xr.Dataset({'u': (['z', 'y', 'x'], data_total_composite[i, 0, :, :, :]), 
                          'v': (['z', 'y', 'x'], data_total_composite[i, 1, :, :, :])}, 
                        coords={'lon': (['x'], lons),
                                'lat': (['y'], lats),
                                'lvl': (['z'], levels)})
    # Save for plotting
    fname = '{}composite_n{}_{}_{}_lag{}.nc'.format(composite_dir, n, variable_name, strat, lags[i])
    comp_ds.to_netcdf(fname)
    print('Data saved in {}.'.format(fname))

# Now save ttest data: Zero out insignificant data points
data_total_composite[np.where(prob > 0.05)] = 0.0
for i in range(len(lags)):
    comp_ds = xr.Dataset({'u': (['z', 'y', 'x'], data_total_composite[i, 0, :, :, :]), 
                          'v': (['z', 'y', 'x'], data_total_composite[i, 1, :, :, :])}, 
                        coords={'lon': (['x'], lons),
                                'lat': (['y'], lats),
                                'lvl': (['z'], levels)})

    # Save for plotting
    fname = '{}composite_n{}_{}_{}_ttest_lag{}.nc'.format(composite_dir, n, variable_name, strat, lags[i])
    comp_ds.to_netcdf(fname)
    print('Data saved in {}.'.format(fname))

print('\nJob complete, terminating program.')
