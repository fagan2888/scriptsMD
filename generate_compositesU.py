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

print('Creating composites: printing to latest out file in {}.\n'.format(os.getcwd()))

# Load genesis point data
print('Loading genesis point data...')
strat = 'weekly1238'
trx_data = np.load(track_dir + 'BoB_genesis_pts_{}.npz'.format(strat))
trx_data = trx_data['arr_0']

#256 lat, 512 lon
dx = 180 / 256 # dy = 360 / 512 = dx

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
variable      = 'U_GDS4_ISBL'
variable_name = 'U'
grb_name      = 'U component of wind'
if variable in ['U_GDS4_ISBL', 'V_GDS4_ISBL']:
    data_type = 'uv'
else:
    data_type = 'sc'

# Pick N evenly spaced genesis pts
N = 150
n = 0
L = len(trx_data)
I = 0
# List of lag days to gather
lags = [-2, -1, 0, 1, 2]


# Pressure levels
levels = [850, 500, 300]

# Put all the data in a giant array for processing
data_total = np.zeros( (len(lags), N, len(levels), 256, 512) )

# Get Harmonics data to use for calculating anomaly
print('Loading seasonal harmonics data...')
harms = {}
for l in levels:
    harms[str(l)] = np.load('{}ei.oper.an.pl.regn128uv.{}_{}mb_all_years_harmonics.npz'.format(harm_dir, variable_name, l))['arr_0']


# Begin looping through genesis points
print('Looping through {} genesis points.'.format(N))
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
    for j in range(len(lags)):
        lag = lags[j]
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
        grb = grbs.select(name=grb_name, level=levels)
        # Transfer data to an array
        genpt_data = np.zeros( (len(levels), 256, 512) )
        for i in range(len(levels)):
            genpt_data[i, :, :] = grb[len(levels) - i - 1].values

        # Find mean from this day and subtract it
        index = nDatesB4Gen + 4 * lags[j] - 1
        for i in range(len(levels)):
            genpt_data[i, :, :] -= harms[str(levels[i])][index, :, :]

        # Shift to have genesis (lon, lat) at (180, 0)
        genpt_data = np.roll(genpt_data, 128-int((90-lat)/dx), axis=1) # move lat to 90, then add 128 (256/2) to put at 0
        genpt_data = np.roll(genpt_data, 256-int(lon/dx), axis=2)      # move lon to 0, then add 256 (512/2) to put at 180

        # Add ERA data at this genesis point to full data array
        data_total[j, n, :, :, :] = genpt_data[:, :, :]

        print('Lag {} complete.'.format(lag))

    # Next
    I += skip
    n += 1
    print('{} genesis points composited.\n'.format(n))
    if I == L:
        print('Ran out of genesis points before getting {} data points!'.format(N))

print('Done. Processing data...')

# Get Lats/Lons from last grib
lats, lons = grb[0].latlons()
lats = lats[:, 0]
lons = lons[0, :]

# t test -- null hypothesis that (data - daily mean) not different from 0.0
# going along axis 1 (the N gen pts) thus outputs prob array of size (lags, z, y, x)
t, prob = ttest_1samp(data_total, popmean=0.0, axis=1)

# First save raw data
data_total_composite = np.mean(data_total, axis=1)
for i in range(len(lags)):
    comp_da = xr.DataArray(data_total_composite[i, :, :, :], 
                coords=[levels, lats, lons], 
                dims=['lvl', 'lat', 'lon'])
    # Save for plotting
    fname = '{}composite_n{}_{}_{}_lag{}.nc'.format(composite_dir, n, variable_name, strat, lags[i])
    comp_da.to_netcdf(fname)
    print('Data saved in {}.'.format(fname))

# Now save ttest data: Zero out insignificant data points
data_total_composite[np.where(prob > 0.05)] = 0.0
for i in range(len(lags)):
    comp_da = xr.DataArray(data_total_composite[i, :, :, :], 
                coords=[levels, lats, lons], 
                dims=['lvl', 'lat', 'lon'])

    # Save for plotting
    fname = '{}composite_n{}_{}_{}_ttest_lag{}.nc'.format(composite_dir, n, variable_name, strat, lags[i])
    comp_da.to_netcdf(fname)
    print('Data saved in {}.'.format(fname))

print('\nJob complete, terminating program.')
