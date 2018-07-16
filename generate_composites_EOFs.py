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
if len(os.sys.argv) != 2:
    os.sys.exit('Error: Must pass one argument for variable name. Exiting now.')
variable_name = os.sys.argv[1] 
if variable_name == 'U':
    variable      = 'U_GDS4_ISBL'
    grb_name      = 'U component of wind'
elif variable_name == 'V':
    variable      = 'V_GDS4_ISBL'
    grb_name      = 'V component of wind'
else:
    os.sys.exit('Error: Variable name not supported. Exiting now.')
print('Variable name given as input: {}\nTranslated to {}, {}.\n'.format(variable_name, variable, grb_name))
if variable in ['U_GDS4_ISBL', 'V_GDS4_ISBL']:
    data_type = 'uv'
else:
    data_type = 'sc'

print('Creating composites: printing to latest out file in {}.\n'.format(os.getcwd()))

# Load EOF dates data
print('Loading EOF dates data...')
strat = 'biweekly1234'
#strat = 'weekly1278'
EOF_dates = np.load(track_dir + 'EOFs_JJAS_{}s.npz'.format(strat))
EOF_dates = EOF_dates['arr_0']

#256 lat, 512 lon
dx = 180 / 256 # dy = 360 / 512 = dx

# Pick N evenly spaced genesis pts
L = len(EOF_dates)
N = L

# Pressure levels
levels = [850, 500, 300]

# keep day in a giant array for ttests 
data_total = np.zeros( (N, len(levels), 256, 512) )
print('data_total size: {}'.format(data_total.shape))

# Get Harmonics data to use for calculating anomaly
print('Loading seasonal harmonics data...')
harms = {}
for l in levels:
    harms[str(l)] = np.load('{}ei.oper.an.pl.regn128uv.{}_{}mb_all_years_harmonics.npz'.format(harm_dir, variable_name, l))['arr_0']

# Begin looping through genesis points
print('Looping through {} dates.'.format(N))
for I in range(N):
    # Read genesis point data
    year, month, day, hour, minute, second, milisecond = EOF_dates[I, :]
    year_s = fix_num(year); month_s = fix_num(month); day_s = fix_num(day); hour_s = fix_num(hour)
    EOF_day = dt.datetime(int(year), int(month), int(day), int(hour))
    # Get number of data points for the year
    dates_b4_EOF_day = pd.date_range(start=fix_num(year)+"-1-1", end=EOF_day, freq='6H')
    nDatesB4EOF = len(dates_b4_EOF_day)

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
    EOF_day_data = np.zeros( (len(levels), 256, 512) )
    for i in range(len(levels)):
        EOF_day_data[i, :, :] = grb[len(levels) - i - 1].values

    # Find mean from this day and subtract it
    index = nDatesB4EOF - 1
    for i in range(len(levels)):
        EOF_day_data[i, :, :] -= harms[str(levels[i])][index, :, :]

    # Add ERA data at this genesis point to full data array
    data_total[I, :, :, :] = EOF_day_data[:, :, :]

    # Next
    print('{} dates composited.\n'.format(I+1))

print('Done. Processing data...')

# Get Lats/Lons from last grib
lats, lons = grb[0].latlons()
lats = lats[:, 0]
lons = lons[0, :]

# t test -- null hypothesis that (data - daily mean) not different from 0.0
# going along axis 1 (the n gen pts) thus outputs prob array of size (lags, z, y, x)
t, prob = ttest_1samp(data_total, popmean=0.0, axis=0)

# First save raw data
data_total_composite = np.mean(data_total, axis=0)
comp_da = xr.DataArray(data_total_composite[:, :, :], 
                        coords=[levels, lats, lons], 
                        dims=['lvl', 'lat', 'lon'])
fname = '{}composite_EOFs_n{}_{}_{}.nc'.format(composite_dir, N, variable_name, strat)
comp_da.to_netcdf(fname)
print('Data saved in {}.'.format(fname))

# Now save ttest data: Zero out insignificant data points
data_total_composite[np.where(prob > 0.05)] = 0.0
comp_da = xr.DataArray(data_total_composite[:, :, :], 
                        coords=[levels, lats, lons], 
                        dims=['lvl', 'lat', 'lon'])
fname = '{}composite_EOFs_n{}_{}_{}_ttest.nc'.format(composite_dir, N, variable_name, strat)
comp_da.to_netcdf(fname)
print('Data saved in {}.'.format(fname))

print('\nJob complete, terminating program.')
