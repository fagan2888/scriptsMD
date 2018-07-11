#!/usr/bin/env python3

################################################################################
# Generate composites for Monsoon Depression genesis points
# in the Bay of Bengal (data from Boos)
################################################################################

import numpy as np
import pandas as pd
import xarray as xr
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
data_gen = np.load(track_dir + 'BoB_gen_pts_JJAS_strat_biweekly_strong_p1234.npz')
data_gen = data_gen['arr_0']

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
variable      = 'V_GDS4_ISBL'
variable_name = 'V'
if variable in ['U_GDS4_ISBL', 'V_GDS4_ISBL']:
    data_type = 'uv'
else:
    data_type = 'sc'

# Pick N evenly spaced genesis pts
N = 150
n = 0
L = len(data_gen)
i = 0
# List of lag days to gather
lags = [-2, -1, 0, 1, 2]

print('Gathering ERA data now.')

# Get Harmonics data to use for calculating anomaly
harms = {}
harms['850'] = np.load('{}ei.oper.an.pl.regn128uv.{}_{}mb_all_years_harmonics.npz'.format(harm_dir, variable_name, 850))['arr_0']
harms['500'] = np.load('{}ei.oper.an.pl.regn128uv.{}_{}mb_all_years_harmonics.npz'.format(harm_dir, variable_name, 500))['arr_0']
harms['300'] = np.load('{}ei.oper.an.pl.regn128uv.{}_{}mb_all_years_harmonics.npz'.format(harm_dir, variable_name, 300))['arr_0']

# Put all the data in a giant array for processing
data_total = np.zeros( (len(lags), N, len(harms.keys()), 256, 512) )

# Begin looping through genesis points
while i < L and n < N:
    # Read genesis point data
    lon, lat, year, month, day, hour, intensity = data_gen[i, :]
    hour_s = fix_num(hour)
    gen_day = dt.datetime(int(year), int(month), int(day), int(hour))
    # Get number of data points for the year
    data_pts = pd.date_range(start=fix_num(year)+"-1-1", end=gen_day, freq='6H')
    nDataPts = len(data_pts)

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
        if os.path.isfile(data_dir + ifile_name + '.nc'):
            print('{} exists.'.format(ifile_name+'.nc'))
        else:
            print('Converting .grb to .nc')
            os.system('ncl_convert2nc {}.grb'.format(ifile))
            os.system('mv {}.nc {}'.format(ifile_name, data_dir))

        # Read the ERA data and pick the variable we want
        ds = xr.open_dataset(data_dir + ifile_name + '.nc')
        # becomes data array:
        genpt_da = ds[variable]
        genpt_da = genpt_da.sel(lv_ISBL0=[850, 500, 300])

        if n == 0:
            # Create an unshifted copy for later
            comp_da = genpt_da

        # Find mean from this day and subtract it
        index = nDataPts + 4 * lags[j] - 1

        genpt_da.values = [genpt_da.sel(lv_ISBL0=850).values - harms['850'][index, :, :],
                           genpt_da.sel(lv_ISBL0=500).values - harms['500'][index, :, :],
                           genpt_da.sel(lv_ISBL0=300).values - harms['300'][index, :, :]]

        # Shift to have genesis (lon, lat) at (180, 0)
        genpt_da = genpt_da.roll(g4_lat_1=128-int((90-lat)/dx)) # move lat to 90, then add 128 (256/2) to put at 0
        genpt_da = genpt_da.roll(g4_lon_2=256-int(lon/dx))      # move lon to 0, then add 256 (512/2) to put at 180

        # Add ERA data at this genesis point to full data array
        data_total[j, n, :, :, :] = genpt_da.values[:, :, :]

        print('Lag {} complete.'.format(lag))

    # Next
    i += int(L/N)
    n += 1
    print('{} genesis points composited.\n'.format(n))
    if i > L:
        print('Ran out of genesis points before getting {} data points!'.format(N))

print('Done. Processing data...')


# t test -- null hypothesis that (data - daily mean) not different from 0.0
# going along axis 1 (the N gen pts) thus outputs prob array of size (lags, z, y, x)
t, prob = ttest_1samp(data_total, popmean=0.0, axis=1)

# First save raw data
data_total_composite = np.mean(data_total, axis=1)
for i in range(len(lags)):
    comp_da.values[:, :, :] = data_total_composite[i, :, :, :]

    # Save for plotting
    fname = '{}composite_n{}_{}_JJAS_strat_harm_lag{}.nc'.format(composite_dir, N, variable_name, lags[i])
    comp_da.to_netcdf(fname)
    print('Data saved in {}.'.format(fname))

# Now save ttest data: Zero out insignificant data points
data_total_composite[np.where(prob > 0.05)] = 0.0
for i in range(len(lags)):
    comp_da.values[:, :, :] = data_total_composite[i, :, :, :]

    # Save for plotting
    fname = '{}composite_n{}_{}_JJAS_strat_harm_ttest_lag{}.nc'.format(composite_dir, N, variable_name, lags[i])
    comp_da.to_netcdf(fname)
    print('Data saved in {}.'.format(fname))

print('\nJob complete, terminating program.')
