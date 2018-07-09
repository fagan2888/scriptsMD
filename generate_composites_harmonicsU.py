#!/usr/bin/env python3

################################################################################
# Generate composites for Monsoon Depression genesis points
# in the Bay of Bengal (data from Boos)
################################################################################

import numpy as np
import pandas as pd
import xarray as xr
from scipy.io import loadmat
from scipy.stats import ttest_1samp
from remove_harmonics import remove_harmonics
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

def log(s):
    '''
    It is much more convenient to have the program write to a file
    rather than standard output. This function makes that process
    easier.
    '''
    with open(out_file, 'a') as f:
        f.write(s + '\n')
    return 0

scratch_dir = '/global/scratch/hpeter/'
data_dir = scratch_dir + 'data/ERA_Interim/to_netCDF/'
composite_dir = scratch_dir + 'composites/'
track_dir = scratch_dir + 'public_trackdata/'

# Create an output file to print to while running.
# This allows for running the job in the background while still monitoring easily.
i = 0
out_file = composite_dir + 'out_file' + str(i)
while os.path.isfile(out_file):
    out_file = composite_dir + 'out_file' + str(i)
    i += 1

log('Creating composites: logging to {}.\n'.format(out_file))

# Load genesis point data
data_gen_df = pd.read_csv(track_dir + 'BoB_genesis_points.dat', delimiter=' ', names=['lon', 'lat', 'year', 'month', 'day', 'hour', 'intensity'], index_col=False)
#data_gen_df = data_gen_df.sort_values(by='year')
data_gen = data_gen_df.values

#256 lat, 512 lon
dx = 180 / 256 # dy = 360 / 512 = dx

'''
ERA-int Variables:
    U_GDS4_ISBL
    V_GDS4_ISBL
    ---
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
variable = 'U_GDS4_ISBL'
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

log('Gathering ERA data now.')
# Put all the data in a giant array for processing
data_total = np.zeros( (len(lags), N, 2, 256, 512) )
prev_year = 0
while i < L and n < N:
    # Read genesis point data
    lon, lat, year, month, day, hour, intensity = data_gen[i, :]
    hour_s = fix_num(hour)
    gen_day = dt.datetime(int(year), int(month), int(day), int(hour))

    # Get first 3 seasonal harmonics for each new year
    log(str(prev_year))
    log(str(year))
    if prev_year != year:
        year850_ds = xr.open_dataset('/global/scratch/williamb/data/eraint/ei.oper.an.pl.regn128uv.850mb{}.nc'.format(int(year)))
        year300_ds = xr.open_dataset('/global/scratch/williamb/data/eraint/ei.oper.an.pl.regn128uv.300mb{}.nc'.format(int(year)))
        if variable == 'U_GDS4_ISBL':
            # save some data arrays
            da850 = year850_ds.u
            da300 = year300_ds.u
        elif variable == 'V_GDS4_ISBL':
            da850 = year850_ds.v
            da300 = year300_ds.v
        else:
            log('Not ready for non-UV data! Exitting')
            os.sys.exit()
        year850 = da850.values
        year300 = da300.values
        sp = year850.shape
        log('Removing Harmonics: 850 mb')
        year850_mean, year850_anom = remove_harmonics(np.reshape(year850, (sp[0], sp[1]*sp[2])), sampfreq=4, n=3)
        year850_mean = np.reshape(year850_mean, sp)
        sp = year300.shape
        log('Removing Harmonics: 300 mb')
        year300_mean, year300_anom = remove_harmonics(np.reshape(year300, (sp[0], sp[1]*sp[2])), sampfreq=4, n=3)
        year300_mean = np.reshape(year300_mean, sp)
        # save means to data arrays to use date searching functionality
        da850.values = year850_mean
        da300.values = year300_mean
        prev_year = year

    # Loop lags
    for j in range(len(lags)):
        lag = lags[j]
        lag_day = gen_day + dt.timedelta(days=lag)
        year_s = fix_num(lag_day.year)
        month_s = fix_num(lag_day.month)
        day_s = fix_num(lag_day.day)
        # Find the corresponding ERA data file
        idir = '/global/scratch/williamb/data/ERA_Interim/{}/'.format(year_s)
        odir = data_dir
        ifile_name = 'ei.oper.an.pl.regn128{}.{}{}{}{}'.format(data_type, year_s, month_s, day_s, hour_s)
        ifile = idir + ifile_name

        # Convert to netCDF if haven't before
        if os.path.isfile(odir + ifile_name + '.nc'):
            log('{} exists.'.format(ifile_name+'.nc'))
        else:
            log('Converting .grb to .nc')
            os.system('ncl_convert2nc {}.grb'.format(ifile))
            os.system('mv {}.nc {}'.format(ifile_name, odir))

        # Read the ERA data and pick the variable we want
        ds = xr.open_dataset(odir + ifile_name + '.nc')
        # becomes data array:
        genpt_da = ds[variable]
        genpt_da = genpt_da.sel(lv_ISBL0=[850, 300])

        if n == 0:
            # Create an unshifted copy for later
            comp_da = genpt_da

        # Subtract out mean from this day
        genpt_da.values = [genpt_da.sel(lv_ISBL0=850).values - da850.sel(time=lag_day).values, 
                            genpt_da.sel(lv_ISBL0=300).values - da300.sel(time=lag_day).values]

        # Shift to have genesis (lon, lat) at (180, 0)
        genpt_da = genpt_da.roll(g4_lat_1=128-int((90-lat)/dx)) # move lat to 90, then add 128 (256/2) to put at 0
        genpt_da = genpt_da.roll(g4_lon_2=256-int(lon/dx))      # move lon to 0, then add 256 (512/2) to put at 180

        # Add ERA data at this genesis point to full data array
        data_total[j, n, :, :, :] = genpt_da.values[:, :, :]

        log('Lag {} complete.'.format(lag))

    # Next
    i += int(L/N)
    n += 1
    log('{} genesis points composited.\n'.format(n))
    if i > L:
        log('Ran out of genesis points before getting {} data points!'.format(N))

log('Done. Processing data...')


# t test -- null hypothesis that (data - daily mean) not different from 0.0
# going along axis 1 (the N gen pts) thus outputs prob array of size (lags, z, y, x)
t, prob = ttest_1samp(data_total, popmean=0.0, axis=1)

# First save raw data
data_total_composite = np.mean(data_total, axis=1)
for i in range(len(lags)):
    comp_da.values[:, :, :] = data_total_composite[i, :, :, :]

    # Save for plotting
    fname = '{}composite_n{}_{}_JJAS_lag{}_harmonics.nc'.format(composite_dir, N, variable, lags[i])
    comp_da.to_netcdf(fname)
    log('Data saved in {}.'.format(fname))

# Now save ttest data: Zero out insignificant data points
data_total_composite[np.where(prob > 0.05)] = 0.0
for i in range(len(lags)):
    comp_da.values[:, :, :] = data_total_composite[i, :, :, :]

    # Save for plotting
    fname = '{}composite_n{}_{}_JJAS_lag{}_ttest_harmonics.nc'.format(composite_dir, N, variable, lags[i])
    comp_da.to_netcdf(fname)
    log('Data saved in {}.'.format(fname))

log('\nTask complete, terminating program.')
