#!/usr/bin/env python3

################################################################################
# This script takes combined ERA-Interim data (see 'combine_data.py') from
#   multiple years and avereges for each calendar day. The data is then fit
#   by 'remove_harmonics.py' for a ~30 year average that can be reused.
################################################################################

import pandas as pd
import xarray as xr
import numpy as np
from remove_harmonics import remove_harmonics

idir = '/global/scratch/users/williamb/data/eraint/'
odir = '/global/scratch/users/hpeter/data/ERA_Interim/'
level = 300
name = 'V'
data_type = 'uv'
#data_type = 'sc'
ofile1 = '{}all_years/ei.oper.an.pl.regn128{}.{}_{}mb_all_years.nc'.format(odir, data_type, name, level)
ofile2 = '{}all_years/ei.oper.an.pl.regn128{}.{}_{}mb_all_years_harmonics.npz'.format(odir, data_type, name, level)

leap_data = np.zeros((4, 256, 512))
not_leap_data = np.zeros((365*4 - 3, 256, 512))

nLeap = 0
nNotLeap = 0

for year in range(1979, 2017):  # 2017 only has half a year, stop at 2016

    #ifile = '{}{}/ei.oper.an.pl.regn128sc.{}_{}mb_{}.nc'.format(odir, year, name, level, year)
    #da = xr.open_dataarray(ifile)
    
    ifile = '{}/ei.oper.an.pl.regn128uv.{}mb{}.nc'.format(idir, level, year)
    ds = xr.open_dataset(ifile)
    da = ds.u

    leap_yr = True
    try: 
        leap_dates = pd.date_range(start=str(year)+'-2-29 00:00', end=str(year)+'-2-29  18:00', freq='6H')
    except:
        leap_yr = False
    
    if leap_yr:
        not_leap_dates1 = pd.date_range(start=str(year)+'-1-1  00:00', end=str(year)+'-2-28  18:00', freq='6H')
        not_leap_dates2 = pd.date_range(start=str(year)+'-3-1  00:00', end=str(year)+'-12-31 00:00', freq='6H')

        leap = da.sel(time=leap_dates).values
        not_leap1 = da.sel(time=not_leap_dates1).values
        not_leap2 = da.sel(time=not_leap_dates2).values
        N = not_leap1.shape[0]

        leap_data += leap
        not_leap_data[:N, :, :] += not_leap1
        not_leap_data[N:, :, :] += not_leap2

        nLeap += 1
        nNotLeap += 1
    else:
        dates = pd.date_range(start=str(year)+'-1-1  00:00', end=str(year)+'-12-31 00:00', freq='6H')

        not_leap = da.sel(time=dates).values

        not_leap_data += not_leap

        nNotLeap += 1

leap_data /= nLeap
not_leap_data /= nNotLeap

full_data = np.zeros( (365*4 - 3 + 4, 256, 512) )
full_data[:N, :, :] = not_leap_data[:N, :, :]
full_data[N:N+4, :, :] = leap_data[:, :, :]
full_data[N+4:, :, :] = not_leap_data[N:, :, :]

# get dummy netCDF file from 2016 since it was a leap year
year = 2016
ifile = '{}{}/ei.oper.an.pl.regn128sc.Geopotential_850mb_{}.nc'.format(odir, year, year)
da_avg = xr.open_dataarray(ifile)

# save the data 
da_avg.values[:, :, :] = full_data[:, :, :]
da_avg.to_netcdf(ofile1)

# do the harmonics calculation
sp = full_data.shape
mean, anom = remove_harmonics(np.reshape(full_data, (sp[0], sp[1]*sp[2])), sampfreq=4, n=3)
mean = np.reshape(mean, sp)

# save the data
np.savez(ofile2, mean)
