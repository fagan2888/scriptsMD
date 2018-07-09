#!/usr/bin/env python3

################################################################################
# This script will likely be more used as a library (hence the 'if' statement
#   below) to calculate seasonal harmonics in a timeseries.
################################################################################
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def remove_harmonics(timeseries, sampfreq, n):
    '''
    IN:
        timeseries: (time x space) data
        sampfreq:   samples per day
        n:          # of harmonics to remove
    OUT:
        dmean: harmonics fit to data
        danom: data minus harmonics fit
    '''

    sp = timeseries.shape
    ntimes = sp[0]

    year = np.arange(ntimes) / (365 * sampfreq)

    # Specify design matrix: Remove first n harmonics
    A = np.zeros((ntimes, 2*n + 1))
    A[:, 0] = np.ones(ntimes)
    for i in range(n):
        A[:, 2*i + 1] = np.cos(2 * np.pi * (i+1) * year)
        A[:, 2*i + 2] = np.sin(2 * np.pi * (i+1) * year)
    dmean = np.zeros(sp)
    for i in range(0, sp[1]):
        # Fit harmonics to each spatial column
        # least squares solving A*x = d with d = ith column of timeseries
        coeff, residuals, rank, s = np.linalg.lstsq(A, timeseries[:, i])
        # A*coeff ~= d
        dmean[:, i] = np.matmul(A, coeff)
    danom = timeseries - dmean

    return dmean, danom

if __name__ == '__main__':
    # year_data.v.values.shape = (time, lat, lon)
    year_data = xr.open_dataset('/home/hpeter/Research2018/MD_files/test_data/ei.oper.an.pl.regn128uv.850mb1997.nc')

    v = year_data.v.values
    sp = v.shape
    v = np.reshape(v, (sp[0], sp[1]*sp[2]))

    vmean, vanom = remove_harmonics(v, 4, 3)

    vmean = np.reshape(vmean, sp)
    vanom = np.reshape(vanom, sp)
    np.savez('vmean.npz', vmean)

    plt.figure()
    plt.contourf(np.roll(vmean[0, :, :], 256, axis=1), origin='upper')
    plt.colorbar()

    plt.figure()
    plt.contourf(np.roll(vmean[700, :, :], 256, axis=1), origin='upper')
    plt.colorbar()

    plt.show()
