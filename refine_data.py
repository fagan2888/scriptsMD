#!/usr/bin/env python3

################################################################################
# Take track data from Boos and just get genesis points from 
#  Bay of Bengal
################################################################################

import numpy as np

# Load data
geninits_NH = np.loadtxt('/global/scratch/hpeter/public_trackdata/geninits_NH.dat',
    dtype='int')
trx_NH = np.loadtxt('/global/scratch/hpeter/public_trackdata/trx_NH.dat')

# get just genesis points
data = trx_NH[geninits_NH - 1, :]
# narrow long
data = data[np.where(np.logical_and(data[:,0] >= 75, data[:,0] <= 95))] 
# narrow lat
data = data[np.where(np.logical_and(data[:,1] >= 10, data[:,1] <= 27))] 
# narrow time
data = data[np.where(np.logical_and(data[:,3] >= 6, data[:,3] <= 9))]

np.savetxt('BoB_genesis_points.dat', data)
