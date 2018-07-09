#!/usr/bin/env python3

################################################################################
# Take track data from Boos and just get genesis points from
#  Bay of Bengal. Then stratify by Kiladis EOF's
################################################################################

import numpy as np

odir = '../MD_files/public_trackdata/'
#odir = '/global/home/users/hpeter/public_trackdata/'

# Load data
geninits_NH_fname = '../MD_files/public_trackdata/geninits_NH.dat'
trx_NH_fname      = '../MD_files/public_trackdata/trx_NH.dat'
PC_dates_fname    = '../MD_files/PC_data/olr.230.westw.2x.5s30n.jja.pc_dates.dat'
PC_vals_fname     = '../MD_files/PC_data/olr.230.westw.2x.5s30n.jja.pc_vals.dat'

#geninits_NH_fname = '/global/home/users/hpeter/public_trackdata/geninits_nh.dat'
#trx_nh_fname      = '/global/home/users/hpeter/public_trackdata/trx_NH.dat'
#PC_dates_fname    = '/global/home/users/hpeter/PC_data/olr.230.westw.2x.5s30n.jja.pc_dates.dat'
#PC_vals_fname     = '/global/home/users/hpeter/PC_data/olr.230.westw.2x.5s30n.jja.pc_vals.dat'

geninits_NH = np.loadtxt(geninits_NH_fname, dtype='int')
trx_NH = np.loadtxt(trx_NH_fname)
PC_dates = np.loadtxt(PC_dates_fname, dtype='int')
PC_vals = np.loadtxt(PC_vals_fname)

# Get just genesis points
gen_pts = trx_NH[geninits_NH - 1, :]

# Narrow lon
llon = 75; rlon = 95
gen_pts = gen_pts[np.where(np.logical_and(gen_pts[:,0] >= llon, gen_pts[:,0] <= rlon))]

# Narrow lat
llat = 10; ulat = 27
gen_pts = gen_pts[np.where(np.logical_and(gen_pts[:,1] >= llat, gen_pts[:,1] <= ulat))]

# Narrow time
m1 = 6;  m2 = 9
gen_pts = gen_pts[np.where(np.logical_and(gen_pts[:,3] >= m1, gen_pts[:,3] <= m2))]

# Calculate phase and amplitude of PCs
strat_indexes = np.zeros(gen_pts.shape[0], dtype=bool)
for i in range(len(gen_pts)):
    lon, lat, year, month, day, hour, intensity = gen_pts[i, :]

    if int(hour) == 6:
        hour = 0
    elif int(hour) == 18:
        hour = 12

    indx = np.where( np.logical_and( np.logical_and( np.logical_and( PC_dates[:, 0] == int(year), PC_dates[:, 1] == int(month)), PC_dates[:, 2] == int(day)), PC_dates[:, 3] == int(hour)) )[0][0]

    eof1, eof2, eof3, eof4, eof5, eof6, eof7, eof8 = PC_vals[indx, :]
    amp = np.sqrt(eof1**2 + eof2**2)

    if amp > 1.0 and eof1 > 0:
        strat_indexes[i] = True

data = gen_pts[strat_indexes, :]
np.savez(odir + 'BoB_genesis_points_stratified.npz', data)
