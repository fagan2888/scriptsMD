#!/usr/bin/env python3

################################################################################
# Take track data from Boos and just look at events from
#  Bay of Bengal. Then stratify by Kiladis EOF's, intensity, etc.
################################################################################

import numpy as np
import matplotlib.pyplot as plt

def get_angle_deg(x, y):
    ''' 
    Calculate angle on unit circle given signed side lengths.
      x: signed length in x dir, float
      y: signed length in y dir, float
    '''
    angle = np.rad2deg(np.arctan(y / x))
    if x < 0:
        angle += 180
    if angle < 0:
        angle += 360
    return angle

def get_phase(angle):
    ''' 
    Get EOF phase 1-8 from angle.
      angle: angle in degrees, float
    ''' 
    for phase in range(1, 9):
        a1 = 270 + (phase - 1) * 45
        a2 = 270 + phase * 45
        if a1 >= 360:
            a1 -= 360
        if a2 > 360:
            a2 -= 360
        if angle > a1 and angle < a2:
            return phase

def stratify_data(data):
    ''' 
    Stratify the Genesis Point data based on EOF strength/phase.
      data:    track data, array
      returns array with each element  [biweekly_phase#, biweekly_strength, weekly_phase#, weekly_strength] for each element of 'data'
    '''
    strat_array = np.zeros( (data.shape[0], 4) )
    for i in range(data.shape[0]):
        lon, lat, year, month, day, hour, intensity = data[i, :]
    
        if int(hour) in [6, 18]:
            # average the PC vals if not on 0 or 12
            hour = 0
            indx = np.where( np.logical_and( np.logical_and( np.logical_and( PC_dates[:, 0] == int(year), PC_dates[:, 1] == int(month)), PC_dates[:, 2] == int(day)), PC_dates[:, 3] == int(hour)) )[0][0]
            eof1, eof2, eof3, eof4, eof5, eof6, eof7, eof8 = ( PC_vals[indx, :] + PC_vals[indx, :]) / 2
        else:
            indx = np.where( np.logical_and( np.logical_and( np.logical_and( PC_dates[:, 0] == int(year), PC_dates[:, 1] == int(month)), PC_dates[:, 2] == int(day)), PC_dates[:, 3] == int(hour)) )[0][0]
            eof1, eof2, eof3, eof4, eof5, eof6, eof7, eof8 = PC_vals[indx, :]
    
        ampBW = np.sqrt(eof1**2 + eof2**2)
        ampW = np.sqrt(eof3**2 + eof4**2)
        angleBW = get_angle_deg(eof1, eof2)
        angleW = get_angle_deg(eof3, eof4)
        phaseBW = get_phase(angleBW)
        phaseW = get_phase(angleW)

        strat_array[i, :] = [ampBW, phaseBW, ampW, phaseW]

    return strat_array

    
# Load data
geninits_NH_fname = '/home/hpeter/Documents/Research2018/MD_files/public_trackdata/geninits_NH.dat'
trx_NH_fname      = '/home/hpeter/Documents/Research2018/MD_files/public_trackdata/trx_NH.dat'
PC_dates_fname    = '/home/hpeter/Documents/Research2018/MD_files/PC_data/olr.230.westw.2x.5s30n.jja.pc_dates.dat'
PC_vals_fname     = '/home/hpeter/Documents/Research2018/MD_files/PC_data/olr.230.westw.2x.5s30n.jja.pc_vals.dat'

#geninits_NH = np.loadtxt(geninits_NH_fname, dtype='int')
trx_NH = np.loadtxt(trx_NH_fname)
PC_dates = np.loadtxt(PC_dates_fname, dtype='int')
PC_vals = np.loadtxt(PC_vals_fname)

# Normalize data
s1 = np.std(PC_vals[:, 0])
s3 = np.std(PC_vals[:, 2])
s5 = np.std(PC_vals[:, 4])
s7 = np.std(PC_vals[:, 6])
PC_vals[:, 0] /= s1
PC_vals[:, 1] /= s1
PC_vals[:, 2] /= s3
PC_vals[:, 3] /= s3
PC_vals[:, 4] /= s5
PC_vals[:, 5] /= s5
PC_vals[:, 6] /= s7
PC_vals[:, 7] /= s7

print('Number of track points: {}'.format(len(trx_NH)))

# Narrow lon
#llon = 75; rlon = 95
llon = 83; rlon = 93
trx_NH = trx_NH[np.where(np.logical_and(trx_NH[:,0] >= llon, trx_NH[:,0] <= rlon))]

# Narrow lat
#llat = 10; ulat = 27
llat = 16; ulat = 21
trx_NH = trx_NH[np.where(np.logical_and(trx_NH[:,1] >= llat, trx_NH[:,1] <= ulat))]

# Narrow time
m1 = 6;  m2 = 9
trx_NH = trx_NH[np.where(np.logical_and(trx_NH[:,3] >= m1, trx_NH[:,3] <= m2))]

# Narrow intensity
min_int = 1
trx_NH = trx_NH[np.where(trx_NH[:,6] >= min_int)]

print('Number of track points in BoB, JJAS, int > 1: {}'.format(len(trx_NH)))

# Stratify
strat_array = stratify_data(trx_NH)

#strat_trx = trx_NH[np.where(np.logical_and( strat_array[:, 0] >= 1.0, strat_array[:, 1] <= 4 ))]
#print('After stratification: {}'.format(strat_trx.shape[0]))

mode = 'weekly'
amp_thresh = 0.5
print(mode)
strong = np.zeros((8))
weak   = np.zeros((8))
for phase in range(1, 9):
    if mode == 'biweekly':
        strong[phase-1] = strat_array[np.where(np.logical_and( strat_array[:, 0] >= amp_thresh, strat_array[:, 1] == phase ))].shape[0]
        weak[phase-1] =   strat_array[np.where(np.logical_and( strat_array[:, 0]  < amp_thresh, strat_array[:, 1] == phase ))].shape[0]
    elif mode == 'weekly':
        strong[phase-1] = strat_array[np.where(np.logical_and( strat_array[:, 2] >= amp_thresh, strat_array[:, 3] == phase ))].shape[0]
        weak[phase-1] =   strat_array[np.where(np.logical_and( strat_array[:, 2]  < amp_thresh, strat_array[:, 3] == phase ))].shape[0]
print('Strong: {}, sum = {}'.format(strong, np.sum(strong)))
print('Weak: {}, sum = {}'.format(weak, np.sum(weak)))
print('total: {}'.format(np.sum(strong+weak)))

#PC_vals = PC_vals[np.where(np.logical_and(PC_dates[:, 1] >= m1, PC_dates[:, 1] <= m2))]
#strong_all = np.zeros((8))
#weak_all   = np.zeros((8))
#for i in range(PC_vals.shape[0]):
#    eof1, eof2, eof3, eof4, eof5, eof6, eof7, eof8 = PC_vals[i, :]
#
#    if mode == 'biweekly':
#        x = eof1; y = eof2
#    elif mode == 'weekly':
#        x = eof3; y = eof4
#
#    amp = np.sqrt(x**2 + y**2)
#    angle = get_angle_deg(x, y)
#
#    if amp > amp_thresh:
#        strong_all[get_phase(angle) - 1] += 1
#    else:
#        weak_all[get_phase(angle) - 1] += 1
#        
##strong /= strong_all
##weak /= weak_all
#strong = strong_all
#weak = weak_all

################################################################################
f = plt.figure(figsize=(8,5))
ax = plt.subplot(111)

#ax.set_title('Distribution of High Intensity MD Tracks\nin BoB JJAS w.r.t. {} EOFs'.format(mode.capitalize()), size=16)
#ax.set_title('Distribution of {} EOFs in JJAS'.format(mode.capitalize()), size=16)
ax.set_xlabel('Phase', size=12)
#ax.set_ylabel('# Track Points / # Total Points in Phase', size=12)
#ax.set_ylabel('# Total Points in Phase', size=12)

centers = np.arange(1, 9)
width = 0.45
ax.bar(centers - width/2, strong, width=width, color='#aa7777', hatch='', edgecolor='k', label='strong phase')
ax.bar(centers + width/2, weak,   width=width, color='#7777aa', hatch='', edgecolor='k', label='weak phase')
ax.legend(loc='lower right')

plt.show()

#f.savefig(mode + '_stratified.png', dpi=120)
#f.savefig(mode + '_stratified_normalized.png', dpi=120)
#f.savefig(mode + '_all.png', dpi=120)
#################################################################################

#f = plt.figure(figsize=(8,5))
#ax = plt.subplot(111)
#
#ax.set_title('Distribution of MD Tracks in BoB w.r.t. Weekly and Biweekly EOFs', size=16)
#ax.set_xlabel('Biweekly Phase', size=12)
#ax.set_ylabel('Weekly Phase', size=12)
#grid = np.zeros( (8, 8) )
#for i in range(trx_NH.shape[0]):
#    lon, lat, year, month, day, hour, intensity = trx_NH[i, :]
#
#    if i % 1000 == 0: print(i)
#
#    if int(hour) in [6, 18]:
#        # average the PC vals if not on 0 or 12
#        hour = 0
#        indx = np.where( np.logical_and( np.logical_and( np.logical_and( PC_dates[:, 0] == int(year), PC_dates[:, 1] == int(month)), PC_dates[:, 2] == int(day)), PC_dates[:, 3] == int(hour)) )[0][0]
#        eof1, eof2, eof3, eof4, eof5, eof6, eof7, eof8 = ( PC_vals[indx, :] + PC_vals[indx, :]) / 2
#    else:
#        indx = np.where( np.logical_and( np.logical_and( np.logical_and( PC_dates[:, 0] == int(year), PC_dates[:, 1] == int(month)), PC_dates[:, 2] == int(day)), PC_dates[:, 3] == int(hour)) )[0][0]
#        eof1, eof2, eof3, eof4, eof5, eof6, eof7, eof8 = PC_vals[indx, :]
#
#    ampBW = np.sqrt(eof1**2 + eof2**2)
#    ampW = np.sqrt(eof3**2 + eof4**2)
#    angleBW = get_angle_deg(eof1, eof2)
#    angleW = get_angle_deg(eof3, eof4)
#    phaseBW = get_phase(angleBW)
#    phaseW = get_phase(angleW)
#
#    grid[phaseBW -1 , phaseW - 1] += 1
#
