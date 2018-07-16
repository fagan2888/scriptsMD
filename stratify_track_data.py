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
            eof1, eof2, eof3, eof4, eof5, eof6, eof7, eof8 = (PC_vals[indx, :] + PC_vals[indx + 1, :]) / 2
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

geninits_NH = np.loadtxt(geninits_NH_fname, dtype='int')
trx_NH = np.loadtxt(trx_NH_fname)
PC_dates = np.loadtxt(PC_dates_fname, dtype='int')
PC_vals = np.loadtxt(PC_vals_fname)

points_type = 'genesis'
#points_type = 'track'
if points_type == 'genesis':
    MD_set = trx_NH[geninits_NH - 1, :]
elif points_type == 'track':
    MD_set = trx_NH

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

print('Total num data points: {}'.format(PC_dates.shape[0]))
m1 = 6; m2 = 9
PC_vals = PC_vals[np.where(np.logical_and(PC_dates[:, 1] >= m1, PC_dates[:, 1] <= m2))]
PC_dates = PC_dates[np.where(np.logical_and(PC_dates[:, 1] >= m1, PC_dates[:, 1] <= m2))]
PC_vals = PC_vals[np.where(PC_dates[:, 0] < 2017)]
PC_dates = PC_dates[np.where(PC_dates[:, 0] < 2017)]
print('Num data points JJAS: {}'.format(PC_dates.shape[0]))
mode = 'biweekly'; phases = '1678'
#mode = 'weekly'; phases = '1238'
amp_thresh = 0.5
keep = np.zeros(PC_dates.shape[0], dtype='bool')
for i in range(PC_dates.shape[0]):
    eof1, eof2, eof3, eof4, eof5, eof6, eof7, eof8 = PC_vals[i, :]
    if mode == 'biweekly':
        amp = np.sqrt(eof1**2 + eof2**2)
        angle = get_angle_deg(eof1, eof2)
    elif mode == 'weekly':
        amp = np.sqrt(eof3**2 + eof4**2)
        angle = get_angle_deg(eof3, eof4)
    phase = get_phase(angle)
    if amp >= amp_thresh and str(phase) in phases:
        keep[i] = True
PC_dates = PC_dates[keep]
print('Num data points JJAS, mode={}, phases={}, amp_thresh={}: {}'.format(mode, phases, amp_thresh, PC_dates.shape[0]))
print(PC_dates)
fname = 'EOFs_JJAS_{}{}s.npz'.format(mode, phases)
np.savez(fname, PC_dates)
print('Saved to {}.'.format(fname))

    

#print('Number of {} points: {}'.format(points_type, MD_set.shape[0]))
#
## Narrow lon
##llon = 75; rlon = 95
#llon = 83; rlon = 93
#MD_set = MD_set[np.where(np.logical_and(MD_set[:,0] >= llon, MD_set[:,0] <= rlon))]
#
## Narrow lat
##llat = 10; ulat = 27
#llat = 16; ulat = 21
#MD_set = MD_set[np.where(np.logical_and(MD_set[:,1] >= llat, MD_set[:,1] <= ulat))]
#
## Narrow time
#m1 = 6;  m2 = 9
#MD_set = MD_set[np.where(np.logical_and(MD_set[:,3] >= m1, MD_set[:,3] <= m2))]
#PC_vals  = PC_vals[ np.where(np.logical_and(PC_dates[:, 1] >= m1, PC_dates[:, 1] <= m2))]
#PC_dates = PC_dates[np.where(np.logical_and(PC_dates[:, 1] >= m1, PC_dates[:, 1] <= m2))]
#
## Narrow intensity
#min_int = 1
#MD_set = MD_set[np.where(MD_set[:,6] >= min_int)]
#
#print('Number of {} points within constraints: {}'.format(points_type, MD_set.shape[0]))
#
## Stratify
#amp_thresh = 0.5
#strat_array = stratify_data(MD_set)
#
##mode = 'biweekly'
##phases = '1678'
#mode = 'weekly'
#phases = '1238'
#if mode == 'biweekly':
#    strat = MD_set[np.where(np.logical_and( strat_array[:, 0] >= amp_thresh, 
#                            np.logical_or(strat_array[:, 1] == int(phases[0]),
#                            np.logical_or(strat_array[:, 1] == int(phases[1]),
#                            np.logical_or(strat_array[:, 1] == int(phases[2]),
#                                          strat_array[:, 1] == int(phases[3]))))))]
#elif mode == 'weekly':
#    strat = MD_set[np.where(np.logical_and( strat_array[:, 2] >= amp_thresh, 
#                            np.logical_or(strat_array[:, 3] == int(phases[0]),
#                            np.logical_or(strat_array[:, 3] == int(phases[1]),
#                            np.logical_or(strat_array[:, 3] == int(phases[2]),
#                                          strat_array[:, 3] == int(phases[3]))))))]
#print('After stratification: {}'.format(strat.shape[0]))
#fname = 'BoB_{}_pts_{}{}.npz'.format(points_type, mode, phases)
#np.savez(fname, strat)
#print('Saved to {}.'.format(fname))
#
###normalize = True
##normalize = False
###show_season_totals = True
##show_season_totals = False
##for mode in ['weekly', 'biweekly']:
##    print('Mode: {}'.format(mode))
##    strong = np.zeros((8))
##    weak   = np.zeros((8))
##    for phase in range(1, 9):
##        if mode == 'biweekly':
##            strong[phase-1] = strat_array[np.where(np.logical_and( strat_array[:, 0] >= amp_thresh, strat_array[:, 1] == phase ))].shape[0]
##            weak[phase-1] =   strat_array[np.where(np.logical_and( strat_array[:, 0]  < amp_thresh, strat_array[:, 1] == phase ))].shape[0]
##        elif mode == 'weekly':
##            strong[phase-1] = strat_array[np.where(np.logical_and( strat_array[:, 2] >= amp_thresh, strat_array[:, 3] == phase ))].shape[0]
##            weak[phase-1] =   strat_array[np.where(np.logical_and( strat_array[:, 2]  < amp_thresh, strat_array[:, 3] == phase ))].shape[0]
##    print('Strong: {}, sum = {}'.format(strong, np.sum(strong)))
##    print('Weak: {}, sum = {}'.format(weak, np.sum(weak)))
##    print('total: {}'.format(np.sum(strong+weak)))
##   
##    if normalize or show_season_totals:
##        strong_all = np.zeros((8))
##        weak_all   = np.zeros((8))
##        for i in range(PC_vals.shape[0]):
##            eof1, eof2, eof3, eof4, eof5, eof6, eof7, eof8 = PC_vals[i, :]
##        
##            if mode == 'biweekly':
##                x = eof1; y = eof2
##            elif mode == 'weekly':
##                x = eof3; y = eof4
##        
##            amp = np.sqrt(x**2 + y**2)
##            angle = get_angle_deg(x, y)
##        
##            if amp > amp_thresh:
##                strong_all[get_phase(angle) - 1] += 1
##            else:
##                weak_all[get_phase(angle) - 1] += 1
##        if normalize: 
##            strong /= strong_all
##            weak /= weak_all
##   
##    if show_season_totals:
##        f = plt.figure(figsize=(8,5))
##        ax = plt.subplot(111)
##        
##        ax.set_title('Distribution of {} EOFs in JJAS'.format(mode.capitalize()), size=16)
##        ax.set_xlabel('Phase', size=12)
##        ax.set_ylabel('# Total Points in Phase', size=12)
##        
##        centers = np.arange(1, 9)
##        width = 0.45
##        ax.bar(centers - width/2, strong_all, width=width, color='#aa7777', hatch='', edgecolor='k', label='strong phase')
##        ax.bar(centers + width/2, weak_all,   width=width, color='#7777aa', hatch='', edgecolor='k', label='weak phase')
##        ax.legend(loc='center right')
##        
##        fname = 'bars_seasonal_totals_{}.png'.format(mode)
##        plt.savefig(fname, dpi=120)
##
##    f = plt.figure(figsize=(8,5))
##    ax = plt.subplot(111)
##    
##    ax.set_title('Distribution of MD {} points\nin BoB JJAS: {} EOFs'.format(points_type.capitalize(), mode.capitalize()), size=16)
##    ax.set_xlabel('Phase', size=12)
##    if normalize:
##        ax.set_ylabel('# Track Points / # Total Points in Phase', size=12)
##    else:
##        ax.set_ylabel('# Total Points in Phase', size=12)
##    
##    centers = np.arange(1, 9)
##    width = 0.45
##    ax.bar(centers - width/2, strong, width=width, color='#aa7777', hatch='', edgecolor='k', label='strong phase')
##    ax.bar(centers + width/2, weak,   width=width, color='#7777aa', hatch='', edgecolor='k', label='weak phase')
##    ax.legend(loc='center right')
##
##    if normalize:
##        fname = 'bars_{}_pts_fraction_{}.png'.format(points_type, mode)
##    else:
##        fname = 'bars_{}_pts_totals_{}.png'.format(points_type, mode)
##    plt.savefig(fname, dpi=120)
##
##plt.show()
