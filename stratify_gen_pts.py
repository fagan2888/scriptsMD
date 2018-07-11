#!/usr/bin/env python3

################################################################################
# Take track data from Boos and just get genesis points from
#  Bay of Bengal. Then stratify by Kiladis EOF's
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

def stratify_data(amp_min, amp_max, phases):
    ''' 
    Stratify the Genesis Point data based on EOF strength/phase.
      amp_min: min threshold amplitude, float
      amp_max: max threshold amplitude, float
      phases:  which phases (1-8) to look for, list
    '''
    strat_indexes = np.zeros(gen_pts.shape[0], dtype=bool)
    for i in range(len(gen_pts)):
        lon, lat, year, month, day, hour, intensity = gen_pts[i, :]
    
        for hour in [0, 12]: 
            indx = np.where( np.logical_and( np.logical_and( np.logical_and( PC_dates[:, 0] == int(year), PC_dates[:, 1] == int(month)), PC_dates[:, 2] == int(day)), PC_dates[:, 3] == int(hour)) )[0][0]

            eof1, eof2, eof3, eof4, eof5, eof6, eof7, eof8 = PC_vals[indx, :]
    
            # Biweekly mode
            x = eof1; y = eof2
            # Weekly mode
            #x = eof3; y = eof4

            amp = np.sqrt(x**2 + y**2)
            angle = get_angle_deg(x, y)
    
            if amp <= amp_max and amp >= amp_min:
                if get_phase(angle) in phases: 
                    strat_indexes[i] = True
                    break

    data = gen_pts[strat_indexes, :]
    return data

    
# Load data
geninits_NH_fname = '/home/hpeter/Documents/Research2018/MD_files/public_trackdata/geninits_NH.dat'
trx_NH_fname      = '/home/hpeter/Documents/Research2018/MD_files/public_trackdata/trx_NH.dat'
PC_dates_fname    = '/home/hpeter/Documents/Research2018/MD_files/PC_data/olr.230.westw.2x.5s30n.jja.pc_dates.dat'
PC_vals_fname     = '/home/hpeter/Documents/Research2018/MD_files/PC_data/olr.230.westw.2x.5s30n.jja.pc_vals.dat'

geninits_NH = np.loadtxt(geninits_NH_fname, dtype='int')
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

print('Number of genesis points in BoB, JJAS: {}'.format(len(gen_pts)))

# Stratify
#strong = np.zeros((8))
#weak   = np.zeros((8))
#for phase in range(1, 9):
#    # strong
#    data = stratify_data(amp_min=1.0, amp_max=np.Inf, phases=[phase])
#    n = data.shape[0]
#    strong[phase-1] = n
#    print('amp_min={}, amp_max={}, phase={}: {}'.format(1.0, np.Inf, phase, n))
#
#    ## weak
#    #data = stratify_data(amp_min=0.0, amp_max=1.0, phases=[phase])
#    #n = data.shape[0]
#    #weak[phase-1] = n
#    #print('amp_min={}, amp_max={}, phase={}: {}'.format(0.0, 1.0, phase, n))
#print('Number of genesis points after stratification: {}'.format(int(np.sum(strong[:4]))))

#data = stratify_data(amp_min=1.0, amp_max=np.Inf, phases=[1,2,3,4])
#print('Number of genesis points after stratification: {}'.format(len(data)))
#np.savez('BoB_gen_pts_JJAS_strat_biweekly_strong_p1234.npz', data)

strong = np.zeros((8))
weak   = np.zeros((8))
for i in range(PC_vals.shape[0]):
    eof1, eof2, eof3, eof4, eof5, eof6, eof7, eof8 = PC_vals[i, :]

    # Biweekly mode
    x = eof1; y = eof2
    # Weekly mode
    #x = eof3; y = eof4

    amp = np.sqrt(x**2 + y**2)
    angle = get_angle_deg(x, y)

    if amp > 1:
        strong[get_phase(angle) - 1] += 1
    else:
        weak[get_phase(angle) - 1] += 1
        

#################################################################################
f = plt.figure(figsize=(8,5))
ax = plt.subplot(111)

#ax.set_title('Density of MD Genesis w.r.t. Biweekly EOFs', size=16)
#ax.set_title('Density of MD Genesis w.r.t. Weekly EOFs', size=16)
ax.set_title('Biweekly EOFs Count', size=16)
#ax.set_title('Weekly EOFs Count', size=16)
ax.set_xlabel('Phase', size=12)
#ax.set_ylabel('# of MD Genesis Points', size=12)
ax.set_ylabel('# of Events', size=12)

centers = np.arange(1, 9)
width = 0.45
ax.bar(centers - width/2, strong, width=width, color='#aa7777', hatch='', edgecolor='k', label='strong phase')
ax.bar(centers + width/2, weak,   width=width, color='#7777aa', hatch='', edgecolor='k', label='weak phase')
ax.legend()

plt.show()

#f.savefig('biweekly_stratified.png', dpi=100)
#f.savefig('weekly_stratified.png', dpi=100)
f.savefig('biweekly_all.png', dpi=100)
#f.savefig('weekly_all.png', dpi=100)
#################################################################################
