#!/usr/bin/env python3

################################################################################
# Plot composites for Monsoon Depression genesis points 
#  in the Bay of Bengal (data from Boos)
################################################################################

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime as dt
import pandas as pd
import os
from stratify_track_data import stratify_data

def plot_vector_indiv(lvl=850, mode='biweekly', phases='1234', date=dt.datetime(year=1994, month=7, day=6, hour=18), zoom=False):
    # setup
    if zoom:
        llon = 55; rlon = 125; blat = -10; tlat = 35
    else:
        llon = 35; rlon = 145; blat = -30; tlat = 55

    # ERA data
    #fname = '{}track_UV_{}_{}{}.nc'.format(idir, date.strftime('%Y%m%d%H'), mode, phases)
    fname = '{}track_noharm_UV_{}_{}{}.nc'.format(idir, date.strftime('%Y%m%d%H'), mode, phases)
    ds = xr.open_dataset(fname)
    u_da = ds.u
    v_da = ds.v

    u_da = u_da.sel(lvl=lvl).sel({'lat' : slice(tlat, blat), 'lon' : slice(llon, rlon)})
    v_da = v_da.sel(lvl=lvl).sel({'lat' : slice(tlat, blat), 'lon' : slice(llon, rlon)})
    x = u_da.coords['lon'].values
    y = u_da.coords['lat'].values
    t = u_da.coords['time'].values
    t = pd.DatetimeIndex(t)

    # track data
    fname = '{}BoB_followed_track_{}_{}{}.npz'.format(track_dir, date.strftime('%Y%m%d'), mode, phases)
    track_pts = np.load(fname)['arr_0']

    strat = stratify_data(track_pts, PC_dates, PC_vals)

    proj = ccrs.PlateCarree()
    for i, time in enumerate(t):
        #title='{} mb Anomalous Wind\nDate: {}'.format(lvl, time.strftime('%Y-%m-%d %H:00:00')) 
        title='{} mb Wind\nDate: {}'.format(lvl, time.strftime('%Y-%m-%d %H:00:00')) 
        print(title)

        U = u_da.sel(time=time).values
        V = v_da.sel(time=time).values
        
        f = plt.figure(figsize=(8,8))
        ax = plt.subplot(111, projection=proj)

        # map
        ax.add_feature(cfeature.OCEAN, zorder=0, color='#7777cc')
        ax.add_feature(cfeature.LAND,  zorder=0, color='#77cc77')
        ax.coastlines(resolution='110m', color='black', linewidth=0.5)
        gl = ax.gridlines(crs=proj, draw_labels=True,
            linewidth=0.5, color='w', alpha=0.6, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        ax.set_xmargin(0)
        ax.set_ymargin(0)

        if zoom:
            skip = 1
        else:
            skip = 3
        if zoom:
            Q = ax.quiver(x[::skip], y[::skip], U[::skip, ::skip], V[::skip, ::skip],
                        pivot='mid', units='inches', scale=40, scale_units='width',
                        headwidth=5, headlength=3, headaxislength=2, lw=0.1)
        else:
            Q = ax.quiver(x[::skip], y[::skip], U[::skip, ::skip], V[::skip, ::skip],
                        pivot='mid', units='inches', scale=300, scale_units='width',
                        headwidth=6, headlength=4, headaxislength=3)
        qk = plt.quiverkey(Q, 0.25, 0.1, 5, '2 m/s', coordinates='figure')

        # bounding box 
        ax.plot([83, 93, 93, 83, 83], [16, 16, 21, 21, 16], 'r-', lw=1.0, transform=proj)
        ax.set_title('{}'.format(title))
        plt.tight_layout(pad=5)

        # EOFs
        ampBW = strat[i][0]
        phaseBW = strat[i][1]
        ampW = strat[i][2]
        phaseW = strat[i][3]
        ax.text(0.6, 0.13, color='k', transform=f.transFigure, size=12, s='Biweekly Phase: {:.0f}'.format(phaseBW))
        ax.text(0.6, 0.10, color='k', transform=f.transFigure, size=12, s='Biweekly Amp: {:.1f}'.format(ampBW))
        ax.text(0.6, 0.07, color='k', transform=f.transFigure, size=12, s='Weekly Phase: {:.0f}'.format(phaseW))
        ax.text(0.6, 0.04, color='k', transform=f.transFigure, size=12, s='Weekly Amp: {:.1f}'.format(ampW))

        if zoom:
            fname = '{}track_following_UV_zoomed_{}_{}{}_{}_{}.png'.format(idir, date.strftime('%Y%m%d%H'), mode, phases, lvl, str(i).zfill(3))
        else:
            #fname = '{}track_following_UV_{}_{}{}_{}_{}.png'.format(idir, date.strftime('%Y%m%d%H'), mode, phases, lvl, str(i).zfill(3))
            fname = '{}track_following_noharm_UV_{}_{}{}_{}_{}.png'.format(idir, date.strftime('%Y%m%d%H'), mode, phases, lvl, str(i).zfill(3))
        plt.savefig(fname, dpi=120)
        plt.close()

    print('\nConverting to .gif')
    command = 'convert {} {}'.format(fname[:-8] + '*', fname[:-8] + '.gif')
    print('"{}"\n'.format(command))
    os.system(command)

idir = '/home/hpeter/Research2018/MD_files/track_following/'
track_dir = '/home/hpeter/Research2018/MD_files/public_trackdata/'

# PC data for showing EOFs
PC_dates_fname    = '/home/hpeter/Documents/Research2018/MD_files/PC_data/olr.230.westw.2x.5s30n.jja.pc_dates.dat'
PC_vals_fname     = '/home/hpeter/Documents/Research2018/MD_files/PC_data/olr.230.westw.2x.5s30n.jja.pc_vals.dat'
PC_dates = np.loadtxt(PC_dates_fname, dtype='int')
PC_vals = np.loadtxt(PC_vals_fname)
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

#phases = '1234'; mode = 'biweekly'; date = dt.datetime(year=1994, month=7, day=6, hour=18)
phases = '1234'; mode = 'biweekly'; date = dt.datetime(year=2010, month=7, day=23, hour=18)
lvls   = [850, 500, 300]

for lvl in lvls:
    plot_vector_indiv(lvl=lvl, mode=mode, phases=phases, date=date, zoom=False)
