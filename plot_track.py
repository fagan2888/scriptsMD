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

def plot_vector_indiv(filename, lvl=850, mode='biweekly', 
                        phases='1234', date=dt.datetime(year=1994, month=7, day=6, hour=18), zoom=False):
    # setup
    if zoom:
        llon = 55; rlon = 125; blat = -10; tlat = 35
    else:
        llon = 35; rlon = 145; blat = -30; tlat = 55

    # data
    filename_full = idir + filename
    ds = xr.open_dataset(filename_full)
    u_da = ds.u
    v_da = ds.v

    u_da = u_da.sel(lvl=lvl).sel({'lat' : slice(tlat, blat), 'lon' : slice(llon, rlon)})
    v_da = v_da.sel(lvl=lvl).sel({'lat' : slice(tlat, blat), 'lon' : slice(llon, rlon)})
    x = u_da.coords['lon'].values
    y = u_da.coords['lat'].values
    t = u_da.coords['time'].values
    t = pd.DatetimeIndex(t)

    proj = ccrs.PlateCarree()
    for i, time in enumerate(t):
        title='{} mb Anomalous Wind\nDate: {}'.format(lvl, time.strftime('%Y-%m-%d %H:00:00')) 
        #title='{} mb Wind\nDate: {}'.format(lvl, time.strftime('%Y-%m-%d %H:00:00')) 
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
        qk = plt.quiverkey(Q, 0.5, 0.05, 5, '2 m/s', coordinates='figure')

        # bounding box 
        ax.plot([83, 93, 93, 83, 83], [16, 16, 21, 21, 16], 'r-', lw=1.0, transform=proj)
        ax.set_title('{}'.format(title))
        plt.tight_layout(pad=5)
        if zoom:
            fname = '{}track_following_UV_zoomed_{}_{}{}_{}_{}.png'.format(idir, date.strftime('%Y%m%d%H'), mode, phases, lvl, str(i).zfill(3))
        else:
            fname = '{}track_following_UV_{}_{}{}_{}_{}.png'.format(idir, date.strftime('%Y%m%d%H'), mode, phases, lvl, str(i).zfill(3))
#            fname = '{}track_following_noharm_UV_{}_{}{}_{}_{}.png'.format(idir, date.strftime('%Y%m%d%H'), mode, phases, lvl, str(i).zfill(3))
        plt.savefig(fname, dpi=120)
        plt.close()

idir = '/home/hpeter/Research2018/MD_files/track_following/'

phases = '1234'; mode = 'biweekly'; date = dt.datetime(year=1994, month=7, day=6, hour=18)
lvls   = [500, 300]

for lvl in lvls:
    plot_vector_indiv('track_UV_{}_{}{}.nc'.format(date.strftime('%Y%m%d%H'), mode, phases), 
    #plot_vector_indiv('track_noharm_UV_{}_{}{}.nc'.format(date.strftime('%Y%m%d%H'), mode, phases), 
                   lvl=lvl, mode=mode, phases=phases, date=date, zoom=False)
