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

#plt.style.use('my_mpl_stylesheet.mplstyle')

def plot1dx(filename):
    filename_full = idir + filename
    llon = 130
    rlon = 230
    f, ax = plt.subplots(1, figsize=(12, 8) )
    da = xr.open_dataarray(filename_full)
    comp = da.sel({'lv_ISBL0' : 850, 'g4_lat_1' : 0}, method='nearest').sel(g4_lon_2=slice(llon, rlon)).values
    xvals = np.linspace(llon-180, rlon-180, len(comp))
    ax.plot(xvals, comp)
    ax.plot([0,0], ax.get_ylim(), 'r-')
    ax.set_xlabel('Longitude ($^\\circ$ from genesis pt)')
    ax.set_title(filename)
    plt.tight_layout()

def plot2dxy_lag(filenames, lags=[-2,-1,0,1,2], ttest=True, lvl=850):
    # setup
    proj = ccrs.PlateCarree()
    f = plt.figure()
    filenames_full = idir + filenames
    if len(lags) == 5:
        nR = 3; nC = 2
        nums = [1, 3, 5, 2, 4]
    elif len(lags) == 3:
        nR = 2; nC = 2
        nums = [1, 3, 2, 4]
    elif len(lags) == 1:
        nR = 1; nC = 1
        nums = [1]
    #llon = 130; rlon = 230; blat = -25; tlat = 25
    llon = 130; rlon = 230; blat = -55; tlat = 30

    # loop
    for i in range(nR * nC):
        if i == len(lags):
            break
        num = str(nR) + str(nC) + str(nums[i])
        ax = plt.subplot(num, projection=proj)
        lag = lags[i]

        # map
        ax.coastlines(resolution='110m', color='black', linewidth=1)
        gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=1, color='k', alpha=0.9, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False

        # data
        if ttest:
            #da = xr.open_dataarray(filenames_full+str(lag)+'_ttest.nc')
            #da = xr.open_dataarray(filenames_full+str(lag)+'_harmonics.nc')
            da = xr.open_dataarray(filenames_full+str(lag)+'_ttest_harmonics.nc')
        else:
            da = xr.open_dataarray(filenames_full+str(lag)+'.nc')
        comp = da.sel({'lv_ISBL0' : lvl}, method='nearest').sel(
            {'g4_lat_1' : slice(tlat, blat), 'g4_lon_2' : slice(llon, rlon)}).values
        cmap = plt.get_cmap('PiYG')
        avg = np.mean(np.abs(comp))
        sd = np.std(np.abs(comp))
        im_extent = (85+180-rlon,85+180-llon,18.5+blat,18.5+tlat)
        im = ax.imshow(comp, cmap=cmap, origin='upper', extent=im_extent,
            vmin=-avg - 3*sd, vmax=avg + 3*sd, interpolation='hermite',
            transform=proj)
        cb = plt.colorbar(im, ax=ax)

        # bounding box and center
        ax.plot([75, 95, 95, 75, 75], [10, 10, 27, 27, 10], 'r-', lw=2.5, transform=proj)
        ax.plot(85, 18.5, 'ro', ms=5, transform=proj)
        if i == 0:
            if ttest:
                ax.set_title('{}: {} mb, ttest\nlag {}'.format(filenames, lvl, lag))
            else:
                ax.set_title('{}: {} mb\nlag {}'.format(filenames, lvl, lag))
        else:
            ax.set_title('lag {}'.format(lag))
    plt.tight_layout()

def plot2dxy_lag_vector(ufilenames, vfilenames, title='Default lags, ttest, 850 mb', lags=[-2,-1,0,1,2], ttest=True, lvl=850):
    # setup
    proj = ccrs.PlateCarree()
    #f = plt.figure(figsize=(16,5))
    f = plt.figure(figsize=(10,10))
    ufilenames_full = idir + ufilenames
    vfilenames_full = idir + vfilenames
    if len(lags) == 5:
        nR = 1; nC = 5
        nums = [1, 2, 3, 4, 5] 
    elif len(lags) == 3:
        nR = 2; nC = 2
        nums = [1, 3, 2, 4]
    elif len(lags) == 1:
        nR = 1; nC = 1
        nums = [1]
    llon = 150; rlon = 210; blat = -25; tlat = 25
    #llon = 130; rlon = 230; blat = -55; tlat = 30

    # loop
    for i in range(nR * nC):
        if i == len(lags):
            break
        num = str(nR) + str(nC) + str(nums[i])
        ax = plt.subplot(num, projection=proj)
        lag = lags[i]

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

        # data
        if ttest:
            #u_da_raw = xr.open_dataarray(ufilenames_full+str(lag)+'.nc')
            #v_da_raw = xr.open_dataarray(vfilenames_full+str(lag)+'.nc')
            #u_da_sig = xr.open_dataarray(ufilenames_full+str(lag)+'_ttest.nc')
            #v_da_sig = xr.open_dataarray(vfilenames_full+str(lag)+'_ttest.nc')
            u_da_raw = xr.open_dataarray(ufilenames_full+str(lag)+'_harmonics.nc')
            v_da_raw = xr.open_dataarray(vfilenames_full+str(lag)+'_harmonics.nc')
            u_da_sig = xr.open_dataarray(ufilenames_full+str(lag)+'_ttest_harmonics.nc')
            v_da_sig = xr.open_dataarray(vfilenames_full+str(lag)+'_ttest_harmonics.nc')
            Uraw = u_da_raw.values
            Vraw = v_da_raw.values
            Usig = u_da_sig.values
            Vsig = v_da_sig.values
            # if either component is significant, draw
            Usig[np.where(Vsig != 0.0)] = Uraw[np.where(Vsig != 0.0)]
            Vsig[np.where(Usig != 0.0)] = Vraw[np.where(Usig != 0.0)]
            # else kill
            Usig[np.where(Usig == 0.0)] = np.nan
            Vsig[np.where(Vsig == 0.0)] = np.nan
            u_da = u_da_raw
            v_da = v_da_raw
            u_da.values = Usig
            v_da.values = Vsig
        else:
            u_da = xr.open_dataarray(ufilenames_full+str(lag)+'.nc')
            v_da = xr.open_dataarray(vfilenames_full+str(lag)+'.nc')

        u_comp = u_da.sel({'lv_ISBL0' : lvl}, method='nearest').sel(
            {'g4_lat_1' : slice(tlat, blat), 'g4_lon_2' : slice(llon, rlon)})
        v_comp = v_da.sel({'lv_ISBL0' : lvl}, method='nearest').sel(
            {'g4_lat_1' : slice(tlat, blat), 'g4_lon_2' : slice(llon, rlon)})
        x = u_comp.coords['g4_lon_2'].values - 180 + 85
        y = u_comp.coords['g4_lat_1'].values + 18.5
        U = u_comp.values; V = v_comp.values
        
        skip = 1
        #skip = 3
        Q = ax.quiver(x[::skip], y[::skip], U[::skip, ::skip], V[::skip, ::skip],
                    #pivot='mid', units='inches', scale=25, scale_units='width',
                    #headwidth=10, headlength=5, headaxislength=3)
                    pivot='mid', units='inches', scale=50, scale_units='width',
                    headwidth=5, headlength=2.5, headaxislength=1.5, lw=0.5)
        #qk = plt.quiverkey(Q, 0.5, 0.15, 2, '2 m/s', coordinates='figure')
        qk = plt.quiverkey(Q, 0.95, 0.95, 2, '2 m/s', coordinates='figure')

        # bounding box and center
        ax.plot([75, 95, 95, 75, 75], [10, 10, 27, 27, 10], 'r-', lw=1.0, transform=proj)
        ax.plot(85, 18.5, 'ro', ms=3, transform=proj)
        #if i == 2:
        if i == 0:
            ax.set_title('{}\nlag ${}$'.format(title, lag))
        else:
            ax.set_title('lag ${}$'.format(lag))
    #plt.tight_layout(pad=3, w_pad=-4)
    plt.tight_layout(pad=3, w_pad=-3)
    #fname = '../images/winds' + str(lvl) + '_harmonics.png'
    fname = '../images/winds' + str(lvl) + '_harmonics_zoomed.png'
    plt.savefig(fname, dpi=300)

def plot2dxz_lag(filenames, title='Default lags, ttest', lags=[-2,-1,0,1,2], ttest=True):
    # setup
    f = plt.figure(figsize=(16,5))
    filenames_full = idir + filenames
    if len(lags) == 5:
        nR = 1; nC = 5
        nums = [1, 2, 3, 4, 5]
    elif len(lags) == 3:
        nR = 2; nC = 2
        nums = [1, 3, 2, 4]
    elif len(lags) == 1:
        nR = 1; nC = 1
        nums = [1]
    llon = 130; rlon = 230; llvl = 1000; ulvl = 1

    # loop
    for i in range(nR * nC):
        if i == len(lags):
            break
        num = str(nR) + str(nC) + str(nums[i])
        ax = plt.subplot(num)
        lag = lags[i]
        if ttest:
            da = xr.open_dataarray(filenames_full+str(lag)+'_ttest.nc')
        else:
            da = xr.open_dataarray(filenames_full+str(lag)+'.nc')
        comp = da.sel({'g4_lat_1' : 0}, method='nearest').sel(
            {'g4_lon_2' : slice(llon, rlon), 'lv_ISBL0' : slice(ulvl, llvl)}).values
        cmap = plt.get_cmap('viridis')
        avg = np.mean(np.abs(comp))
        sd = np.std(np.abs(comp))
        im = ax.imshow(comp, cmap=cmap, extent=(llon-180, rlon-180, llvl, ulvl), aspect=0.08,
            vmin=-1.5, vmax=1.5, interpolation='hermite')
            #vmin=-avg - 3*sd, vmax=avg + 3*sd, interpolation='hermite')
        #levels = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5]
        levels = [-1.5, -0.75, 0.0, 0.75, 1.5]
        cs = ax.contour(comp, extent=(llon-180, rlon-180, llvl, ulvl), origin='upper',
            levels=levels, colors='w', linewidths=1.0)

        if i in [4]:
            cb = plt.colorbar(im, ax=ax, label='m/s', fraction=0.04, pad=0.04)
            cb.add_lines(cs)

        if i in [0]:
            ax.set_ylabel('Pressure (hPa)')

        if i in [2]:
            ax.set_title('{}\nlag {}'.format(title, lag))
            ax.set_xlabel('Longitude (degrees from genesis pt)')
        else:
            ax.set_title('lag {}'.format(lag))
    plt.tight_layout(pad=1, w_pad=-5)
    fname = 'verticalV_ttest.png'
    #plt.savefig(fname, dpi=300)

lags0 = [0]
lags02 = [-2, 0, 2]
lags012 = [-2, -1, 0, 1, 2]

idir = '/global/scratch/hpeter/composites/'

#plot1dx('composite_n150_Vanom_JJAS_lag0.nc')

#plot2dxy_lag('composite_n150_Vanom_JJAS_lag', lags=lags0, ttest=True, lvl=850)
#plot2dxy_lag('composite_n150_Vanom_JJAS_lag', lags=lags02, ttest=True, lvl=300)
#plot2dxy_lag('composite_n150_Uanom_JJAS_lag', lags=lags02, ttest=True, lvl=850)
#plot2dxy_lag('composite_n150_Uanom_JJAS_lag', lags=lags02, ttest=True, lvl=300)

#plot2dxy_lag_vector('composite_n150_Uanom_JJAS_lag', 'composite_n150_Vanom_JJAS_lag', 
#    title='Significant Anomalous Wind (850 mb, n = 150)', lags=lags012, ttest=True, lvl=850)
#plot2dxy_lag_vector('composite_n150_Uanom_JJAS_lag', 'composite_n150_Vanom_JJAS_lag', 
#    title='Significant Anomalous Wind (300 mb, n = 150)', lags=lags012, ttest=True, lvl=300)
#
#plot2dxz_lag('composite_n150_Vanom_JJAS_lag', 
#    title='Anomalous V (n = 150, genesis latitude)', lags=lags012, ttest=False)
#plot2dxz_lag('composite_n150_Vanom_JJAS_lag', 
#    title='Significant Anomalous V (n = 150, genesis latitude)', lags=lags012, ttest=True)

#plot2dxy_lag_vector('composite_n150_Uanom_JJAS_lag', 'composite_n150_Vanom_JJAS_lag', 
#    title='Significant Anomalous Wind (850 mb, n = 150)', lags=lags0, ttest=True, lvl=850)
#plot2dxy_lag_vector('composite_n150_Uanom_JJAS_lag', 'composite_n150_Vanom_JJAS_lag', 
#    title='Significant Anomalous Wind (300 mb, n = 150)', lags=lags0, ttest=True, lvl=300)

#plot2dxy_lag('composite_n150_U_JJAS_lag', lags=lags012, ttest=True, lvl=850)
#plot2dxy_lag('composite_n150_U_JJAS_lag', lags=lags012, ttest=True, lvl=300)
#plot2dxy_lag('composite_n150_V_JJAS_lag', lags=lags012, ttest=True, lvl=850)
#plot2dxy_lag('composite_n150_V_JJAS_lag', lags=lags012, ttest=True, lvl=300)

plot2dxy_lag_vector('composite_n150_U_JJAS_lag', 'composite_n150_V_JJAS_lag', 
    title='Significant Anomalous Wind (850 mb, n = 150)', lags=lags0, ttest=True, lvl=850)
plot2dxy_lag_vector('composite_n150_U_JJAS_lag', 'composite_n150_V_JJAS_lag', 
    title='Significant Anomalous Wind (300 mb, n = 150)', lags=lags0, ttest=True, lvl=300)
#plot2dxy_lag_vector('composite_n150_U_JJAS_lag', 'composite_n150_V_JJAS_lag', 
#    title='Significant Anomalous Wind (850 mb, n = 150)', lags=lags012, ttest=True, lvl=850)
#plot2dxy_lag_vector('composite_n150_U_JJAS_lag', 'composite_n150_V_JJAS_lag', 
#    title='Significant Anomalous Wind (300 mb, n = 150)', lags=lags012, ttest=True, lvl=300)

plt.show()
