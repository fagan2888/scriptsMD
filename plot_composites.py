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
    comp = da.sel({'lvl' : 850, 'lat' : 0}, method='nearest').sel(lon=slice(llon, rlon)).values
    xvals = np.linspace(llon-180, rlon-180, len(comp))
    ax.plot(xvals, comp)
    ax.plot([0,0], ax.get_ylim(), 'r-')
    ax.set_xlabel('Longitude ($^\\circ$ from genesis pt)')
    ax.set_title(filename)
    plt.tight_layout()

def plot2dxy_lag(filenames, lags=[-2,-1,0,1,2], lvl=850):
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
    llon = 35; rlon = 145; blat = -30; tlat = 55

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
        da = xr.open_dataarray(filenames_full+str(lag)+'.nc')
        comp = da.sel({'lvl' : lvl}, method='nearest').sel(
            {'lat' : slice(tlat, blat), 'lon' : slice(llon, rlon)}).values
        cmap = plt.get_cmap('PiYG')
        avg = np.mean(np.abs(comp))
        sd = np.std(np.abs(comp))
        im_extent = (llon,rlon,blat,tlat)
        im = ax.imshow(comp, cmap=cmap, origin='upper', extent=im_extent,
            vmin=-avg - 3*sd, vmax=avg + 3*sd, interpolation='hermite',
            transform=proj)
        cb = plt.colorbar(im, ax=ax)

        # bounding box and center
        ax.plot([75, 95, 95, 75, 75], [10, 10, 27, 27, 10], 'r-', lw=2.5, transform=proj)
        ax.plot(85, 18.5, 'ro', ms=5, transform=proj)
        if i == 0:
            ax.set_title('{}: {} mb\nlag ${}$'.format(filenames, lvl, lag))
        else:
            ax.set_title('lag ${}$'.format(lag))
    plt.tight_layout()

def plot2dxy_lag_vector_indiv(ufilenames, vfilenames,
                        ufilenames_ttest=None, 
                        vfilenames_ttest=None, 
                        lags=[-2,-1,0,1,2], 
                        lvl=850, mode='biweekly', 
                        phases='1234', zoom=False):
    # setup
    title='Significant Anomalous Wind\n({} mb, n = {}, strong {} mode phases {})'.format(lvl, n, mode, phases) 
    proj = ccrs.PlateCarree()
    if zoom:
        llon = 55; rlon = 125; blat = -10; tlat = 35
    else:
        llon = 35; rlon = 145; blat = -30; tlat = 55

    # loop
    for i in range(len(lags)):
        f = plt.figure(figsize=(8,8))
        ax = plt.subplot(111, projection=proj)
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
        ufilenames_full = idir + ufilenames
        vfilenames_full = idir + vfilenames
        if ufilenames_ttest:
            ufilenames_ttest_full = idir + ufilenames_ttest
            vfilenames_ttest_full = idir + vfilenames_ttest
            u_da_raw = xr.open_dataarray(ufilenames_full+str(lag)+'.nc')
            v_da_raw = xr.open_dataarray(vfilenames_full+str(lag)+'.nc')
            u_da_sig = xr.open_dataarray(ufilenames_ttest_full+str(lag)+'.nc')
            v_da_sig = xr.open_dataarray(vfilenames_ttest_full+str(lag)+'.nc')
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

        u_comp = u_da.sel({'lvl' : lvl}, method='nearest').sel(
            {'lat' : slice(tlat, blat), 'lon' : slice(llon, rlon)})
        v_comp = v_da.sel({'lvl' : lvl}, method='nearest').sel(
            {'lat' : slice(tlat, blat), 'lon' : slice(llon, rlon)})
        x = u_comp.coords['lon'].values
        y = u_comp.coords['lat'].values
        U = u_comp.values; V = v_comp.values
        
        if zoom:
            skip = 1
        else:
            skip = 2
        if zoom:
            Q = ax.quiver(x[::skip], y[::skip], U[::skip, ::skip], V[::skip, ::skip],
                        pivot='mid', units='inches', scale=40, scale_units='width',
                        headwidth=5, headlength=3, headaxislength=2, lw=0.1)
        else:
            Q = ax.quiver(x[::skip], y[::skip], U[::skip, ::skip], V[::skip, ::skip],
                        pivot='mid', units='inches', scale=25, scale_units='width',
                        headwidth=6, headlength=4, headaxislength=3)
        qk = plt.quiverkey(Q, 0.5, 0.05, 2, '2 m/s', coordinates='figure')

        # bounding box and center
        ax.plot([83, 93, 93, 83, 83], [16, 16, 21, 21, 16], 'r-', lw=1.0, transform=proj)
        ax.plot((83+93)/2, (16+21)/2, 'ro', ms=3, transform=proj)
        ax.set_title('{}\nlag ${}$'.format(title, lag))
        plt.tight_layout(pad=5)
        #fname = '/global/home/users/hpeter/images/winds' + str(lvl) + '_biweekly1234s.png'
        if zoom:
            fname = '/home/hpeter/Research2018/images/winds{}mb_{}{}_zoomed_lag{}.png'.format(lvl, mode, phases, lag)
        else:
            # fname = 'genesis_{}mb_n{}_int{}_{}{}_vort850_{}.png'.format(lvl, n, intensities, mode, phases, lag+2)
            fname = 'tracks_{}mb_n{}_int{}_{}{}_vort850_{}.png'.format(lvl, n, intensities, mode, phases, lag+2)
        print(fname)
        plt.savefig(fname, dpi=100)

def plot2dxy_lag_vector(ufilenames, vfilenames,
                        ufilenames_ttest=None, 
                        vfilenames_ttest=None, 
                        lags=[-2,-1,0,1,2], 
                        lvl=850, mode='biweekly', 
                        phases='1234'): 
    # setup
    title='Significant Anomalous Wind\n({} mb, n = {}, {} strong mode phases {})'.format(lvl, n, mode, phases) 
    proj = ccrs.PlateCarree()
    f = plt.figure(figsize=(16,5))
    if len(lags) == 5:
        nR = 1; nC = 5
        nums = [1, 2, 3, 4, 5] 
    elif len(lags) == 3:
        nR = 2; nC = 2
        nums = [1, 3, 2, 4]
    elif len(lags) == 1:
        nR = 1; nC = 1
        nums = [1]
    llon = 35; rlon = 145; blat = -30; tlat = 55

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
        ufilenames_full = idir + ufilenames
        vfilenames_full = idir + vfilenames
        if ufilenames_ttest:
            ufilenames_ttest_full = idir + ufilenames_ttest
            vfilenames_ttest_full = idir + vfilenames_ttest
            u_da_raw = xr.open_dataarray(ufilenames_full+str(lag)+'.nc')
            v_da_raw = xr.open_dataarray(vfilenames_full+str(lag)+'.nc')
            u_da_sig = xr.open_dataarray(ufilenames_ttest_full+str(lag)+'.nc')
            v_da_sig = xr.open_dataarray(vfilenames_ttest_full+str(lag)+'.nc')
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

        u_comp = u_da.sel({'lvl' : lvl}, method='nearest').sel(
            {'lat' : slice(tlat, blat), 'lon' : slice(llon, rlon)})
        v_comp = v_da.sel({'lvl' : lvl}, method='nearest').sel(
            {'lat' : slice(tlat, blat), 'lon' : slice(llon, rlon)})
        x = u_comp.coords['lon'].values
        y = u_comp.coords['lat'].values
        U = u_comp.values; V = v_comp.values
        
        skip = 3
        Q = ax.quiver(x[::skip], y[::skip], U[::skip, ::skip], V[::skip, ::skip],
                    pivot='mid', units='inches', scale=25, scale_units='width',
                    headwidth=10, headlength=5, headaxislength=4)
        qk = plt.quiverkey(Q, 0.5, 0.15, 2, '2 m/s', coordinates='figure')

        # bounding box and center
        #ax.plot([75, 95, 95, 75, 75], [10, 10, 27, 27, 10], 'r-', lw=1.0, transform=proj)
        #ax.plot(85, 18.5, 'ro', ms=3, transform=proj)
        ax.plot([83, 93, 93, 83, 83], [16, 16, 21, 21, 16], 'r-', lw=1.0, transform=proj)
        ax.plot((83+93)/2, (16+21)/2, 'ro', ms=3, transform=proj)
        if i == 2:
            ax.set_title('{}\nlag ${}$'.format(title, lag))
        else:
            ax.set_title('lag ${}$'.format(lag))
    plt.tight_layout(pad=3, w_pad=-6)
    #fname = '/home/hpeter/Research2018/images/call20180712/winds{}mb_{}{}.png'.format(lvl, mode, phases)
    fname = 'winds{}mb_alltracks_int23_{}{}_olr.png'.format(lvl, mode, phases)
    plt.savefig(fname, dpi=120)

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
        comp = da.sel({'lat' : 0}, method='nearest').sel(
            {'lon' : slice(llon, rlon), 'lvl' : slice(ulvl, llvl)}).values
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
    #plt.savefig(fname, dpi=120)

def plot2dxy_vector_EOF(ufilenames, vfilenames,
                        ufilenames_ttest=None, 
                        vfilenames_ttest=None, 
                        lvl=850, mode='biweekly', 
                        phases='1234', 
                        amp='strong', zoom=False):
    # setup
    title='Significant Anomalous Wind\n({} mb, n = {}, {} {} mode phase {})'.format(lvl, n, amp, mode, phases) 
    proj = ccrs.PlateCarree()
    if zoom:
        llon = 20; rlon = 170; blat = -15; tlat = 40
    else:
        llon = 35; rlon = 145; blat = -30; tlat = 55

    if zoom:
        f = plt.figure(figsize=(12,5))
    else:
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

    # data
    ufilenames_full = idir + ufilenames
    vfilenames_full = idir + vfilenames
    if ufilenames_ttest:
        ufilenames_ttest_full = idir + ufilenames_ttest
        vfilenames_ttest_full = idir + vfilenames_ttest
        u_da_raw = xr.open_dataarray(ufilenames_full)
        v_da_raw = xr.open_dataarray(vfilenames_full)
        u_da_sig = xr.open_dataarray(ufilenames_ttest_full)
        v_da_sig = xr.open_dataarray(vfilenames_ttest_full)
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
        u_da = xr.open_dataarray(ufilenames_full)
        v_da = xr.open_dataarray(vfilenames_full)

    u_comp = u_da.sel({'lvl' : lvl}, method='nearest').sel(
        {'lat' : slice(tlat, blat), 'lon' : slice(llon, rlon)})
    v_comp = v_da.sel({'lvl' : lvl}, method='nearest').sel(
        {'lat' : slice(tlat, blat), 'lon' : slice(llon, rlon)})
    x = u_comp.coords['lon'].values 
    y = u_comp.coords['lat'].values 
    U = u_comp.values; V = v_comp.values
    
    if zoom:
        skip = 2
    else:
        skip = 2
    if zoom:
        Q = ax.quiver(x[::skip], y[::skip], U[::skip, ::skip], V[::skip, ::skip],
                    pivot='mid', units='inches', scale=80, scale_units='width',
                    headwidth=6, headlength=4, headaxislength=2.5, width=0.01)
    else:
        Q = ax.quiver(x[::skip], y[::skip], U[::skip, ::skip], V[::skip, ::skip],
                    pivot='mid', units='inches', scale=25, scale_units='width',
                    headwidth=6, headlength=4, headaxislength=3)
    qk = plt.quiverkey(Q, 0.5, 0.05, 2, '2 m/s', coordinates='figure')

    ax.set_title(title)
    plt.tight_layout(pad=5)
    if zoom:
        fname = 'EOFs_vort850_UV_{}mb_{}{}.png'.format(lvl, mode, phases)
    else:
        fname = 'EOFs_vort850_UV_{}mb_{}{}.png'.format(lvl, mode, phases)
    plt.savefig(fname, dpi=120)
    plt.close()

lags0 = [0]
lags02 = [-2, 0, 2]
lags012 = [-2, -1, 0, 1, 2]

#idir = '/global/scratch/hpeter/composites/'
idir = '/home/hpeter/Research2018/MD_files/composites/'

#plot1dx('composite_n150_Vanom_JJAS_lag0.nc')

#plot2dxy_lag('composite_n150_V_JJAS_strat_harm_ttest_lag', lags=lags012, lvl=850)
#plot2dxy_lag('composite_n150_V_JJAS_strat_harm_ttest_lag', lags=lags012, lvl=300)


# n = 97; intensities = '123'; mode ='weekly'; phases ='1678'
# n = 47; intensities = '123'; mode ='weekly'; phases ='2345'

# n = 500; intensities = '123'; mode ='weekly'; phases ='1238'
n = 500; intensities = '123'; mode ='weekly'; phases ='4567'
lvls = [850, 500]
for lvl in lvls:
    plot2dxy_lag_vector_indiv('composite_n{}_int{}_U_{}{}_vort850_lag'.format(n, intensities, mode, phases), 
                              'composite_n{}_int{}_V_{}{}_vort850_lag'.format(n, intensities, mode, phases), 
                              'composite_n{}_int{}_U_{}{}_vort850_ttest_lag'.format(n, intensities, mode, phases),
                              'composite_n{}_int{}_V_{}{}_vort850_ttest_lag'.format(n, intensities, mode, phases),
                               lags=lags012, lvl=lvl, mode=mode, phases=phases, zoom=False)

## genesis points
##phases = '1678'; mode = 'biweekly'; n = '91'
##phases = '1238'; mode = 'weekly'; n = '87'

## all tracks
## phases = '1234'; mode = 'biweekly'; n = '2359'
## phases = '1234'; mode = 'biweekly'; n = '1129'

## phases = '1278'; mode = 'weekly'; n = '2383'
## phases = '1278'; mode = 'weekly'; n = '1222'

## phases = '5678'; mode = 'biweekly'; n = '1465'
## phases = '3456'; mode = 'weekly'; n = '1416'

#lvls = [850]
#for lvl in lvls:
#    #plot2dxy_lag_vector_indiv('composite_n{}_U_{}{}_lag'.format(n, mode, phases), 
#    #                    'composite_n{}_V_{}{}_lag'.format(n, mode, phases), 
#    #              'composite_n{}_U_{}{}_ttest_lag'.format(n, mode, phases),
#    #              'composite_n{}_V_{}{}_ttest_lag'.format(n, mode, phases),
#    #               lags=lags012, lvl=lvl, mode=mode, phases=phases, zoom=True)
#    # plot2dxy_lag_vector('composite_n{}_U_{}{}_lag'.format(n, mode, phases), 
#    #                     'composite_n{}_V_{}{}_lag'.format(n, mode, phases), 
#    #               'composite_n{}_U_{}{}_ttest_lag'.format(n, mode, phases),
#    #               'composite_n{}_V_{}{}_ttest_lag'.format(n, mode, phases),
#    #                lags=lags012, lvl=lvl, mode=mode, phases=phases)
#    plot2dxy_lag_vector('composite_n{}_U_int3_{}{}_lag'.format(n, mode, phases), 
#                        'composite_n{}_V_int3_{}{}_lag'.format(n, mode, phases), 
#                  'composite_n{}_U_int3_{}{}_ttest_lag'.format(n, mode, phases),
#                  'composite_n{}_V_int3_{}{}_ttest_lag'.format(n, mode, phases),
#                   lags=lags012, lvl=lvl, mode=mode, phases=phases)

# all EOFs
#phases = '1678'; mode = 'biweekly'; n = '4201'
#phases = '1238'; mode = 'weekly'; n = '4271'
#phases = '1278'; mode = 'weekly'; n = '4271'
#phases = '1234'; mode = 'biweekly'; n = '4318'

#mode = 'weekly'
#lvls = [850]
#ns = ['2303', '2198', '1939', '1832', '1843','2108', '2274', '2370']
#ps = ['1', '2', '3', '4', '5', '6', '7', '8']
#for n, phases in zip(ns, ps):
#    for lvl in lvls:
#        plot2dxy_vector_EOF('composite_EOFs_n{}_U_vort850_{}{}.nc'.format(n, mode, phases), 
#                            'composite_EOFs_n{}_V_vort850_{}{}.nc'.format(n, mode, phases), 
#                      'composite_EOFs_n{}_U_vort850_{}{}_ttest.nc'.format(n, mode, phases),
#                      'composite_EOFs_n{}_V_vort850_{}{}_ttest.nc'.format(n, mode, phases),
#                      lvl=lvl, mode=mode, phases=phases, zoom=True)

# plt.show()
