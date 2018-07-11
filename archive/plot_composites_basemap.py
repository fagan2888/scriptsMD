#!/usr/bin/env python3

################################################################################
# Plot composites for Monsoon Depression genesis points 
#  in the Bay of Bengal (data from Boos)
################################################################################

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import warnings
warnings.filterwarnings('ignore') #basemap 'hold' depreciation

#plt.style.use('my_mpl_stylesheet.mplstyle')

def plot1dx(da, title):
    llon = 130
    rlon = 230
    f, ax = plt.subplots(1, figsize=(12, 8) )
    comp = da.sel({'lv_ISBL0' : 850, 'g4_lat_1' : 0}, method='nearest').sel(g4_lon_2=slice(llon, rlon)).values
    xvals = np.linspace(llon-180, rlon-180, len(comp))
    ax.plot(xvals, comp)
    ax.plot([0,0], ax.get_ylim(), 'r-')
    ax.set_xlabel('Longitude ($^\\circ$ from genesis pt)')
    ax.set_title(title)
    plt.tight_layout(pad=4)

def plot2dxy(filename, lvl):
    f, ax = plt.subplots(1, figsize=(12, 8) )
    #llon = 130; rlon = 230; blat = -25; tlat = 25
    llon = 130; rlon = 230; blat = -55; tlat = 30
    #draw map
    m = Basemap(ax=ax, llcrnrlon=85+180-rlon,llcrnrlat=18.5+blat,urcrnrlon=85+180-llon,urcrnrlat=18.5+tlat, 
        resolution='c', projection='merc', lat_0 = 18.5, lon_0 = 85)
    m.drawcoastlines()
    parallels = np.arange(-90,90,15.)
    m.drawparallels(parallels,labels=[1,0,0,0])
    meridians = np.arange(0.,360.,15.)
    m.drawmeridians(meridians,labels=[0,0,0,1])
    #imshow on top
    da = xr.open_dataarray(filename)
    comp = da.sel({'lv_ISBL0' : lvl}, method='nearest').sel({'g4_lat_1' : slice(tlat, blat), 'g4_lon_2' : slice(llon, rlon)}).values
    cmap = plt.get_cmap('PiYG')
    avg = np.mean(np.abs(comp))
    sd = np.std(np.abs(comp))
    im = m.imshow(comp, cmap=cmap, origin='upper',
        vmin=-avg - 3*sd, vmax=avg + 3*sd, interpolation='hermite')
    cb = plt.colorbar(im, ax=ax)
    #bounding box and center
    x, y = m([75, 95, 95, 75, 75], [10, 10, 27, 27, 10])
    m.plot(x, y, 'r-', lw=2.5)
    x, y = m([85], [18.5])
    m.plot(x, y, 'ro', ms=5)
    ax.set_title('{}: {} mb'.format(filename, lvl))
    plt.tight_layout(pad=4)
    
def plot2dxz(da, title):
    f, ax = plt.subplots(1, figsize=(12, 8) )
    llon = 130
    rlon = 230
    llvl = 1000
    ulvl = 1
    comp = da.sel({'g4_lat_1' : 0}, method='nearest').sel({'g4_lon_2' : slice(llon, rlon), 'lv_ISBL0' : slice(ulvl, llvl)}).values
    cmap = plt.get_cmap('PiYG')
    avg = np.mean(np.abs(comp))
    sd = np.std(np.abs(comp))
    im = ax.imshow(comp, cmap=cmap, extent=(llon-180, rlon-180, llvl, ulvl), aspect=0.05,
        vmin=-avg - 3*sd, vmax=avg + 3*sd, interpolation='hermite')
    cb = plt.colorbar(im, ax=ax)
    ax.set_xlabel('Longitude ($^\\circ$ from genesis pt)')
    ax.set_title(title)
    ax.margins(0, 0)
    plt.tight_layout(pad=4)

def plot2dxy_lag(filenames, lags=[-2,-1,0,1,2], ttest=True, lvl=850):
    if len(lags) == 5:
        nR = 3; nC = 2
        f, axes = plt.subplots(nR, nC, figsize=(16, 10) )
    elif len(lags) == 3:
        nR = 2; nC = 2
        f, axes = plt.subplots(nR, nC, figsize=(16, 10) )
    #llon = 130; rlon = 230; blat = -25; tlat = 25
    llon = 130; rlon = 230; blat = -55; tlat = 30
    lag_i = 0
    for j in range(nC):
        for i in range(nR):
            if nC > 1 and nR > 1:
                ax = axes[i, j]
            elif nC == 1:
                ax = axes[i]
            else:
                ax = axes[j]
            if lag_i == len(lags):
                f.delaxes(ax)
                plt.draw()
                break
            lag = lags[lag_i]
            #draw map
            m = Basemap(ax=ax, llcrnrlon=85+180-rlon,llcrnrlat=18.5+blat,urcrnrlon=85+180-llon,urcrnrlat=18.5+tlat, 
                resolution='c', projection='merc', lat_0 = 18.5, lon_0 = 85)
            m.drawcoastlines()
            parallels = np.arange(-90,90,15.)
            m.drawparallels(parallels,labels=[1,0,0,0])
            meridians = np.arange(0.,360.,15.)
            m.drawmeridians(meridians,labels=[0,0,0,1])
            #imshow on top
            if ttest:
                da = xr.open_dataarray(filenames+str(lag)+'_ttest.nc')
            else:
                da = xr.open_dataarray(filenames+str(lag)+'.nc')
            comp = da.sel({'lv_ISBL0' : lvl}, method='nearest').sel(
                {'g4_lat_1' : slice(tlat, blat), 'g4_lon_2' : slice(llon, rlon)}).values
            cmap = plt.get_cmap('PiYG')
            avg = np.mean(np.abs(comp))
            sd = np.std(np.abs(comp))
            im = m.imshow(comp, cmap=cmap, origin='upper', 
                vmin=-avg - 3*sd, vmax=avg + 3*sd, interpolation='hermite')
            cb = plt.colorbar(im, ax=ax)
            #bounding box and center
            x, y = m([75, 95, 95, 75, 75], [10, 10, 27, 27, 10])
            m.plot(x, y, 'r-', lw=2.5)
            x, y = m([85], [18.5])
            m.plot(x, y, 'ro', ms=5)
            if i == 0 and j == 0:
                if ttest:
                    ax.set_title('{}: {} mb, ttest\nlag {}'.format(filenames, lvl, lag))
                else:
                    ax.set_title('{}: {} mb\nlag {}'.format(filenames, lvl, lag))
            else:
                ax.set_title('lag {}'.format(lag))
            lag_i += 1
    plt.tight_layout(pad=4)

def plot2dxy_vector(ufilename, vfilename, ttest=True, lvl=850):
    f, ax = plt.subplots(1, figsize=(16, 10) )
    #llon = 130; rlon = 230; blat = -25; tlat = 25
    llon = 130; rlon = 230; blat = -55; tlat = 30
    #draw map
    m = Basemap(ax=ax, llcrnrlon=85+(llon-180),llcrnrlat=18.5+blat,urcrnrlon=85+(rlon-180),urcrnrlat=18.5+tlat, 
        resolution='c', projection='merc', lat_0 = 18.5, lon_0 = 85)
    m.drawcoastlines()
    parallels = np.arange(-90,90,15.)
    m.drawparallels(parallels,labels=[1,0,0,0])
    meridians = np.arange(0.,360.,15.)
    m.drawmeridians(meridians,labels=[0,0,0,1])
    #quiver on top
    if ttest:
        u_da_raw = xr.open_dataarray(ufilename+'.nc')
        v_da_raw = xr.open_dataarray(vfilename+'.nc')
        u_da_sig = xr.open_dataarray(ufilename+'_ttest.nc')
        v_da_sig = xr.open_dataarray(vfilename+'_ttest.nc')
        Uraw = u_da_raw.values
        Vraw = v_da_raw.values
        Usig = u_da_sig.values
        Vsig = v_da_sig.values
        Usig[np.where(Vsig != 0.0)] = Uraw[np.where(Vsig != 0.0)]
        Vsig[np.where(Usig != 0.0)] = Vraw[np.where(Usig != 0.0)]
        u_da = u_da_raw
        v_da = v_da_raw
        u_da.values = Usig
        v_da.values = Vsig
    else:
        u_da = xr.open_dataarray(ufilename+'.nc')
        v_da = xr.open_dataarray(vfilename+'.nc')

    u_comp = u_da.sel({'lv_ISBL0' : lvl}, method='nearest').sel(
        {'g4_lat_1' : slice(tlat, blat), 'g4_lon_2' : slice(llon, rlon)})
    v_comp = v_da.sel({'lv_ISBL0' : lvl}, method='nearest').sel(
        {'g4_lat_1' : slice(tlat, blat), 'g4_lon_2' : slice(llon, rlon)})
    x = u_comp.coords['g4_lon_2'].values - 180 + 85
    y = u_comp.coords['g4_lat_1'].values + 18.5
    x, dummy = m(x, x)
    dummy, y = m(y, y)
    U = u_comp.values; V = v_comp.values

    #Q = m.quiver(x, y, U, V)
    Q = m.quiver(x[::3], y[::3], U[::3, ::3], V[::3, ::3],
            pivot='mid', units='inches', scale=75, scale_units='width',
            headwidth=5, headlength=2, headaxislength=5)
    qk = plt.quiverkey(Q, 0.9, 0.9, 2, '2 m/s', labelpos='E',
                       coordinates='figure')
    #bounding box and center
    x, y = m([75, 95, 95, 75, 75], [10, 10, 27, 27, 10])
    m.plot(x, y, 'r-', lw=2.5)
    x, y = m([85], [18.5])
    m.plot(x, y, 'ro', ms=5)
    if ttest:
        ax.set_title('{}: {} mb, ttest'.format(ufilename, lvl))
    else:
        ax.set_title('{}: {} mb'.format(ufilename, lvl))
    plt.tight_layout(pad=4)

def plot2dxy_lag_vector(ufilenames, vfilenames, lags=[-2,-1,0,1,2], ttest=True, lvl=850):
    if len(lags) == 5:
        nR = 3; nC = 2
        f, axes = plt.subplots(nR, nC, figsize=(16, 10) )
    elif len(lags) == 3:
        nR = 2; nC = 2
        f, axes = plt.subplots(nR, nC, figsize=(16, 10) )
    llon = 130; rlon = 230; blat = -25; tlat = 25
    #llon = 130; rlon = 230; blat = -55; tlat = 30
    lag_i = 0
    for j in range(nC):
        for i in range(nR):
            if nC > 1 and nR > 1:
                ax = axes[i, j]
            elif nC == 1:
                ax = axes[i]
            else:
                ax = axes[j]
            if lag_i == len(lags):
                f.delaxes(ax)
                plt.draw()
                break
            lag = lags[lag_i]
            #draw map
            m = Basemap(ax=ax, llcrnrlon=85+(llon-180),llcrnrlat=18.5+blat,urcrnrlon=85+(rlon-180),urcrnrlat=18.5+tlat, 
                resolution='c', projection='merc', lat_0 = 18.5, lon_0 = 85)
            m.drawcoastlines()
            parallels = np.arange(-90,90,15.)
            m.drawparallels(parallels,labels=[1,0,0,0])
            meridians = np.arange(0.,360.,15.)
            m.drawmeridians(meridians,labels=[0,0,0,1])
            #quiver on top
            if ttest:
                u_da = xr.open_dataarray(ufilenames+str(lag)+'_ttest.nc')
                v_da = xr.open_dataarray(vfilenames+str(lag)+'_ttest.nc')
            else:
                u_da = xr.open_dataarray(ufilenames+str(lag)+'.nc')
                v_da = xr.open_dataarray(vfilenames+str(lag)+'.nc')
            u_comp = u_da.sel({'lv_ISBL0' : lvl}, method='nearest').sel(
                {'g4_lat_1' : slice(tlat, blat), 'g4_lon_2' : slice(llon, rlon)})
            v_comp = v_da.sel({'lv_ISBL0' : lvl}, method='nearest').sel(
                {'g4_lat_1' : slice(tlat, blat), 'g4_lon_2' : slice(llon, rlon)})
            x = u_comp.coords['g4_lon_2'].values - 180 + 85
            y = u_comp.coords['g4_lat_1'].values + 18.5
            x, dummy = m(x, x)
            dummy, y = m(y, y)
            m.quiver(x, y, u_comp.values, v_comp.values)
            #bounding box and center
            x, y = m([75, 95, 95, 75, 75], [10, 10, 27, 27, 10])
            m.plot(x, y, 'r-', lw=2.5)
            x, y = m([85], [18.5])
            m.plot(x, y, 'ro', ms=5)
            if i == 0 and j == 0:
                if ttest:
                    ax.set_title('{}: {} mb, ttest\nlag {}'.format(ufilenames, lvl, lag))
                else:
                    ax.set_title('{}: {} mb\nlag {}'.format(ufilenames, lvl, lag))
            else:
                ax.set_title('lag {}'.format(lag))
            lag_i += 1
    plt.tight_layout(pad=4)

def plot2dxz_lag(filenames, lags=[-2,-1,0,1,2], ttest=True):
    if len(lags) == 5:
        nR = 3; nC = 2
        f, axes = plt.subplots(nR, nC, figsize=(16, 10) )
    elif len(lags) == 3:
        nR = 3; nC = 1
        f, axes = plt.subplots(nR, nC, figsize=(8, 10) )
    llon = 130
    rlon = 230
    llvl = 1000
    ulvl = 1
    lag_i = 0
    for j in range(nC):
        for i in range(nR):
            if nC > 1 and nR > 1:
                ax = axes[i, j]
            elif nC == 1:
                ax = axes[i]
            else:
                ax = axes[j]
            if lag_i == len(lags):
                f.delaxes(ax)
                plt.draw()
                break
            lag = lags[lag_i]
            if ttest:
                da = xr.open_dataarray(filenames+str(lag)+'_ttest.nc')
            else:
                da = xr.open_dataarray(filenames+str(lag)+'.nc')
            comp = da.sel({'g4_lat_1' : 0}, method='nearest').sel({'g4_lon_2' : slice(llon, rlon), 'lv_ISBL0' : slice(ulvl, llvl)}).values
            cmap = plt.get_cmap('PiYG')
            avg = np.mean(np.abs(comp))
            sd = np.std(np.abs(comp))
            im = ax.imshow(comp, cmap=cmap, extent=(llon-180, rlon-180, llvl, ulvl), aspect=0.05,
                vmin=-avg - 3*sd, vmax=avg + 3*sd, interpolation='hermite')
            cb = plt.colorbar(im, ax=ax)
            ax.set_xlabel('Longitude ($^\\circ$ from genesis pt)')
            if i == 0 and j == 0:
                if ttest:
                    ax.set_title('{}: ttest\nlag {}'.format(filenames, lag))
                else:
                    ax.set_title('{}\nlag {}'.format(filenames, lag))
            else:
                ax.set_title('lag {}'.format(lag))
            lag_i += 1
    plt.tight_layout(pad=4)

lags202 = [-2, 0, 2]
lags21012 = [-2, -1, 0, 1, 2]

#comp_da = xr.open_dataarray('composite_n150_Vanom_JJAS.nc')
#plot1dx(comp_da,  'x, 850 mb, n = 150, anomalous V, JJAS')
#plot2dxy(comp_da, 'xy, 850 mb, n = 150, anomalous V, JJAS')
#plot2dxz(comp_da, 'xz, n = 150, anomalous V, JJAS')

#comp_da = xr.open_dataarray('composite_n10_V_JJAS.nc')
#plot2dxy(comp_da, 'xy, 850 mb, n = 10, V, JJAS')

#plot2dxy('composite_n150_Vanom_JJAS_lag0_ttest.nc', lvl=850)
#plot2dxy('composite_n150_Vanom_JJAS_lag0_ttest.nc', lvl=300)

#plot2dxy_lag('composite_n150_Vanom_JJAS_lag', lvl=850)
#plot2dxy_lag('composite_n150_Vanom_JJAS_lag', lvl=400)
#plot2dxy_lag('composite_n150_Vanom_JJAS_lag', lvl=300)
#plot2dxz_lag('composite_n150_Vanom_JJAS_lag')

#plot2dxy_lag('composite_n150_Vanom_JJAS_lag', lags=lags202, ttest=True, lvl=850)
#plot2dxy_lag('composite_n150_Vanom_JJAS_lag', lags=lags202, ttest=True, lvl=300)
#plot2dxy_lag('composite_n150_Uanom_JJAS_lag', lags=lags202, ttest=True, lvl=850)
#plot2dxy_lag('composite_n150_Uanom_JJAS_lag', lags=lags202, ttest=True, lvl=300)

#plot2dxy_lag_vector('composite_n150_Uanom_JJAS_lag', 'composite_n150_Vanom_JJAS_lag', lags=lags202, ttest=True, lvl=850)

#plot2dxy_vector('composite_n150_Uanom_JJAS_lag-2', 
#                'composite_n150_Vanom_JJAS_lag-2', ttest=False, lvl=850)
#plot2dxy_vector('composite_n150_Uanom_JJAS_lag0', 
#                'composite_n150_Vanom_JJAS_lag0', ttest=False, lvl=850)
#plot2dxy_vector('composite_n150_Uanom_JJAS_lag2', 
#                'composite_n150_Vanom_JJAS_lag2', ttest=False, lvl=850)
#plot2dxy_vector('composite_n150_Uanom_JJAS_lag-2', 
#                'composite_n150_Vanom_JJAS_lag-2',ttest=True, lvl=850)
plot2dxy_vector('composite_n150_Uanom_JJAS_lag0', 
                'composite_n150_Vanom_JJAS_lag0', ttest=True, lvl=850)
#plot2dxy_vector('composite_n150_Uanom_JJAS_lag2', 
#               'composite_n150_Vanom_JJAS_lag2', ttest=True, lvl=850)

#plot2dxy_lag_2im('composite_n150_Uanom_JJAS_lag', 'composite_n150_Uanom_JJAS_lag', lags=lags202, ttest=True, lvl=850)

#plot2dxz_lag('composite_n150_Vanom_JJAS_lag', ttest=True)
#plot2dxz_lag('composite_n150_Vanom_JJAS_lag', ttest=False)
plt.show()
