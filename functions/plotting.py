#! /usr/bin/env python3

"""
Author: Lori Garzio on 5/12/2021
Last modified: 5/12/2021
"""
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean as cmo
plt.rcParams.update({'font.size': 12})


def add_map_features(axis, extent, edgecolor=None, landcolor=None):
    edgecolor = edgecolor or 'black'
    landcolor = landcolor or 'tan'

    land = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor=edgecolor, facecolor=landcolor)

    state_lines = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')

    # Axes properties and features
    axis.set_extent(extent)
    axis.add_feature(land)
    axis.add_feature(cfeature.RIVERS)
    axis.add_feature(cfeature.LAKES)
    axis.add_feature(cfeature.BORDERS)
    axis.add_feature(state_lines, zorder=11, edgecolor=edgecolor)

    # Gridlines and grid labels
    gl = axis.gridlines(
        draw_labels=True,
        linewidth=.5,
        color='black',
        alpha=0.25,
        linestyle='--'
    )

    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = {'size': 11, 'color': 'black'}
    gl.ylabel_style = {'size': 11, 'color': 'black'}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER


def glider_track(fig, ax, ds, region, bathy=None, landcolor=None, title=None, current_glider_loc=None,
                 timevariable=None):
    bathy = bathy or None
    landcolor = landcolor or 'tan'
    title = title or None
    current_glider_loc = current_glider_loc or None
    timevariable = timevariable or 'time'

    extent = region['extent']

    if bathy:
        #levels = np.arange(-9000, 9100, 100)
        levels = np.arange(-5000, 5100, 50)
        bath_lat = bathy.variables['lat'][:]
        bath_lon = bathy.variables['lon'][:]
        try:
            bath_elev = bathy.variables['elevation'][:]
        except KeyError:
            bath_elev = bathy.variables['altitude'][:]
        plt.contourf(bath_lon, bath_lat, bath_elev,  levels, cmap=cmo.cm.topo, transform=ccrs.PlateCarree())

        levels = np.arange(-100, 0, 50)
        CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                         transform=ccrs.PlateCarree())
        ax.clabel(CS, [-100], inline=True, fontsize=7, fmt='%d')

    margs = dict()
    margs['landcolor'] = landcolor
    add_map_features(ax, extent, **margs)

    # plot full glider track
    try:
        lon = ds.longitude.values
        lat = ds.latitude.values
    except AttributeError:
        lon = ds.Longitude.values
        lat = ds.Latitude.values
    ax.scatter(lon, lat, color='k', marker='.', s=60, transform=ccrs.PlateCarree(), zorder=10)
    sct = ax.scatter(lon, lat, c=ds[timevariable].values, marker='.', s=15, cmap='rainbow', transform=ccrs.PlateCarree(),
                     zorder=10)
    if current_glider_loc:
        ax.plot(lon[-1], lat[-1], color='white', marker='^', markeredgecolor='black',
                markersize=8.5, transform=ccrs.PlateCarree())

    # Plot title
    if title:
        plt.title(title)

    # Set colorbar height equal to plot height
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
    fig.add_axes(cax)

    # generate colorbar
    cbar = plt.colorbar(sct, cax=cax)
    cbar.ax.set_yticklabels(pd.to_datetime(cbar.ax.get_yticks()).strftime(date_format='%Y-%m-%d'))

    return fig, ax


def profile_yos(fig, ax, da, yo_idx, colors, xlabel=None, ylabel=None, title=None, grid=None):
    """
    Plot profiles with the paired down-up profiles (yos) the same color
    :param fig:
    :param ax:
    :param da:
    :param segment_idx:
    :param colors:
    :param xlabel:
    :param ylabel:
    :param title:
    :param grid:
    :return:
    """
    xlabel = xlabel or None
    ylabel = ylabel or 'Depth (m)'
    title = title or None
    grid = grid or False

    # loop through each profile to plot each down-up segment with the same color
    for i, ii in enumerate(yo_idx):  # each yo
        c = colors[i]
        for j, jj in enumerate(ii):  # each profile
            try:
                x = da.values[jj:ii[j + 1]]
                y = da.depth.values[jj:ii[j + 1]]
                df = pd.DataFrame({'x': x, 'y': y})
                # interpolate depth
                df['y'] = df['y'].interpolate(method='linear', limit_direction='both')
                df.dropna(subset=['x'], inplace=True)

                ax.plot(df.x, df.y, c=c)  # plot lines
                ax.scatter(df.x, df.y, color=c, s=30, edgecolor='None')
            except IndexError:
                continue

    ax.invert_yaxis()
    ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)

    if title:
        ax.set_title(title, fontsize=14)

    if grid:
        ax.grid(ls='--', lw=.5)


def xsection(fig, ax, x, y, z, xlabel=None, ylabel=None, clabel=None, cmap=None, title=None, date_fmt=None,
             grid=None, extend=None):
    xlabel = xlabel or 'Time'
    ylabel = ylabel or 'Depth (m)'
    clabel = clabel or None
    cmap = cmap or 'jet'
    title = title or None
    date_fmt = date_fmt or None
    grid = grid or False
    extend = extend or 'both'

    xc = ax.scatter(x, y, c=z, cmap=cmap, s=10, edgecolor='None')

    ax.invert_yaxis()
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    if title:
        ax.set_title(title, fontsize=16)

    # format colorbar
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
    fig.add_axes(cax)
    if clabel:
        cb = plt.colorbar(xc, cax=cax, label=clabel, extend=extend)
    else:
        cb = plt.colorbar(xc, cax=cax, extend=extend)

    # format x-axis
    if date_fmt:
        xfmt = mdates.DateFormatter(date_fmt)
        ax.xaxis.set_major_formatter(xfmt)

    if grid:
        ax.grid(ls='--', lw=.5)
