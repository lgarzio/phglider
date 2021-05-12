#! /usr/bin/env python3

"""
Author: Lori Garzio on 5/12/2021
Last modified: 5/12/2021
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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


def glider_track(fig, ax, ds, region, bathy=None, landcolor=None, title=None, current_glider_loc=None):
    bathy = bathy or None
    landcolor = landcolor or 'tan'
    title = title or None
    current_glider_loc = current_glider_loc or None

    extent = region['extent']

    if bathy:
        #levels = np.arange(-9000, 9100, 100)
        levels = np.arange(-5000, 5100, 50)
        bath_lat = bathy.variables['lat'][:]
        bath_lon = bathy.variables['lon'][:]
        bath_elev = bathy.variables['elevation'][:]
        plt.contourf(bath_lon, bath_lat, bath_elev,  levels, cmap=cmo.cm.topo, transform=ccrs.PlateCarree())

        levels = np.arange(-100, 0, 50)
        CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                         transform=ccrs.PlateCarree())
        ax.clabel(CS, [-100], inline=True, fontsize=7, fmt='%d')

    margs = dict()
    margs['landcolor'] = landcolor
    add_map_features(ax, extent, **margs)

    # plot full glider track
    ax.scatter(ds.longitude.values, ds.latitude.values, color='white', marker='.', s=40, transform=ccrs.PlateCarree(),
               zorder=10)
    sct = ax.scatter(ds.longitude.values, ds.latitude.values, c=ds.time.values, marker='.', s=15, cmap='rainbow',
                     transform=ccrs.PlateCarree(), zorder=10)
    if current_glider_loc:
        ax.plot(ds.longitude.values[-1], ds.latitude.values[-1], color='white', marker='^', markeredgecolor='black',
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
