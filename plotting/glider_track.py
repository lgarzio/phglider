#!/usr/bin/env python

"""
Author: Lori Garzio on 4/6/2021
Last modified: 4/6/2021
"""

import numpy as np
import os
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy.feature as cfeature
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console
plt.rcParams.update({'font.size': 13})


def add_map_features(ax, axes_limits, land_color=None):
    """
    Adds latitude and longitude gridlines and labels, coastlines, and statelines to a cartopy map object
    :param ax: plotting axis object
    :param axes_limits: list of axis limits [min lon, max lon, min lat, max lat]
    :param land_color: optional color for land
    """
    #gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='dotted', x_inline=False)
    gl = ax.gridlines(draw_labels=True, linewidth=0)
    gl.top_labels = False
    gl.right_labels = False
    #gl.rotate_labels = False
    gl.xpadding = 12
    gl.ypadding = 12
    ax.set_extent(axes_limits)

    if land_color:
        lc = land_color
    else:
        lc = cfeature.COLORS['land_alt1']

    land = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='black', facecolor=lc)
    ax.add_feature(land, zorder=1)

    state_lines = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')

    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(state_lines, edgecolor='black')

    feature = cfeature.NaturalEarthFeature(
        name='coastline', category='physical',
        scale='10m',
        edgecolor='black', facecolor='none')
    ax.add_feature(feature, zorder=8)


def main(deploy, coord_lims, fname):
    savedir = os.path.join(os.path.dirname(fname), 'plots')
    os.makedirs(savedir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 9), subplot_kw=dict(projection=ccrs.PlateCarree()))
    plt.subplots_adjust(right=0.82)
    plt.title('{}'.format(deploy))
    axlims = [coord_lims['longitude']['min'], coord_lims['longitude']['max'],
              coord_lims['latitude']['min'], coord_lims['latitude']['max']]
    add_map_features(ax, axlims, 'darkgray')

    # add bathymetry
    bathydir = '/Users/garzio/Documents/repo/lgarzio/phglider/bathymetry'
    gf = os.path.join(bathydir, 'GMRTv3_9_20210406topo-NJ.grd')
    grid_file = xr.open_dataset(gf)
    gf_lon = grid_file['lon']
    gf_lat = grid_file['lat']
    gf_lon_ind = np.logical_and(gf_lon > coord_lims['longitude']['min']-.05, gf_lon < coord_lims['longitude']['max']+.05)
    gf_lat_ind = np.logical_and(gf_lat > coord_lims['latitude']['min']-.05, gf_lat < coord_lims['latitude']['max']+.05)
    bathy = grid_file['altitude'][gf_lat_ind, gf_lon_ind].values

    land_mask = np.logical_and(bathy >= -1, bathy >= -1)
    bathy[land_mask] = np.nan  # turn values over land to nans
    depth = abs(bathy)
    levs = [25, 50, 100, 1000, 2000, 3000, 4000]
    cs = ax.contour(gf_lon[gf_lon_ind], gf_lat[gf_lat_ind], depth, levs, colors='gray', linestyles='solid',
                    linewidths=.5, alpha=.5)
    ax.clabel(cs, inline=True, fontsize=9, fmt='%d')
    # ln_depth = np.log(depth)  # for filled contours
    # levs = [1, 25, 50, 100, 500, 1000, 2000, 3000, 4000]
    # ln_levs = np.log(levs)  # for filled contours
    # cs = ax.contourf(gf_lon[gf_lon_ind], gf_lat[gf_lat_ind], ln_depth, ln_levs, cmap='Blues', alpha=.75)

    # add glider track
    ds = xr.open_dataset(fname)
    loc_idx = np.logical_and(ds.qartod_longitude_location_flag == 1, ds.qartod_latitude_location_flag == 1)
    gllon = ds.longitude[loc_idx]
    gllat = ds.latitude[loc_idx]
    gltm = ds.time[loc_idx]
    sct = plt.scatter(gllon, gllat, c=gltm, marker='.', s=5, cmap='rainbow', transform=ccrs.PlateCarree())

    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
    fig.add_axes(cax)
    cbar = plt.colorbar(sct, cax=cax)
    cbar.ax.set_yticklabels(pd.to_datetime(cbar.ax.get_yticks()).strftime(date_format='%Y-%m-%d'))

    sname = os.path.join(savedir, '{}_glider-track.png'.format(deploy))
    plt.savefig(sname, dpi=200)
    plt.close()


if __name__ == '__main__':
    deployment = 'ru30-20210226T1647'
    coordinate_lims = {'latitude': {'min': 37, 'max': 41.5}, 'longitude': {'min': -76, 'max': -71}}
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed_shifted_qc.nc'
    main(deployment, coordinate_lims, ncfile)
