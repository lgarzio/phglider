#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2021
Last modified: 9/18/2024
Plot locations of vessel-based water samples for TA-salinity regressions
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import functions.plotting as pf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
plt.rcParams.update({'font.size': 15})


def main(file, savedir):
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-76, -70, 37.5, 41.5]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    ds = xr.open_dataset(file)

    # combine datasets

    years = np.unique(pd.to_datetime(ds.time).year)
    min_year = np.min(years)
    max_year = np.max(years)
    year_list = [y for y in range(min_year, max_year + 1)]
    cruise_list = np.unique(ds.cruise.values)

    # plot map of summer bottom pH data for entire dataset
    fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
    plt.subplots_adjust(top=.92, bottom=0.08, right=.96, left=0.06)

    # define bathymetry levels and data
    bath_lat = bathy.variables['lat'][:]
    bath_lon = bathy.variables['lon'][:]
    bath_elev = bathy.variables['elevation'][:]

    levels = [-3000, -1000, -100]
    CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                     transform=ccrs.PlateCarree())
    ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

    pf.add_map_features(ax, extent)

    sct = ax.scatter(ds.lon.values, ds.lat.values, c='k', marker='.',
                     s=50, transform=ccrs.PlateCarree(), zorder=10)

    plt.title(f'Vessel-based TA sampling: {min_year} to {max_year}')

    sfile = os.path.join(savedir, f'vessel_based_TA_salinity-NYB.png')

    plt.savefig(sfile, dpi=200)
    plt.close()


if __name__ == '__main__':
    file = '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/vessel_based_TA_salinity_NYB_20241002.nc'
    save_directory = '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data'
    main(file, save_directory)