#!/usr/bin/env python

"""
Author: Lori Garzio on 4/6/2021
Last modified: 8/10/2021
"""

import numpy as np
import os
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import functions.common as cf
import functions.plotting as pf
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console
plt.rcParams.update({'font.size': 13})


def main(deploy, fname):
    # initialize keyword arguments for plotting function
    kwargs = dict()
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    # bathymetry = '/home/lgarzio/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'  # on server
    if fname.endswith('_to_dac.nc'):
        savedir = os.path.join(os.path.dirname(fname), 'dac_plots')
        timevar = 'UnixTime'
        kwargs['timevariable'] = 'UnixTime'
    else:
        savedir = os.path.join(os.path.dirname(fname), 'plots')
        timevar = 'time'
    os.makedirs(savedir, exist_ok=True)

    ds = xr.open_dataset(fname)
    glider_region = cf.glider_region(ds)  # define the glider region

    fig, ax = plt.subplots(figsize=(11, 8), subplot_kw=dict(projection=ccrs.Mercator()))
    plt.subplots_adjust(right=0.82)

    t0 = pd.to_datetime(np.nanmin(ds[timevar])).strftime('%Y-%m-%dT%H:%M')
    tf = pd.to_datetime(np.nanmax(ds[timevar])).strftime('%Y-%m-%dT%H:%M')

    title = f'{deploy}\n{t0} to {tf}'

    extent = glider_region['extent']
    bathy = xr.open_dataset(bathymetry)

    kwargs['bathy'] = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                                lat=slice(extent[2] - .1, extent[3] + .1))

    kwargs['title'] = title
    kwargs['landcolor'] = 'none'
    pf.glider_track(fig, ax, ds, glider_region, **kwargs)

    sname = os.path.join(savedir, '{}_glider-track.png'.format(deploy))
    plt.savefig(sname, dpi=200)
    plt.close()


if __name__ == '__main__':
    deployment = 'ru30-20210226T1647'
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647_to_dac.nc'
    main(deployment, ncfile)
