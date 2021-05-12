#!/usr/bin/env python

"""
Author: Lori Garzio on 5/12/2021
Last modified: 5/12/2021
Plot realtime glider track
"""

import argparse
import sys
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import functions.common as cf
import functions.plotting as pf
plt.rcParams.update({'font.size': 12})


#def main(deploy, sfilename):
def main(args):
    ru_server = 'http://slocum-data.marine.rutgers.edu//erddap'
    # bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    bathymetry = '/home/lgarzio/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'  # on server
    deploy = args.deployment
    glider_id = '{}-profile-sci-rt'.format(deploy)
    glider_vars = ['latitude', 'longitude']
    ds = cf.get_erddap_dataset(ru_server, glider_id, glider_vars)
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)

    glider_region = cf.glider_region(ds)  # define the glider region

    fig, ax = plt.subplots(figsize=(11, 8), subplot_kw=dict(projection=ccrs.Mercator()))
    plt.subplots_adjust(right=0.82)

    t0 = pd.to_datetime(np.nanmin(ds.time)).strftime('%Y-%m-%dT%H:%M')
    tf = pd.to_datetime(np.nanmax(ds.time)).strftime('%Y-%m-%dT%H:%M')
    now = dt.datetime.now().strftime('%d-%b-%Y %H:%M')

    title = f'{deploy.split("-")[0]} track: {t0} to {tf}\nUpdated: {now}'

    extent = glider_region['extent']
    bathy = xr.open_dataset(bathymetry)

    kwargs = dict()
    kwargs['bathy'] = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                                lat=slice(extent[2] - .1, extent[3] + .1))

    kwargs['title'] = title
    kwargs['landcolor'] = 'none'

    fig, ax = pf.glider_track(fig, ax, ds, glider_region, **kwargs)

    plt.savefig(args.sfilename, dpi=200)
    plt.close()


if __name__ == '__main__':
    # deployment = 'ru30-20210503T1929'
    # savefile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210503T1929/rt_plotting/phoxy_live_track.png'
    # main(deployment, savefile)
    arg_parser = argparse.ArgumentParser(description='Plot real time glider track',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-deploy',
                            dest='deployment',
                            type=str,
                            help='Glider deployment to plot (glider-yyyymmddTHHMM), e.g. ru30-20210226T1647')

    arg_parser.add_argument('-s', '--save_file',
                            dest='sfilename',
                            type=str,
                            help='Full file path to save directory and save filename')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
