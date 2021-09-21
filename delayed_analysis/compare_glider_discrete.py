#!/usr/bin/env python

"""
Author: Lori Garzio on 9/16/2021
Last modified: 9/16/2021
Compare glider data to discrete water samples collected during glider deployment/recovery
"""

from pathlib import Path
from geographiclib.geodesic import Geodesic
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import os
import seawater as sw
import PyCO2SYS as pyco2
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})


def main(fname):
    basedir = Path().absolute().parent
    discrete_samples = os.path.join(basedir, 'water_sampling', 'NOAA_OA_pH_discrete_samples.csv')
    df = pd.read_csv(discrete_samples)

    # drop rows without pH values
    df = df.dropna(axis=0, how='all', subset=['AveragepH_25degC'])

    # create date times
    df['date_time'] = df['date'] + 'T' + df['gmt_time']
    df['date_time'] = df['date_time'].map(lambda t: dt.datetime.strptime(t, '%m/%d/%yT%H:%M'))

    # convert lat/lon to decimal degrees  (negative longitude because these are degrees W)
    df['lat_decimal'] = df['lat'].map(lambda x: float(str(x).split(' ')[-1])/60 + float(str(x).split(' ')[0]))
    df['lon_decimal'] = df['lon'].map(lambda x: -(float(str(x).split(' ')[-1]) / 60 + float(str(x).split(' ')[0])))

    # calculate pressure from depth
    df['pressure_dbar'] = sw.pres(df['depth_m'], df['lat_decimal'])

    # calculate pH corrected for temperature pressure salinity
    # pyCO2SYS needs two parameters to calculate pH corrected, so if AverageTA isn't available, fill with 2200
    df['AverageTA'] = df['AverageTA'].fillna(2200)
    par1 = df['AveragepH_25degC']
    par1_type = 3  # parameter type (pH)
    par2 = df['AverageTA']
    par2_type = 1

    kwargs = dict(salinity=df['sal'],
                  temperature=25,
                  temperature_out=df['temp_c'],
                  pressure=0,
                  pressure_out=df['pressure_dbar'],
                  opt_pH_scale=1,
                  opt_k_carbonic=4,
                  opt_k_bisulfate=1,
                  opt_total_borate=1,
                  opt_k_fluoride=2)


    results = pyco2.sys(par1, par2, par1_type, par2_type, **kwargs)
    df['ph_corrected'] = results['pH_out']

    ds = xr.open_dataset(fname)
    ds = ds.sortby(ds.time)
    ds = ds.swap_dims({'time': 'profile_lat'})
    deploy = ds.deployment_name
    df2 = df.loc[df['deployment'] == deploy]

    plt_vars = ['ph', 'salinity', 'temperature']
    for dr in np.unique(df2['type']):
        df_dr = df2.loc[df2['type'] == dr]
        for sample_time in np.unique(df_dr['date_time']):
            df_drt = df_dr.loc[df_dr['date_time'] == sample_time]
            sample_lat = np.unique(df_drt['lat_decimal'])
            sample_lon = np.unique(df_drt['lon_decimal'])
            if len(sample_lat) == 1 and len(sample_lon) == 1:
                # discrete water sample metadata
                st_str = pd.to_datetime(sample_time).strftime('%Y-%m-%dT%H:%M')
                st_save_str = pd.to_datetime(sample_time).strftime('%Y%m%dT%H%M')
                slat = np.round(sample_lat[0], 4)
                slon = np.round(sample_lon[0], 4)
                sample_meta = f'Sample: {st_str}, location {[slat, slon]}'

                # find the closest glider profile to the sample collection
                a = abs(sample_lat[0] - ds.profile_lat.values) + abs(sample_lon[0] - ds.profile_lon.values)
                idx = np.unravel_index(a.argmin(), a.shape)
                glider_profile_lat = ds.profile_lat.values[idx]
                dss = ds.sel(profile_lat=glider_profile_lat)
                dss_t0 = pd.to_datetime(np.nanmin(dss.time.values))
                dss_tf = pd.to_datetime(np.nanmax(dss.time.values))

                # glider metadata
                dss_t0str = pd.to_datetime(dss_t0).strftime('%Y-%m-%dT%H:%M')
                gllat = np.round(np.unique(dss.profile_lat.values[0]), 4)[0]
                gllon = np.round(np.unique(dss.profile_lon.values[0]), 4)[0]
                glider_meta = f'Glider profile: {dss_t0str}, location {[gllat, gllon]}'

                # calculate distance between glider and samples
                geod = Geodesic.WGS84
                g = geod.Inverse(gllat, gllon, slat, slon)
                diff_loc_meters = np.round(g['s12'], 2)
                diff_seconds = np.round(abs(dss_t0 - pd.to_datetime(sample_time)).seconds / 60, 2)
                diff_meta = f'Difference: {diff_seconds} seconds, {diff_loc_meters} meters'

                for pv in plt_vars:
                    # get the discrete water sample data
                    sample_depth = df_drt['depth_m']
                    if pv == 'ph':
                        sample = df_drt['ph_corrected']
                    elif pv == 'salinity':
                        sample = df_drt['sal']
                    elif pv == 'temperature':
                        sample = df_drt['temp_c']

                    fig, ax = plt.subplots(figsize=(8, 10))

                    # plot glider data
                    if pv == 'ph':
                        glider_vars = ['ph_total', 'ph_total_shifted']
                    else:
                        glider_vars = [pv]
                    colors = ['tab:blue', 'tab:green']
                    for i, variable in enumerate(glider_vars):
                        glider_data = dss[variable]
                        ax.scatter(glider_data, glider_data.depth, c=colors[i], label=variable)

                    # plot water sample data
                    ax.scatter(sample.astype(float), sample_depth.astype(float), c='tab:red', ec='k', s=70,
                               label='water samples')

                    ax.legend()
                    if pv == 'ph':
                        xlims = [7.6, 8.15]
                        xlab = 'pH'
                    elif pv == 'salinity':
                        xlims = [28.5, 32.5]
                        xlab = 'Salinity'
                    elif pv == 'temperature':
                        xlims = [14.75, 25.25]
                        xlab = 'Temperature (degrees C)'

                    plt.ylim(0, 16)
                    plt.xlim(xlims)
                    ax.set_xlabel(xlab)
                    ax.invert_yaxis()
                    ax.set_ylabel('Depth (m)')
                    ax.set_title(f'Glider {dr.capitalize()}\n{sample_meta}\n{glider_meta}\n{diff_meta}')

                    sfilename = f'{deploy}_discrete_comparison_{st_save_str}_{pv}.png'
                    sfile = os.path.join(os.path.dirname(fname), 'compare_glider_discrete', sfilename)
                    plt.savefig(sfile, dpi=300)
                    plt.close()

            else:
                raise ValueError(f'too many lat/lons for samples collected at {st_str}')


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210716T1804/delayed/ru30-20210716T1804-profile-sci-delayed_shifted.nc'
    main(ncfile)
