#!/usr/bin/env python

"""
Author: Lori Garzio on 9/16/2021
Last modified: 7/19/2022
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
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console


def main(fname):
    save_dir = os.path.join(os.path.dirname(fname), 'compare_glider_discrete')
    os.makedirs(save_dir, exist_ok=True)

    summary_headers = ['deployment_recovery', 'glider_date', 'discrete_date', 'time_difference_minutes',
                       'glider_pressure_dbar', 'discrete_pressure_dbar', 'glider_lon', 'glider_lat',
                       'discrete_lon', 'discrete_lat', 'distance_m', 'glider_ph', 'discrete_ph', 'diff_ph',
                       'glider_ta', 'discrete_ta', 'diff_ta', 'glider_temp', 'discrete_temp', 'glider_sal',
                       'discrete_sal']
    summary_rows = []

    basedir = Path().absolute().parent
    discrete_samples = os.path.join(basedir, 'water_sampling', 'NOAA_OA_pH_discrete_samples.csv')
    df = pd.read_csv(discrete_samples)

    # drop rows without pH values or lat/lon
    df = df.dropna(axis=0, how='all', subset=['AveragepH_25degC'])
    df = df.dropna(axis=0, how='all', subset=['lat', 'lon'])

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
    filename = fname.split('/')[-1]
    deploy = f'{filename.split("-")[0]}-{filename.split("-")[1]}'
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
                data = dict(glat=ds.profile_lat.values, glon=ds.profile_lon.values,
                            slat=np.repeat(slat, len(ds.profile_lat.values)),
                            slon=np.repeat(slon, len(ds.profile_lat.values)))
                df = pd.DataFrame(data)
                df.drop_duplicates(inplace=True)

                # calculate distances from each glider profile to sample location
                geod = Geodesic.WGS84
                distances = np.array([])
                for ii, row in df.iterrows():
                    g = geod.Inverse(row.glat, row.glon, row.slat, row.slon)
                    distances = np.append(distances, g['s12'])

                # find the smallest distance
                idx = distances.argmin()
                plat = df.iloc[idx].glat
                plat_idx = np.where(ds.profile_lat.values == plat)[0]

                # subset the data
                dss = ds.isel(time=plat_idx)

                # determine if the glider profile is ascending or descending to grab the corresponding profile
                d0 = np.nanmedian(dss.depth.values[0:10])
                d1 = np.nanmedian(dss.depth.values[-10:-1])

                if d0 < d1:  # descending profile, grab the next profile to complete the pair
                    idx2 = idx + 1
                else:  # ascending profile, grab the previous profile to complete the pair
                    idx2 = idx - 1
                plat2 = df.iloc[idx2].glat
                plat_idx2 = np.where(np.logical_or(ds.profile_lat.values == plat, ds.profile_lat.values == plat2))[0]

                # subset the data again to get the down-up profile pair
                dss = ds.isel(time=plat_idx2)

                dss_t0 = pd.to_datetime(np.nanmin(dss.time.values))

                # glider metadata
                dss_t0str = pd.to_datetime(dss_t0).strftime('%Y-%m-%dT%H:%M')
                gllat = np.round(np.unique(dss.profile_lat.values[0]), 4)[0]
                gllon = np.round(np.unique(dss.profile_lon.values[0]), 4)[0]
                glider_meta = f'Glider profile: {dss_t0str}, location {[gllat, gllon]}'

                # differences between glider and samples - distance and time
                diff_loc_meters = np.round(distances[idx], 2)
                diff_mins = np.round(abs(dss_t0 - pd.to_datetime(sample_time)).seconds / 60, 2)
                diff_meta = f'Difference: {diff_mins} minutes, {diff_loc_meters} meters'

                for pv in plt_vars:
                    # get the discrete water sample data
                    sample_pressure = df_drt['pressure_dbar']
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
                        ax.scatter(glider_data, dss.pressure_interpolated, c=colors[i], label=f'glider {variable}')

                    # plot water sample data
                    ax.scatter(sample.astype(float), sample_pressure.astype(float), c='tab:red', ec='k', s=70,
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
                    #plt.xlim(xlims)
                    ax.set_xlabel(xlab)
                    ax.invert_yaxis()
                    ax.set_ylabel('Pressure (dbar)')
                    ax.set_title(f'Glider {dr.capitalize()}\n{sample_meta}\n{glider_meta}\n{diff_meta}')

                    sfilename = f'{deploy}_discrete_comparison_{dr}_{st_save_str}_{pv}.png'
                    sfile = os.path.join(save_dir, sfilename)
                    plt.savefig(sfile, dpi=300)
                    plt.close()

                # write summary file
                discrete_pressure_unique = np.unique(df_drt.pressure_dbar)
                for dpu in discrete_pressure_unique:
                    sdf = df_drt.loc[df_drt['pressure_dbar'] == dpu]
                    discrete_pressure = np.nanmedian(sdf.pressure_dbar)
                    if discrete_pressure < 1:
                        gl_pressure_idx = np.where(dss.pressure_interpolated < 3)[0]
                    else:
                        press1 = discrete_pressure - 1.5
                        press2 = discrete_pressure + 1.5
                        gl_pressure_idx = np.where(np.logical_and(dss.pressure_interpolated > press1,
                                                                  dss.pressure_interpolated < press2))[0]
                    gl_pressure = np.round(np.nanmedian(dss.pressure_interpolated[gl_pressure_idx]), 3)

                    discrete_ph = np.round(np.nanmedian(sdf.ph_corrected), 3)
                    glider_ph = np.round(np.nanmedian(dss.ph_total_shifted[gl_pressure_idx]), 3)
                    ph_diff = np.round(glider_ph - discrete_ph, 3)
                    discrete_ta = np.round(np.nanmedian(sdf.AverageTA), 3)
                    glider_ta = np.round(np.nanmedian(dss.total_alkalinity[gl_pressure_idx]), 3)
                    ta_diff = np.round(glider_ta - discrete_ta, 3)
                    discrete_temp = np.round(np.nanmedian(sdf.temp_c), 3)
                    glider_temp = np.round(np.nanmedian(dss.temperature_interpolated[gl_pressure_idx]), 3)
                    discrete_sal = np.round(np.nanmedian(sdf.sal), 3)
                    glider_sal = np.round(np.nanmedian(dss.salinity_interpolated[gl_pressure_idx]), 3)

                    summary_data = [dr, dss_t0str, st_str, diff_mins, gl_pressure, discrete_pressure, gllon, gllat,
                                    slon, slat, diff_loc_meters, glider_ph, discrete_ph, ph_diff, glider_ta, discrete_ta,
                                    ta_diff, glider_temp, discrete_temp, glider_sal, discrete_sal]

                    summary_rows.append(summary_data)

            else:
                raise ValueError(f'too many lat/lons for samples collected at {st_str}')

    summary_df = pd.DataFrame(summary_rows, columns=summary_headers)
    summary_df.to_csv(os.path.join(save_dir, f'{deploy}_compare_glider_discrete.csv'), index=False)


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed_qc.nc'
    main(ncfile)
