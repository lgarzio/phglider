#!/usr/bin/env python

"""
Author: Lori Garzio on 4/2/2021
Last modified: 4/2/2021
"""

import os
import numpy as np
import pandas as pd
import datetime as dt
import json
import xarray as xr
import time
from erddapy import ERDDAP
import cmocean as cmo
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import functions.phcalc as phcalc
import functions.common as cf
plt.rcParams.update({'font.size': 16})
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console



def main(deploy, sensor_sn, shifts, fname):
    ds = xr.open_dataset(fname)
    ds = ds.drop_vars(names=['profile_id', 'rowSize'])
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)

    # read cal file
    calfile = cf.find_calfile(deploy, sensor_sn)
    with open(calfile) as json_file:
        cc = json.load(json_file)

    # drop timestamps with missing data
    df = ds.to_dataframe()

    # convert profile_time to seconds since 1970-01-01 so you don't lose the column in resampling
    df['profile_time'] = df['profile_time'].astype('int64')//1e9
    df = df.resample('1s').mean()
    col_subset = ['conductivity', 'salinity', 'sci_water_pressure', 'temperature', 'sbe41n_ph_ref_voltage',
                           'chlorophyll_a', 'oxygen_concentration']
    for key, val in shifts.items():
        if val > 0:  # if there is a shift calculated
            if key == 'ph':
                varname = 'sbe41n_ph_ref_voltage'
                ph_shiftname = '{}_shifted_{}s'.format(varname, val)
            elif key == 'oxygen':
                varname = 'oxygen_concentration'
                oxy_shiftname = '{}_shifted_{}s'.format(varname, val)
            da = pd.DataFrame(df[varname])
            tm_shift = df.index - dt.timedelta(seconds=val)
            da['time_shift'] = tm_shift
            da.reset_index(drop=True, inplace=True)  # drop the original time index
            da = da.rename(columns={varname: '{}_shifted_{}s'.format(varname, val), 'time_shift': 'time'})
            da.set_index('time', inplace=True)
            df = df.merge(da, how='left', left_index=True, right_index=True)
            col_subset.append('{}_shifted_{}s'.format(varname, val))

    df = df.dropna(axis=0, how='all', subset=col_subset)

    # calculate pressure_dbar and drop any rows where pressure is <1 dbar
    df['sci_water_pressure_dbar'] = df['sci_water_pressure'] * 10
    df = df[df.sci_water_pressure_dbar > 1]

    # calculate pH and add to dataframe
    df['f_p'] = np.polyval([cc['f6'], cc['f5'], cc['f4'], cc['f3'], cc['f2'], cc['f1'], 0], df.sci_water_pressure_dbar)
    phfree, phtot = phcalc.phcalc(df.sbe41n_ph_ref_voltage, df.sci_water_pressure_dbar, df.temperature, df.salinity,
                                 cc['k0'], cc['k2'], df.f_p)
    df['ph_total'] = phtot
    df[df['ph_total'] > 14] = np.nan

    try:
        phfree_sh, phtot_sh = phcalc.phcalc(df[ph_shiftname], df.sci_water_pressure_dbar, df.temperature, df.salinity,
                                            cc['k0'], cc['k2'], df.f_p)
        df['ph_total_shifted'] = phtot_sh
        df[df['ph_total_shifted'] > 14] = np.nan
    except NameError:
        print('no pH shift defined')

    ds_shifted = df.to_xarray()

    # assign attributes from the original ds
    for varname, da in ds_shifted.data_vars.items():
        try:
            ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds[varname].attrs)
        except KeyError:
            if 'sbe41n_ph_ref_voltage' in varname:
                ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds['sbe41n_ph_ref_voltage'].attrs)
                ds_shifted[varname].attrs['comment'] = 'Values time shifted by -{} seconds to correct for thermal lag.'.format(sensor_shifts['ph'])
            if 'oxygen_concentration' in varname:
                ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds['oxygen_concentration'].attrs)
                ds_shifted[varname].attrs['comment'] = 'Values time shifted by -{} seconds to correct for thermal lag.'.format(sensor_shifts['oxygen'])
            if 'sci_water_pressure' in varname:
                ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds['sci_water_pressure'].attrs)
                for item in ['actual_range', 'colorBarMaximum', 'colorBarMinimum', 'valid_max', 'valid_min']:
                    del ds_shifted[varname].attrs[item]
                ds_shifted[varname].attrs['observation_type'] = 'calculated'
                ds_shifted[varname].attrs['units'] = 'dbar'
            if 'f_p' in varname:
                ds_shifted[varname].attrs['observation_type'] = 'calculated'
                ds_shifted[varname].attrs['comment'] = 'Polynomial evaluation of sensor pressure response polynomial coefficients (f6-f1) and pressure_dbar'
                ds_shifted[varname].attrs['f6'] = cc['f6']
                ds_shifted[varname].attrs['f5'] = cc['f5']
                ds_shifted[varname].attrs['f4'] = cc['f4']
                ds_shifted[varname].attrs['f3'] = cc['f3']
                ds_shifted[varname].attrs['f2'] = cc['f2']
                ds_shifted[varname].attrs['f1'] = cc['f1']
                ds_shifted[varname].attrs['units'] = 1
            if varname == 'ph_total':
                ds_shifted[varname].attrs['observation_type'] = 'calculated'
                ds_shifted[varname].attrs['units'] = 1
                ds_shifted[varname].attrs['k0'] = cc['k0']
                ds_shifted[varname].attrs['k2'] = cc['k2']
                ds_shifted[varname].attrs['comment'] = 'pH calculated on the total scale using pH reference voltage, pressure_dbar ,' \
                                                       'temperature, salinity, sensor calibration coefficients, and pressure ' \
                                                       'polynomial coefficients.'
            if varname == 'ph_total_shifted':
                ds_shifted[varname].attrs['observation_type'] = 'calculated'
                ds_shifted[varname].attrs['units'] = 1
                ds_shifted[varname].attrs['k0'] = cc['k0']
                ds_shifted[varname].attrs['k2'] = cc['k2']
                ds_shifted[varname].attrs['comment'] = 'pH calculated on the total scale using time shifted pH reference voltage, pressure_dbar, ' \
                                                       'temperature, salinity, sensor calibration coefficients, and pressure ' \
                                                       'polynomial coefficients.'

    ds_shifted = ds_shifted.assign_coords(
        {'latitude': ds_shifted.latitude, 'longitude': ds_shifted.longitude, 'depth': ds_shifted.depth})

    ds_global_attrs = ds.attrs
    ds_shifted = ds_shifted.assign_attrs(ds_global_attrs)
    ds_shifted.to_netcdf('{}_shifted.nc'.format(fname.split('.')[0]))

    # plot
    # df.reset_index(inplace=True)
    # daterange = pd.date_range(df.time[0], df.time[len(df) - 1], freq='3H')
    # for dt_idx, dr in enumerate(daterange):
    #     if dt_idx > 0:
    #         df_dr = df[(pd.to_datetime(df.time) < dr) & (pd.to_datetime(df.time) >= daterange[dt_idx - 1])]
    #         for pv in ['ph', 'oxygen']:
    #             if pv == 'ph':
    #                 colname = 'ph_total'
    #                 colname_shift = 'ph_total_shifted'
    #             else:
    #                 colname = 'oxygen_concentration'
    #                 colname_shift = oxy_shiftname
    #
    #             fig, (ax1, ax2) = plt.subplots(2, figsize=(12, 12), sharex=True, sharey=True)
    #             # xc = ax1.scatter(tm, depth, c=oxy, cmap=cmo.cm.oxy, s=10, edgecolor='None')
    #             # if pv == 'ph':
    #             #     xc = ax1.scatter(df_dr.time, df_dr.depth, c=df_dr[colname], cmap=cm.jet, s=10, edgecolor='None',
    #             #                      vmin=7.92, vmax=8.19)
    #             # else:
    #             #     xc = ax1.scatter(df_dr.time, df_dr.depth, c=df_dr[colname], cmap=cm.jet, s=10, edgecolor='None')
    #             xc = ax1.scatter(df_dr.time, df_dr.depth, c=df_dr[colname], cmap=cm.jet, s=10, edgecolor='None')
    #
    #             # format colorbar
    #             divider = make_axes_locatable(ax1)
    #             cax = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
    #             fig.add_axes(cax)
    #             cb = plt.colorbar(xc, cax=cax)
    #
    #             ax1.set_title('{} - no shift'.format(pv), fontsize=18)
    #             ax1.invert_yaxis()
    #
    #             # xc = ax2.scatter(tm, depth, c=oxy_shift, cmap=cmo.cm.oxy, s=10, edgecolor='None')
    #             # if pv == 'ph':
    #             #     xc = ax2.scatter(df_dr.time, df_dr.depth, c=df_dr[colname_shift], cmap=cm.jet, s=10, edgecolor='None',
    #             #                      vmin=7.92, vmax=8.19)
    #             # else:
    #             #     xc = ax2.scatter(df_dr.time, df_dr.depth, c=df_dr[colname_shift], cmap=cm.jet, s=10, edgecolor='None')
    #             xc = ax2.scatter(df_dr.time, df_dr.depth, c=df_dr[colname_shift], cmap=cm.jet, s=10, edgecolor='None')
    #
    #             # format colorbar
    #             divider = make_axes_locatable(ax2)
    #             cax = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
    #             fig.add_axes(cax)
    #             cb = plt.colorbar(xc, cax=cax)
    #
    #             ax2.set_title('{} - {} second shift'.format(pv, sensor_shifts[pv]), fontsize=18)
    #
    #             xfmt = mdates.DateFormatter('%m-%dT%H:%M')
    #             ax2.xaxis.set_major_formatter(xfmt)
    #             ax2.xaxis.set_tick_params(labelsize=15)
    #
    #             sname = os.path.join(sdir, '{}_{}.png'.format(pv, daterange[dt_idx - 1].strftime('%Y%m%dT%H%M')))
    #             plt.savefig(sname, dpi=300)
    #             plt.close()


if __name__ == '__main__':
    deployment = 'ru30-20210226T1647'
    ph_sn = 'sbe10344'
    sensor_shifts = {'ph': 31, 'oxygen': 36}  # sensor_shifts = {'ph': 31, 'oxygen': 36}
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed.nc'
    main(deployment, ph_sn, sensor_shifts, ncfile)
