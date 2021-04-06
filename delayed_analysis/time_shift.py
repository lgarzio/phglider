#!/usr/bin/env python

"""
Author: Lori Garzio on 4/2/2021
Last modified: 4/3/2021
Apply a best time shift for pH and oxygen glider data. Before a shift is applied, the dataset is recalculated to
one second averages. The time shifted data are added back to the dataset which is saved as a new netCDF file.
"""

import numpy as np
import pandas as pd
import datetime as dt
import json
import xarray as xr
import functions.phcalc as phcalc
import functions.common as cf
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
            if key == 'ph':
                varname = 'sbe41n_ph_ref_voltage'
            elif key == 'oxygen':
                varname = 'oxygen_concentration'
            da = pd.DataFrame(df[varname])
            if val > 0:  # if there is a time shift calculated
                tm_shift = df.index - dt.timedelta(seconds=val)
                da['time_shift'] = tm_shift
                da.reset_index(drop=True, inplace=True)  # drop the original time index
                da = da.rename(columns={varname: '{}_shifted'.format(varname), 'time_shift': 'time'})
                da = da.set_index('time')
            else:  # if there is no time shift calculated, convert the array to nans
                da[varname] = np.nan
                da = da.rename(columns={varname: '{}_shifted'.format(varname)})
            df = df.merge(da, how='left', left_index=True, right_index=True)
            col_subset.append('{}_shifted'.format(varname))

    df = df.dropna(axis=0, how='all', subset=col_subset)

    # calculate pressure_dbar and drop any rows where pressure is <1 dbar
    df['sci_water_pressure_dbar'] = df['sci_water_pressure'] * 10
    df = df[df.sci_water_pressure_dbar > 1]

    # calculate pH and add to dataframe
    df['f_p'] = np.polyval([cc['f6'], cc['f5'], cc['f4'], cc['f3'], cc['f2'], cc['f1'], 0], df.sci_water_pressure_dbar)
    phfree, phtot = phcalc.phcalc(df.sbe41n_ph_ref_voltage, df.sci_water_pressure_dbar, df.temperature, df.salinity,
                                 cc['k0'], cc['k2'], df.f_p)
    df['ph_total'] = phtot

    if np.sum(~np.isnan(df['sbe41n_ph_ref_voltage_shifted'].values)) > 0:
        phfree_sh, phtot_sh = phcalc.phcalc(df['sbe41n_ph_ref_voltage_shifted'], df.sci_water_pressure_dbar,
                                            df.temperature, df.salinity, cc['k0'], cc['k2'], df.f_p)
    else:
        ph_empty = np.array(np.empty(len(phtot)))
        ph_empty[:] = np.nan
        phtot_sh = ph_empty

    # add shifted pH (calculated) to dataframe
    df['ph_total_shifted'] = phtot_sh

    # convert the dataframe to an xarray dataset
    ds_shifted = df.to_xarray()

    # assign attributes from the original ds
    for varname, da in ds_shifted.data_vars.items():
        try:
            ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds[varname].attrs)
        except KeyError:
            if varname == 'sbe41n_ph_ref_voltage_shifted':
                ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds['sbe41n_ph_ref_voltage'].attrs)
                if np.sum(~np.isnan(ds_shifted[varname].values)) > 0:
                    com = 'Values time shifted by -{} seconds to correct for thermal lag.'.format(sensor_shifts['ph'])
                else:
                    del ds_shifted[varname].attrs['actual_range']
                    com = 'No pH time shift calculated - use sbe41n_ph_ref_voltage variable'
                ds_shifted[varname].attrs['comment'] = com
            if varname == 'oxygen_concentration_shifted':
                ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds['oxygen_concentration'].attrs)
                if np.sum(~np.isnan(ds_shifted[varname].values)) > 0:
                    com = 'Values time shifted by -{} seconds to correct for thermal lag.'.format(sensor_shifts['oxygen'])
                else:
                    com = 'No oxygen time shift calculated - use oxygen_concentration variable'
                ds_shifted[varname].attrs['comment'] = com
            if 'sci_water_pressure' in varname:
                ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds['sci_water_pressure'].attrs)
                for item in ['actual_range', 'colorBarMaximum', 'colorBarMinimum', 'valid_max', 'valid_min']:
                    del ds_shifted[varname].attrs[item]
                ds_shifted[varname].attrs['observation_type'] = 'calculated'
                ds_shifted[varname].attrs['units'] = 'dbar'
            if varname == 'f_p':
                ds_shifted[varname].attrs['observation_type'] = 'calculated'
                ds_shifted[varname].attrs['comment'] = 'Polynomial evaluation of sensor pressure response polynomial coefficients (f6-f1) and pressure_dbar'
                ds_shifted[varname].attrs['f6'] = cc['f6']
                ds_shifted[varname].attrs['f5'] = cc['f5']
                ds_shifted[varname].attrs['f4'] = cc['f4']
                ds_shifted[varname].attrs['f3'] = cc['f3']
                ds_shifted[varname].attrs['f2'] = cc['f2']
                ds_shifted[varname].attrs['f1'] = cc['f1']
                ds_shifted[varname].attrs['units'] = '1'
            if varname == 'ph_total':
                ds_shifted[varname].attrs['observation_type'] = 'calculated'
                ds_shifted[varname].attrs['units'] = '1'
                ds_shifted[varname].attrs['k0'] = cc['k0']
                ds_shifted[varname].attrs['k2'] = cc['k2']
                ds_shifted[varname].attrs['comment'] = 'pH calculated on the total scale using pH reference voltage, pressure_dbar ,' \
                                                       'temperature, salinity, sensor calibration coefficients, and pressure ' \
                                                       'polynomial coefficients.'
            if varname == 'ph_total_shifted':
                ds_shifted[varname].attrs['observation_type'] = 'calculated'
                ds_shifted[varname].attrs['units'] = '1'
                ds_shifted[varname].attrs['k0'] = cc['k0']
                ds_shifted[varname].attrs['k2'] = cc['k2']
                if np.sum(~np.isnan(ds_shifted[varname].values)) > 0:
                    com = 'pH calculated on the total scale using time shifted pH reference voltage, pressure_dbar, ' \
                          'temperature, salinity, sensor calibration coefficients, and pressure polynomial coefficients.'
                else:  # if there is no pH shift
                    com = 'No pH time shift calculated - use ph_total variable'
                ds_shifted[varname].attrs['comment'] = com

    ds_shifted = ds_shifted.assign_coords(
        {'latitude': ds_shifted.latitude, 'longitude': ds_shifted.longitude, 'depth': ds_shifted.depth})

    ds_global_attrs = ds.attrs
    ds_shifted = ds_shifted.assign_attrs(ds_global_attrs)
    ds_shifted.to_netcdf('{}_shifted.nc'.format(fname.split('.')[0]))


if __name__ == '__main__':
    deployment = 'ru30-20210226T1647'
    ph_sn = 'sbe10344'
    sensor_shifts = {'ph': 31, 'oxygen': 36}  # sensor_shifts = {'ph': 31, 'oxygen': 36}
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed.nc'
    main(deployment, ph_sn, sensor_shifts, ncfile)
