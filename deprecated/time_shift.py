#!/usr/bin/env python

"""
Author: Lori Garzio on 4/2/2021
Last modified: 8/15/2021
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


def main(sensor_sn, shifts, fname):
    deploy = '-'.join(fname.split('/')[-1].split('-')[0:2])
    ds = xr.open_dataset(fname)
    ds = ds.drop_vars(names=['profile_id', 'rowSize'])
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)

    # convert pH voltages of 0.0 to nan
    ds['sbe41n_ph_ref_voltage'][ds['sbe41n_ph_ref_voltage'] == 0.0] = np.nan

    # apply QC from files
    ctd_vars = ['conductivity', 'temperature', 'salinity', 'density', 'instrument_ctd']
    for cv in ctd_vars:
        qc_vars = [x for x in ds.data_vars if np.logical_and(f'{cv}_' in x, x.endswith('_test'))]
        for qv in qc_vars:
            flag_vals = ds[qv].values
            if 'hysteresis' in qv:
                qc_idx = np.where(np.logical_or(flag_vals == 3, flag_vals == 4))[0]
            else:
                qc_idx = np.where(flag_vals == 4)[0]
            if len(qc_idx) > 0:
                if cv == 'instrument_ctd':
                    ds['conductivity'].values[qc_idx] = np.nan
                    ds['temperature'].values[qc_idx] = np.nan
                    ds['salinity'].values[qc_idx] = np.nan
                    ds['density'].values[qc_idx] = np.nan
                else:
                    ds[cv].values[qc_idx] = np.nan  # convert QC'd values to nan

    # read cal file
    calfile = cf.find_calfile(deploy, sensor_sn)
    with open(calfile) as json_file:
        cc = json.load(json_file)

    df = ds.to_dataframe()

    try:
        df['chlorophyll_a']
        chlvar = 'chlorophyll_a'
    except KeyError:
        chlvar = 'sci_flntu_chlor_units'

    # convert profile_time to seconds since 1970-01-01 so you don't lose the column in resampling
    df['profile_time'] = df['profile_time'].astype('int64')//1e9
    df = df.resample('1s').mean()
    col_subset = ['conductivity', 'salinity', 'sci_water_pressure', 'temperature', 'sbe41n_ph_ref_voltage',
                  chlvar, 'oxygen_concentration']
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

    # drop timestamps with duplicated and missing data
    df = df.dropna(axis=0, how='all', subset=col_subset)
    df.drop_duplicates(inplace=True)

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
                    com = 'values shifted by -{} seconds before calculating pH to correct for thermal lag.'.format(sensor_shifts['ph'])
                else:
                    del ds_shifted[varname].attrs['actual_range']
                    com = 'No pH time shift calculated - use sbe41n_ph_ref_voltage variable'
                ds_shifted[varname].attrs['comment'] = com
            if varname == 'oxygen_concentration_shifted':
                ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds['oxygen_concentration'].attrs)
                if np.sum(~np.isnan(ds_shifted[varname].values)) > 0:
                    com = 'values shifted by -{} seconds to correct for thermal lag.'.format(sensor_shifts['oxygen'])
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
                ds_shifted[varname].attrs['calibration'] = str(cc)
                ds_shifted[varname].attrs['units'] = '1'
            if varname == 'ph_total':
                ds_shifted[varname].attrs['observation_type'] = 'calculated'
                ds_shifted[varname].attrs['units'] = '1'
                ds_shifted[varname].attrs['calibration'] = str(cc)
                ds_shifted[varname].attrs['comment'] = 'calculated from instrument calibration coefficents, pressure, salinity, temperature, and measured reference voltage'
            if varname == 'ph_total_shifted':
                ds_shifted[varname].attrs['observation_type'] = 'calculated'
                ds_shifted[varname].attrs['units'] = '1'
                ds_shifted[varname].attrs['calibration'] = str(cc)
                if np.sum(~np.isnan(ds_shifted[varname].values)) > 0:
                    com = 'calculated from instrument calibration coefficents, pressure, salinity, temperature, and ' \
                          'measured reference voltage (shifted -{} seconds).'.format(sensor_shifts['ph'])
                else:  # if there is no pH shift
                    com = 'No pH time shift calculated - use ph_total variable'
                ds_shifted[varname].attrs['comment'] = com

    ds_shifted['time'] = ds_shifted['time'].assign_attrs(ds['time'].attrs)

    ds_shifted = ds_shifted.assign_coords(
        {'latitude': ds_shifted.latitude, 'longitude': ds_shifted.longitude, 'depth': ds_shifted.depth})

    ds_global_attrs = ds.attrs
    ds_shifted = ds_shifted.assign_attrs(ds_global_attrs)
    ds_shifted.to_netcdf('{}_shifted.nc'.format(fname.split('.')[0]))


if __name__ == '__main__':
    ph_sn = 'sbeC18'  # 'sbeC18' 'sbe10344'
    sensor_shifts = {'ph': 23, 'oxygen': 0}  # sensor_shifts = {'ph': 31, 'oxygen': 36}
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210716T1804/delayed/testserver_qc/ru30-20210716T1804-profile-sci-delayed-qc.nc'
    # ph_sn = 'sbe10528'  # 'sbeC18' 'sbe10344'  'sbe10490'
    # sensor_shifts = {'ph': 17, 'oxygen': 0}  # sensor_shifts = {'ph': 31, 'oxygen': 36}
    # ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/sbu01-20210720T1628/delayed/sbu01-20210720T1628-profile-sci-delayed-qc.nc'
    main(ph_sn, sensor_shifts, ncfile)
