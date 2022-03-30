#!/usr/bin/env python

"""
Author: Lori Garzio on 10/6/2021
Last modified: 10/6/2021
Apply multiple time shifts for pH glider data. Before a shift is applied, the dataset is recalculated to
one second averages. The time shifted data are added back to the dataset which is saved as a new netCDF file.
"""

import numpy as np
import pandas as pd
import datetime as dt
import json
import pytz
import xarray as xr
import functions.phcalc as phcalc
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def shift_df(data, vname, seconds_shift):
    da = pd.DataFrame(data[vname])
    tm_shift = data.index - dt.timedelta(seconds=seconds_shift)
    da['time_shift'] = tm_shift
    da.reset_index(drop=True, inplace=True)  # drop the original time index
    da = da.rename(columns={vname: '{}_shifted'.format(vname), 'time_shift': 'time'})
    da = da.set_index('time')

    return da


def main(sensor_sn, ph_shifts, ph_shift_times, fname):
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
    df['profile_time'] = df['profile_time'].astype('int64') / 1e9
    df = df.resample('1s').mean()

    col_subset = ['conductivity', 'salinity', 'density', 'sci_water_pressure', 'temperature', 'sbe41n_ph_ref_voltage',
                  chlvar, 'oxygen_concentration']

    # drop timestamps with duplicated and missing data
    df = df.dropna(axis=0, how='all', subset=col_subset)
    df.drop_duplicates(inplace=True)

    # change profile time back to datetime
    df['profile_time'] = df['profile_time'].map(lambda t: dt.datetime.utcfromtimestamp(t).replace(tzinfo=pytz.UTC))

    # shift each section of the dataset based on times provided
    # select the section of data based on profile time rather than time
    varname = 'sbe41n_ph_ref_voltage'
    shifted_df = pd.DataFrame()
    for i, st in enumerate(ph_shift_times):   ###START HERE
        try:
            st = st.replace(tzinfo=pytz.UTC)
        except TypeError:
            st = st
        if i == 0:
            dfs = df.loc[df.profile_time < st]
        elif i == len(ph_shift_times) - 1:
            dfs = df.loc[df.profile_time > ph_shift_times[i - 1].replace(tzinfo=pytz.UTC)]
        else:
            dfs = df.loc[(df.profile_time > ph_shift_times[i - 1].replace(tzinfo=pytz.UTC)) & (df.profile_time < st)]

        shifted_data = shift_df(dfs, varname, ph_shifts[i])
        shifted_df = shifted_df.append(shifted_data)

    # format shift information for metadata
    ph_shifts = [-x for x in ph_shifts]
    ph_shift_times_str = [x.strftime('%Y-%m-%dT%H:%M:%S') for x in ph_shift_times[0:-1]]

    # merge the shifted dataset back into the original dataframe
    df = df.merge(shifted_df, how='left', left_index=True, right_index=True)
    col_subset.append('{}_shifted'.format(varname))

    # add empty variable for oxygen_concentration_shifted
    da = pd.DataFrame(df['oxygen_concentration'])
    da['oxygen_concentration'] = np.nan
    da = da.rename(columns={'oxygen_concentration': 'oxygen_concentration_shifted'})
    df = df.merge(da, how='left', left_index=True, right_index=True)
    col_subset.append('oxygen_concentration_shifted')

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

    # format profile_time (again)
    ds_shifted['profile_time'].values = np.array(pd.to_datetime(ds_shifted.profile_time.values).strftime('%Y-%m-%dT%H:%M:%S').astype('datetime64[ns]'))

    # assign attributes from the original ds
    for varname, da in ds_shifted.data_vars.items():
        try:
            ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds[varname].attrs)
        except KeyError:
            if varname == 'sbe41n_ph_ref_voltage_shifted':
                ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds['sbe41n_ph_ref_voltage'].attrs)
                if np.sum(~np.isnan(ds_shifted[varname].values)) > 0:
                    com = 'values shifted by {} seconds for time divisions: {} before calculating pH to correct for ' \
                          'thermal lag.'.format(ph_shifts, ph_shift_times_str)
                else:
                    del ds_shifted[varname].attrs['actual_range']
                    com = 'No pH time shift calculated - use sbe41n_ph_ref_voltage variable'
                ds_shifted[varname].attrs['comment'] = com
            if varname == 'oxygen_concentration_shifted':
                ds_shifted[varname] = ds_shifted[varname].assign_attrs(ds['oxygen_concentration'].attrs)
                if np.sum(~np.isnan(ds_shifted[varname].values)) > 0:
                    #com = 'values shifted by -{} seconds to correct for thermal lag.'.format(sensor_shifts['oxygen'])
                    com = 'values shifted by -xx seconds to correct for thermal lag.'
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
                          'measured reference voltage (shifted {} seconds for time divisions: ' \
                          '{}).'.format(ph_shifts, ph_shift_times_str)
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
    ph_shifts = [7, 20, 11, 30]
    ph_times = [dt.datetime(2021, 7, 17, 19, 11, 19), dt.datetime(2021, 7, 22, 19, 58, 57),
                dt.datetime(2021, 7, 27, 11, 41, 45), 'end']
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210716T1804/delayed/testserver_qc/ru30-20210716T1804-profile-sci-delayed-qc.nc'
    main(ph_sn, ph_shifts, ph_times, ncfile)
