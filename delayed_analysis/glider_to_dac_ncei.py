#!/usr/bin/env python

"""
Author: Lori Garzio on 4/28/2021
Last modified: 5/17/2023
Process final glider dataset to upload to the IOOS glider DAC (https://gliders.ioos.us/) and
NCEI OA data portal (https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system-portal/)
Modified from code written by Leila Belabbassi.
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import json
import ast
import csv
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def make_encoding(dataset, time_start="seconds since 1970-01-01 00:00:00", fillvalue=-9999.0, datatype=np.float32):
    encoding = dict()

    for k in dataset.data_vars:
        if 'profile_time' in k:
            # add the encoding for profile_time so xarray exports the proper time
            encoding[k] = dict(units=time_start, calendar="gregorian", zlib=False, _FillValue=-9999.0, dtype=np.double)
        elif 'profile_id' in k:
            encoding[k] = dict(_FillValue=2147483647, dtype=np.int32)
        elif 'trajectory' in k:
            encoding[k] = dict(dtype=object)
        else:
            encoding[k] = dict(zlib=True, _FillValue=np.float32(fillvalue), dtype=datatype)

    # add the encoding for time so xarray exports the proper time
    encoding["time"] = dict(units=time_start, calendar="gregorian", zlib=False, _FillValue=None, dtype=np.double)

    return encoding


def main(fname):
    savedir = os.path.join(os.path.dirname(fname), 'ngdac')
    os.makedirs(savedir, exist_ok=True)

    deploy = '-'.join(fname.split('/')[-1].split('-')[0:2])
    ds = xr.open_dataset(fname)

    proc_vars = ['time', 'profile_time', 'profile_lon', 'profile_lat', 'latitude', 'longitude', 'depth', 'pressure',
                 'temperature', 'density', 'salinity', 'conductivity',
                 'sbe41n_ph_ref_voltage', 'sbe41n_ph_ref_voltage_shifted', 'ph_total', 'ph_total_shifted',
                 'pressure_interpolated', 'temperature_interpolated', 'salinity_interpolated', 'depth_interpolated']

    search_vars = ['chlorophyll_a', 'sci_flntu_chlor_units', 'total_alkalinity', 'aragonite_saturation_state',
                   'beta_700nm', 'cdom', 'oxygen_concentration', 'oxygen_concentration_shifted', 'oxygen_saturation',
                   'oxygen_saturation_shifted', 'conductivity_combined', 'temperature_combined', 'salinity_combined',
                   'density_combined']
    for sv in search_vars:
        if sv in ds.data_vars:
            proc_vars.append(sv)

    proc_vars = list(np.unique(proc_vars))

    # get the pH sensor calibration coefficients
    phsensorcals = ast.literal_eval(ds.polynomial_coefficients.calibration_coefficients)

    # find deployment configuration files
    ga, va, inst = cf.find_configs(deploy)
    with open(va) as json_file:
        var_attrs = json.load(json_file)

    # add additional attribute information to total_alkalinity
    try:
        var_attrs['total_alkalinity']['comment'] = ds.total_alkalinity.comment
    except AttributeError:
        print('no TA in file')

    # build a dictionary with the data and appropriate variable attributes
    data_dict = {'data_vars': {}}
    for key in proc_vars:
        data_dict['data_vars'].update({key: {'dims': 'time',
                                             'data': ds[key].values,
                                             'attrs': var_attrs[key]}})

    # add global attributes
    with open(ga) as json_file:
        global_attrs = json.load(json_file)

    global_attrs['title'] = deploy
    global_attrs['format_version'] = deploy
    tnow = dt.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:00Z')
    global_attrs['date_created'] = tnow
    global_attrs['date_issued'] = tnow
    global_attrs['date_modified'] = tnow

    sorted_global_attributes = {k: v for k, v in sorted(global_attrs.items(), key=lambda item: item[1])}
    data_dict.update({'attrs': sorted_global_attributes})

    # add dimensions
    data_dict.update({'dims': {'time': {'dims': 'time', 'data': len(proc_vars) - 1}}})

    # add coordinates
    data_dict.update({'coords': {'time': {'dims': 'time', 'data': ds['time'].values}}})

    # convert to dataset
    ds = xr.Dataset.from_dict(data_dict)

    # add profile_id
    attributes = dict(
        ancillary_variables='profile_time',
        cf_role='profile_id',
        comment='Unique identifier of the profile. The profile ID is the mean profile timestamp.',
        ioos_category='Identifier',
        long_name='Profile ID'
    )
    name = 'profile_id'
    pid = ds.profile_time.values.astype('datetime64[s]').astype('int')
    da = xr.DataArray(pid, coords=ds.profile_time.coords, dims=ds.profile_time.dims, name=name, attrs=attributes)
    ds[name] = da

    # add trajectory
    attributes = dict(
        cf_role='trajectory_id',
        comment='A trajectory is a single deployment of a glider and may span multiple data files.',
        ioos_category='Identifier',
        long_name='Trajectory/Deployment Name'
    )
    name = 'trajectory'
    da = xr.DataArray(np.repeat(ds.title, len(ds.time)).astype(object), coords=ds.profile_time.coords,
                      dims=ds.profile_time.dims, name=name, attrs=attributes)
    ds[name] = da

    # add instrument, platform, and constant attributes
    with open(inst) as json_file:
        instruments = json.load(json_file)

    for key, values in instruments.items():
        ds[key] = np.nan
        ds[key].attrs = values

    # add sensor calibration information to instrument_ph
    ph_cal_attrs = 'pH sensor calibration coefficents: '
    for k, v in phsensorcals.items():
        ph_cal_attrs += f'{k}: {v}, '

    # format string to remove final comma-space
    ph_cal_attrs = ph_cal_attrs[0:len(ph_cal_attrs)-2]

    ds['instrument_ph'].attrs['calibration_coefficients'] = ph_cal_attrs

    # rename variables
    rename_dict = {'latitude': 'lat', 'longitude': 'lon',
                   'sbe41n_ph_ref_voltage': 'pH_reference_voltage_raw',
                   'sbe41n_ph_ref_voltage_shifted': 'pH_reference_voltage_corrected',
                   'ph_total': 'pH_raw',
                   'ph_total_shifted': 'pH_corrected',
                   'oxygen_concentration': 'oxygen_concentration_raw',
                   'oxygen_concentration_shifted': 'oxygen_concentration_corrected',
                   'oxygen_saturation': 'oxygen_saturation_raw',
                   'oxygen_saturation_shifted': 'oxygen_saturation_corrected',
                   'conductivity_combined': 'conductivity_lag_shifted',
                   'temperature_combined': 'temperature_lag_shifted',
                   'salinity_combined': 'salinity_lag_shifted',
                   'density_combined': 'density_lag_shifted'
                   }
    ds = ds.rename(rename_dict)

    # make sure valid_min and valid_max are the same data type as the variables
    for v in ds.data_vars:
        try:
            ds[v].attrs['valid_min'] = np.float32(ds[v].valid_min)
            ds[v].attrs['valid_max'] = np.float32(ds[v].valid_max)
        except AttributeError:
            continue

    # specify the variable encoding
    var_encoding = make_encoding(ds)

    # split up the dataset into profile files
    ds_dac = ds.swap_dims({'time': 'profile_time'})
    for pt in np.unique(ds_dac.profile_time.values):
        dspt = ds_dac.sel(profile_time=pt)
        try:
            dspt = dspt.swap_dims({'profile_time': 'time'})
            dspt = dspt.reset_coords(names='profile_time')
            nc_fname = '{}_{}_delayed.nc'.format(deploy.split('-')[0], pd.to_datetime(pt).strftime('%Y%m%dT%H%MZ'))
            nc_filename = os.path.join(savedir, nc_fname)
            dspt.to_netcdf(nc_filename, encoding=var_encoding, format="netCDF4", engine="netcdf4",
                           unlimited_dims=["time"])
        except ValueError:
            continue

    # format file for NCEI submission
    savedir = os.path.join(os.path.dirname(fname), 'ncei')
    os.makedirs(savedir, exist_ok=True)

    # export unique lonlat.csv file NCEI data submission
    lonlat = dict(lon=list(np.round(ds.lon.values, 4)),
                  lat=list(np.round(ds.lat.values, 4)))
    df = pd.DataFrame(lonlat)
    df.drop_duplicates(inplace=True)

    df.to_csv(os.path.join(savedir, f'lonlat-{deploy}.csv'), header=None, index=None)

    savefile = os.path.join(savedir, f'{deploy}-delayed.nc')
    ds.to_netcdf(savefile, encoding=var_encoding, format="netCDF4", engine="netcdf4", unlimited_dims=["time"])


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed_qc.nc'
    main(ncfile)
