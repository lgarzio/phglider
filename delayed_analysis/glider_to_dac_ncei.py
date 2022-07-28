#!/usr/bin/env python

"""
Author: Lori Garzio on 4/28/2021
Last modified: 7/28/2022
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
import csv
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def make_encoding(dataset, time_start="seconds since 1970-01-01 00:00:00", fillvalue=-9999.0, datatype=np.float32):
    encoding = dict()

    for k in dataset.data_vars:
        if 'time' in k:
            # add the encoding for profile_time so xarray exports the proper time
            encoding[k] = dict(units=time_start, calendar="gregorian", zlib=False, _FillValue=None, dtype=np.double)
        else:
            encoding[k] = dict(zlib=True, _FillValue=np.float32(fillvalue), dtype=datatype)

    # add the encoding for profile_time so xarray exports the proper time
    encoding["time"] = dict(units=time_start, calendar="gregorian", zlib=False, _FillValue=None, dtype=np.double)

    return encoding


def main(fname):
    savedir = os.path.join(os.path.dirname(fname), 'ngdac')
    os.makedirs(savedir, exist_ok=True)

    deploy = '-'.join(fname.split('/')[-1].split('-')[0:2])
    proc_vars = ['time', 'profile_time', 'profile_lon', 'profile_lat', 'latitude', 'longitude', 'depth', 'pressure',
                 'temperature', 'density', 'salinity', 'chlorophyll_a', 'total_alkalinity',
                 'aragonite_saturation_state', 'conductivity',
                 'sbe41n_ph_ref_voltage', 'sbe41n_ph_ref_voltage_shifted', 'ph_total', 'ph_total_shifted',
                 'oxygen_concentration', 'oxygen_concentration_shifted', 'oxygen_saturation',
                 'oxygen_saturation_shifted', 'pressure_interpolated', 'temperature_interpolated',
                 'salinity_interpolated']
    ds = xr.open_dataset(fname)

    # find deployment configuration files
    ga, va, inst = cf.find_configs(deploy)
    with open(va) as json_file:
        var_attrs = json.load(json_file)

    # add additional attribute information to total_alkalinity and instrument_pH
    var_attrs['total_alkalinity']['comment'] = ds.total_alkalinity.comment
    phsensorcals = ds.polynomial_coefficients.calibration_coefficients
    var_attrs['ph_total_shifted']['comment'] = ' sensor cals: '.join((var_attrs['ph_total_shifted']['comment'],
                                                                      phsensorcals))

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

    # add instrument, platform, and constant attributes
    with open(inst) as json_file:
        instruments = json.load(json_file)

    for key, values in instruments.items():
        ds[key] = np.nan
        ds[key].attrs = values

    # add sensor calibration information to instrument_ph
    ds['instrument_ph'].attrs['calibration_coefficients'] = f'pH sensor calibration coefficents: {phsensorcals}'

    # rename variables
    rename_dict = {'latitude': 'lat', 'longitude': 'lon',
                   'sbe41n_ph_ref_voltage': 'pH_reference_voltage_raw',
                   'sbe41n_ph_ref_voltage_shifted': 'pH_reference_voltage_corrected',
                   'ph_total': 'pH_raw',
                   'ph_total_shifted': 'pH_corrected',
                   'oxygen_concentration': 'oxygen_concentration_raw',
                   'oxygen_concentration_shifted': 'oxygen_concentration_corrected',
                   'oxygen_saturation': 'oxygen_saturation_raw',
                   'oxygen_saturation_shifted': 'oxygen_saturation_corrected'}
    ds = ds.rename(rename_dict)

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

    # export tab-delimited lonlat.txt file NCEI data submission
    lon = list(np.round(ds.lon.values, 4))
    lat = list(np.round(ds.lat.values, 4))
    txt_file = os.path.join(savedir, f'lonlat-{deploy}.txt')
    with open(txt_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['lon', 'lat'], delimiter='\t')
        for a, b in zip(lon, lat):
            writer.writerow({'lon': a, 'lat': b})

    savefile = os.path.join(savedir, f'{fname.split("/")[-1].split("_qc.nc")[0]}_ncei.nc')
    ds.to_netcdf(savefile, encoding=var_encoding, format="netCDF4", engine="netcdf4", unlimited_dims=["time"])



if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed_qc.nc'
    main(ncfile)
