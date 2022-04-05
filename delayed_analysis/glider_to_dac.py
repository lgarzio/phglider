#!/usr/bin/env python

"""
Author: Lori Garzio on 4/28/2021
Last modified: 4/5/2021
Process final glider dataset to upload to the IOOS glider DAC (https://gliders.ioos.us/).
Modified from code written by Leila Belabbassi.
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import json
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def main(fname):
    deploy = '-'.join(fname.split('/')[-1].split('-')[0:2])
    proc_vars = ['unix_time', 'pressure', 'temperature', 'salinity', 'chlorophyll_a',
                 'total_alkalinity', 'saturation_aragonite', 'latitude', 'longitude', 'conductivity',
                 'sbe41n_ph_ref_voltage_shifted', 'ph_total_shifted', 'oxygen_concentration_shifted',
                 'oxygen_concentration_shifted_mgL']
    ds = xr.open_dataset(fname)

    # add unix time
    dates = pd.to_datetime(ds.time.values)
    unix = (dates - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')
    da = xr.DataArray(unix, coords=ds.time.coords, dims=ds.time.dims, name='unix_time')
    ds['unix_time'] = da

    # find deployment configuration files
    ga, va, inst = cf.find_configs(deploy)
    with open(va) as json_file:
        var_attrs = json.load(json_file)

    # add additional attribute information to total_alkalinity and pH
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
    data_dict.update({'coords': {'time': {'dims': 'time', 'data': ds['unix_time'].values}}})

    # convert to dataset
    ds = xr.Dataset.from_dict(data_dict)

    # add instrument, platform, and constant attributes
    with open(inst) as json_file:
        instruments = json.load(json_file)

    for key, values in instruments.items():
        ds[key] = np.nan
        ds[key].attrs = values

    # rename variables so they match the previous versions that were sent to the DAC
    rename_dict = {'unix_time': 'UnixTime', 'pressure': 'Pressure', 'temperature': 'Temperature',
                   'salinity': 'Salinity', 'chlorophyll_a': 'Chlorophyll', 'total_alkalinity': 'TotalAlkalinity',
                   'saturation_aragonite': 'AragoniteSaturationState', 'latitude': 'Latitude',
                   'longitude': 'Longitude', 'conductivity': 'Conductivity',
                   'sbe41n_ph_ref_voltage_shifted': 'pHReferenceVoltage', 'ph_total_shifted': 'pH',
                   'oxygen_concentration_shifted': 'Oxygen_molar', 'oxygen_concentration_shifted_mgL': 'Oxygen_mgL'}
    ds = ds.rename(rename_dict)

    nc_filename = '{}/{}_to_dac.nc'.format(os.path.dirname(fname), deploy)
    ds.to_netcdf(nc_filename)


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed_qc.nc'
    main(ncfile)
