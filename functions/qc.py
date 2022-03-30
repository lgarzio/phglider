#!/usr/bin/env python

"""
Author: Lori Garzio on 4/5/2021
Last modified: 3/29/2022
Modified from MATLAB code written by Liza Wright-Fairbanks, Spring 2020.
Based on IOOS manual for real-time quality control of pH data observations.
"""

import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def assign_qartod_attrs(test, varname, thresholds):
    """
    Define the QARTOD QC variable attributes
    :param test: QARTOD QC test
    :param varname: sensor variable name (e.g. conductivity)
    :param thresholds: flag thresholds from QC configuration files
    """

    flag_meanings = 'GOOD NOT_EVALUATED SUSPECT FAIL MISSING'
    flag_values = [1, 2, 3, 4, 9]
    standard_name = f'{test}_quality_flag'  # 'gross_range_test_quality_flag'
    long_name = f'{" ".join([x.capitalize() for x in test.split("_")])} Quality Flag'

    # Defining gross/flatline QC variable attributes
    attrs = {
        'standard_name': standard_name,
        'long_name': long_name,
        'flag_values': np.byte(flag_values),
        'flag_meanings': flag_meanings,
        'flag_configurations': str(thresholds),
        'valid_min': np.byte(min(flag_values)),
        'valid_max': np.byte(max(flag_values)),
        'ioos_qc_module': 'qartod',
        'ioos_qc_test': f'{test}',
        'ioos_qc_target': varname,
    }

    return attrs


def location_test(dataset, coordlims):
    """
    Test 3: QARTOD test for the global and local coordinate locations of the latitude and longitude variables
    :param dataset: xarray dataset with latitude and longitude coordinates
    :param coordlims: dictionary containing coordinate limits
    :return: dataset with lat/lon location flag variables added
    """
    global_lims = {'latitude': [-90, 90], 'longitude': [-180, 180]}
    global_values = dict(latitude=[global_lims['latitude'][0], global_lims['latitude'][0],
                                   global_lims['latitude'][1], global_lims['latitude'][1]],
                         longitude=[global_lims['longitude'][0], global_lims['longitude'][1],
                                    global_lims['longitude'][1], global_lims['longitude'][0]])
    global_polygon = Polygon(list(zip(global_values['longitude'], global_values['latitude'])))

    local_values = dict(latitude=[coordlims['latitude'][0], coordlims['latitude'][0],
                                  coordlims['latitude'][1], coordlims['latitude'][1]],
                        longitude=[coordlims['longitude'][0], coordlims['longitude'][1],
                                   coordlims['longitude'][1], coordlims['longitude'][0]])
    local_polygon = Polygon(list(zip(local_values['longitude'], local_values['latitude'])))

    flag_vals = np.full(len(dataset.longitude), 2, dtype='int32')

    # check for missing values
    nan_idx = np.logical_or(np.isnan(dataset.longitude.values), np.isnan(dataset.latitude.values))
    flag_vals[nan_idx] = 9

    # make sure each lat/lon pair is within global and local ranges
    for i, lon in enumerate(dataset.longitude.values):
        lat = dataset.latitude.values[i]

        if not global_polygon.contains(Point(lon, lat)):
            flag_vals[i] = 4
        elif not local_polygon.contains(Point(lon, lat)):
            flag_vals[i] = 3
        else:
            flag_vals[i] = 1

    # made a data array for the qc flag and add it to the dataset
    varname = 'location_qartod_test'
    vmin, vmax = np.nanmin(flag_vals), np.nanmax(flag_vals)
    if vmin == vmax:
        vrange = '[{}]'.format(vmin)
    else:
        vrange = '[{} {}]'.format(vmin, vmax)

    thresholds = dict(suspect=coordlims,
                      fail=global_lims)

    attrs = assign_qartod_attrs('location', 'latitude longitude', thresholds)

    da = xr.DataArray(flag_vals, coords=dataset.longitude.coords, dims=dataset.longitude.dims, name=varname,
                      attrs=attrs)
    da.attrs['ancillary_variables'] = 'latitude longitude'
    da.attrs['actual_range'] = vrange

    dataset[varname] = da

    return dataset


def spike_test(dataset, vname, thresholds):
    """
    Test 6: QARTOD spike test
    :param dataset: xarray dataset
    :param vname: variable name (e.g. ph_total)
    :param thresholds: user-defined suspect and fail thresholds for spikes, and time threshold in seconds
    :return: dataset with lat/lon location flag variables added
    """
    print('\nRunning spike test on {}'.format(vname))
    qc_varname = f'{vname}_qartod_spike_test'

    # grab data for sensor
    data = dataset[vname]
    # identify where not nan
    non_nan_ind = np.invert(np.isnan(data))
    # non_nan_ind[50:] = False  # for debugging

    data_nonan = data[non_nan_ind]

    flag_vals = 2 * np.ones(np.shape(data))
    flag_vals[np.invert(non_nan_ind)] = 9

    # define spike reference as median of 10 surrounding points
    # or if at the very beginning(end) of the array, spike reference is the median of the 10 subsequent(previous) points
    test_vals = 2 * np.ones(np.shape(data_nonan))
    for i, value in enumerate(data_nonan):
        if i < 5:
            data_ref = data_nonan[(i + 1):(i + 11)]
        elif i > len(data_nonan) - 6:
            data_ref = data_nonan[(i - 10):i]
        else:
            data_ref = xr.concat([data_nonan[(i - 5):i], data_nonan[(i + 1):(i + 6)]], dim='time')

        time_diff_seconds = np.diff(data_ref.time).astype('timedelta64[s]').astype(float)
        if np.max(time_diff_seconds) > thresholds['seconds']:
            # if the time interval is longer than the defined threshold,
            # try taking the median of the 4 surrounding data points
            if i < 2:
                data_ref = data_nonan[(i + 1):(i + 5)]
            elif i > len(data_nonan) - 3:
                data_ref = data_nonan[(i - 4):i]
            else:
                data_ref = xr.concat([data_nonan[(i - 2):i], data_nonan[(i + 1):(i + 3)]], dim='time')

            time_diff_seconds = np.diff(data_ref.time).astype('timedelta64[s]').astype(float)
            if np.max(time_diff_seconds) > thresholds['seconds']:
                # if the time interval is still longer than the defined threshold,
                # try taking the median of the 4 subsequent points
                data_ref = data_nonan[(i + 1):(i + 5)]
                time_diff_seconds = np.diff(data_ref.time).astype('timedelta64[s]').astype(float)
                if np.max(time_diff_seconds) > thresholds['seconds']:
                    # try taking the median of the 4 subsequent or previous points
                    data_ref = data_nonan[(i - 4):i]
                    time_diff_seconds = np.diff(data_ref.time).astype('timedelta64[s]').astype(float)
                    if np.max(time_diff_seconds) > thresholds['seconds']:
                        # if the time interval is still longer than the defined threshold, leave flag values as
                        # 2/NOT_EVALUATED
                        continue

        spike_ref = np.median(data_ref)

        # check the difference vs thresholds
        diff = abs(spike_ref - value.values)
        if diff > thresholds['fail']:
            qc_value = 4
        elif diff > thresholds['suspect']:
            qc_value = 3
        else:
            qc_value = 1
        test_vals[i] = qc_value

    # insert test results in the original data array with missing values
    flag_vals[non_nan_ind] = test_vals

    attrs = assign_qartod_attrs('spike_test', vname, thresholds)
    da = xr.DataArray(flag_vals.astype('int32'), coords=data.coords, dims=data.dims, attrs=attrs, name=qc_varname)

    return da
