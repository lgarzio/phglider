#!/usr/bin/env python

"""
Author: Lori Garzio on 4/5/2021
Last modified: 8/10/2021
Modified from MATLAB code written by Liza Wright-Fairbanks, Spring 2020.
Based on IOOS manual for real-time quality control of pH data observations.
"""

import numpy as np
import pandas as pd
import xarray as xr
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def assign_qartod_attrs(data_array, data_range, test):
    data_array.attrs['actual_range'] = data_range
    data_array.attrs['comment'] = 'QARTOD {} test'.format(test)
    data_array.attrs['flag_meanings'] = 'GOOD NOT_EVALUATED SUSPECT BAD MISSING'
    data_array.attrs['flag_values'] = '[1 2 3 4 9]'
    data_array.attrs['ioos_category'] = 'Quality'
    data_array.attrs['long_name'] = 'QARTOD {} test'.format(test)
    data_array.attrs['valid_min'] = '1'
    data_array.attrs['valid_max'] = '9'

    return data_array


def bad_profiles():
    """
    Check for bad CTD profile pairs
    """


def gross_range(dataset, vname, sensor_lims, user_lims=None):
    """
    Test 4: QARTOD test for sensor and user pH value ranges
    :param dataset: xarray dataset
    :param vname: variable name (e.g. ph_total)
    :param sensor_lims: minimum and maximum sensor limits (e.g. for pH min: 0, max: 14)
    :param user_lims: optional user-defined minimum and maximum data limits (e.g. for pH min: 6.5, max : 9)
    :return: dataset with the pH gross range flag variable added
    """
    val = dataset[vname]
    gross_range_flag = np.full(len(val), 1, dtype='int8')

    # check for nans
    nan_idx = np.isnan(val.values)
    gross_range_flag[nan_idx] = 9

    # user range
    if user_lims:
        user_idx = np.logical_or(val.values > user_lims['max'], val.values < user_lims['min'])
        gross_range_flag[user_idx] = 3
        add_com_user = 'SUSPECT: values < {} or > {}.'.format(user_lims['min'], user_lims['max'])

    # sensor range
    if np.logical_and(~np.isnan(sensor_lims['min']), ~np.isnan(sensor_lims['max'])):
        sensor_idx = np.logical_or(val.values >= sensor_lims['max'], val.values <= sensor_lims['min'])
        add_com_sen = 'BAD: values <= {} or >= {}.'.format(sensor_lims['min'], sensor_lims['max'])
    elif ~np.isnan(sensor_lims['max']):
        sensor_idx = val.values >= sensor_lims['max']
        add_com_sen = 'BAD: values >= {}.'.format(sensor_lims['max'])
    elif ~np.isnan(sensor_lims['min']):
        sensor_idx = val.values <= sensor_lims['min']
        add_com_sen = 'BAD: values <= {}.'.format(sensor_lims['min'])
    gross_range_flag[sensor_idx] = 4

    varname = 'qartod_{}_gross_range_flag'.format(val.name)
    da = xr.DataArray(gross_range_flag, coords=val.coords, dims=val.dims, name=varname)
    vmin, vmax = np.nanmin(gross_range_flag), np.nanmax(gross_range_flag)
    if vmin == vmax:
        vrange = '[{}]'.format(vmin)
    else:
        vrange = '[{} {}]'.format(vmin, vmax)
    da = assign_qartod_attrs(da, vrange, 'gross range')
    if user_lims:
        com = '{}. {} {}'.format(da.attrs['comment'], add_com_user, add_com_sen)
    else:
        com = '{}. {}'.format(da.attrs['comment'], add_com_sen)
    da.attrs['comment'] = com
    dataset[varname] = da

    return dataset


def location_test(dataset, coordlims):
    """
    Test 3: QARTOD test for the global and local coordinate locations of the latitude and longitude variables
    :param dataset: xarray dataset with latitude and longitude coordinates
    :param coordlims: dictionary containing coordinate limits
    :return: dataset with lat/lon location flag variables added
    """
    global_values = {'latitude': 90, 'longitude': 180}
    locs = ['latitude', 'longitude']
    for loc in locs:
        val = dataset[loc]
        location_flag = np.full(len(val), 1, dtype='int8')

        # check for nans
        nan_idx = np.isnan(val.values)
        location_flag[nan_idx] = 9

        # global check
        gl_idx = np.argwhere(abs(val.values) > global_values[loc])
        location_flag[gl_idx] = 4

        # local check
        idx = np.logical_or(val.values > coordlims[loc]['max'], val.values < coordlims[loc]['min'])
        location_flag[idx] = 3

        # made a data array for the qc flag and add it to the dataset
        varname = 'qartod_{}_location_flag'.format(loc)
        da = xr.DataArray(location_flag, coords=val.coords, dims=val.dims, name=varname)
        vmin, vmax = np.nanmin(location_flag), np.nanmax(location_flag)
        if vmin == vmax:
            vrange = '[{}]'.format(vmin)
        else:
            vrange = '[{} {}]'.format(vmin, vmax)
        da = assign_qartod_attrs(da, vrange, 'location')
        dataset[varname] = da

    return dataset


def spike_test(dataset, vname, thresholds):
    """
    Test 6: QARTOD spike test
    :param dataset: xarray dataset
    :param vname: variable name (e.g. ph_total)
    :param thresholds: user-defined low and high thresholds for spikes
    :return: dataset with lat/lon location flag variables added
    """
    print('\nRunning spike test on {}'.format(vname))
    val = dataset[vname]
    #val = val[0:50]  # for debugging
    spike_flag = np.empty(len(val))
    spike_flag[:] = np.nan

    # TODO - check number of nans in spike_ref

    # define spike reference as mean of 10 surrounding points
    # or if at the very beginning(end) of the array, spike reference is the mean of the 10 subsequent(previous) points
    for i, value in enumerate(val.values):
        #print('{} of {}'.format(i, len(val)-1))
        if i < 5:
            spike_ref = np.nanmean(val.values[(i + 1):(i + 11)])
        elif i > len(val) - 6:
            spike_ref = np.nanmean(val.values[(i - 10):i])
        else:
            spike_ref = np.nanmean(np.append(val.values[(i - 5):i], val.values[(i + 1):(i + 6)]))

        # check the difference vs thresholds
        diff = abs(spike_ref - value)
        if np.isnan(diff):  # check for nans
            qc = 9
        elif diff > thresholds['high']:
            qc = 4
        elif diff > thresholds['low']:
            qc = 3
        else:
            qc = 1
        spike_flag[i] = qc

    spike_flag = spike_flag.astype('int8')
    varname = 'qartod_{}_spike_flag'.format(val.name)
    da = xr.DataArray(spike_flag, coords=val.coords, dims=val.dims, name=varname)
    vmin, vmax = np.nanmin(spike_flag), np.nanmax(spike_flag)
    if vmin == vmax:
        vrange = '[{}]'.format(vmin)
    else:
        vrange = '[{} {}]'.format(vmin, vmax)
    da = assign_qartod_attrs(da, vrange, 'spike')
    da.attrs['comment'] = '{}. SUSPECT = threshold {}. BAD = threshold {}'.format(da.attrs['comment'],
                                                                                  thresholds['low'],
                                                                                  thresholds['high'])
    dataset[varname] = da

    return dataset
