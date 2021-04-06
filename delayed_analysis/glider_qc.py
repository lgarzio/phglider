#!/usr/bin/env python

"""
Author: Lori Garzio on 4/5/2021
Last modified: 4/5/2021
Modified from MATLAB code writted by Liza Wright-Fairbanks, Spring 2020.
Based on IOOS manual for real-time quality control of pH data observations.
Tests to be applied as a post-processing QA/QC measure.
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def assign_attrs(data_array, data_range, test):
    data_array.attrs['actual_range'] = data_range
    data_array.attrs['comment'] = 'QARTOD {} test'.format(test)
    data_array.attrs['flag_meanings'] = 'GOOD NOT_EVALUATED SUSPECT BAD MISSING'
    data_array.attrs['flag_values'] = '[1 2 3 4 9]'
    data_array.attrs['ioos_category'] = 'Quality'
    data_array.attrs['long_name'] = 'QARTOD {} test'.format(test)
    data_array.attrs['valid_min'] = '1'
    data_array.attrs['valid_max'] = '9'

    return data_array


def gross_range(dataset, vname, user_lims, sensor_lims):
    """
    Test 4: QARTOD test for sensor and user pH value ranges
    :param dataset: xarray dataset
    :param vname: variable name (e.g. ph_total)
    :param user_lims: user-defined minimum and maximum data limits (e.g. for pH min: 6.5, max : 9)
    :param sensor_lims: minimum and maximum sensor limits (e.g. for pH min: 0, max: 14)
    :return: dataset with the pH gross range flag variable added
    """
    val = dataset[vname]
    gross_range_flag = np.full(len(val), 1, dtype='int8')

    # check for nans
    nan_idx = np.isnan(val.values)
    gross_range_flag[nan_idx] = 9

    # user range
    user_idx = np.logical_or(val.values > user_lims['max'], val.values < user_lims['min'])
    gross_range_flag[user_idx] = 3

    # sensor range
    sensor_idx = np.logical_or(val.values >= sensor_lims['max'], val.values <= sensor_lims['min'])
    gross_range_flag[sensor_idx] = 4

    varname = 'qartod_{}_gross_range_flag'.format(val.name)
    da = xr.DataArray(gross_range_flag, coords=val.coords, dims=val.dims, name=varname)
    vmin, vmax = np.nanmin(gross_range_flag), np.nanmax(gross_range_flag)
    if vmin == vmax:
        vrange = '[{}]'.format(vmin)
    else:
        vrange = '[{} {}]'.format(vmin, vmax)
    da = assign_attrs(da, vrange, 'gross range')
    da.attrs['comment'] = '{}. SUSPECT = values < {} or > {}. BAD = values <= {} or >= {}.'.format(da.attrs['comment'],
                                                                                                   user_lims['min'],
                                                                                                   user_lims['max'],
                                                                                                   sensor_lims['min'],
                                                                                                   sensor_lims['max'])
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
        da = assign_attrs(da, vrange, 'location')
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
    da = assign_attrs(da, vrange, 'spike')
    da.attrs['comment'] = '{}. SUSPECT = threshold {}. BAD = threshold {}'.format(da.attrs['comment'],
                                                                                  thresholds['low'],
                                                                                  thresholds['high'])
    dataset[varname] = da

    return dataset


def time_test(da, sdir):
    """
    Test 1: IDs data gaps of pH data >1 hour
    :param da: xarray data array dimensioned by time
    :param sdir: full path to save file directory
    """
    # drop nan values
    da = da.dropna(dim='time')
    time_df = pd.DataFrame(da.time.values, columns=['time'])

    # calculate time difference between each timestamp and add to dictionary
    gap_list = []
    time_df['diff'] = time_df['time'].diff()
    idx_gap = time_df['diff'][time_df['diff'] > pd.Timedelta(hours=1)].index.tolist()
    if len(idx_gap) > 0:
        for i in idx_gap:
            t0 = time_df['time'][i - 1]
            tf = time_df['time'][i]
            gaplen_totalhours = ((tf - t0).days / 24) + np.round((tf - t0).seconds/3600, 2)
            gap_list.append([pd.to_datetime(t0).strftime('%Y-%m-%dT%H:%M:%S'),
                             pd.to_datetime(tf).strftime('%Y-%m-%dT%H:%M:%S'), gaplen_totalhours])
        ttest_df = pd.DataFrame(gap_list, columns=['gap_start', 'gap_end', 'gap_length_hours'])
        ttest_df.to_csv(os.path.join(sdir, '{}_timegaps.csv'.format(da.name)), index=False)

    print('\n{} data gaps >1 hour: {}'.format(da.name, len(gap_list)))


def main(coord_lims, phu_lims, fname):
    ds = xr.open_dataset(fname)
    savedir = os.path.join(os.path.dirname(fname), 'qc')
    os.makedirs(savedir, exist_ok=True)

    # Test 1: QARTOD Time Test
    time_test(ds['ph_total'], savedir)
    time_test(ds['sbe41n_ph_ref_voltage'], savedir)
    if np.sum(~np.isnan(ds.ph_total_shifted.values)) > 0:
        time_test(ds['ph_total_shifted'], savedir)
        time_test(ds['sbe41n_ph_ref_voltage_shifted'], savedir)
    time_test(ds['oxygen_concentration'], savedir)
    if np.sum(~np.isnan(ds.oxygen_concentration_shifted.values)) > 0:
        time_test(ds['oxygen_concentration_shifted'], savedir)
    time_test(ds.chlorophyll_a, savedir)
    time_test(ds.temperature, savedir)
    time_test(ds.salinity, savedir)

    # Test 3: QARTOD Location Test
    ds = location_test(ds, coord_lims)

    # Test 4: QARTOD Gross Range Test
    ph_sensor_lims = {'min': 0, 'max': 14}
    ds = gross_range(ds, 'ph_total', phu_lims, ph_sensor_lims)
    if np.sum(~np.isnan(ds.ph_total_shifted.values)) > 0:
        ds = gross_range(ds, 'ph_total_shifted', phu_lims, ph_sensor_lims)

    # Test 6: QARTOD Spike Test
    ph_thresholds = {'low': 0.1, 'high': 0.2}
    ds = spike_test(ds, 'ph_total', ph_thresholds)
    if np.sum(~np.isnan(ds.ph_total_shifted.values)) > 0:
        ds = spike_test(ds, 'ph_total_shifted', ph_thresholds)

    # Test 7: Rate of Change Test

    ds.to_netcdf('{}_qc.nc'.format(fname.split('.')[0]))


if __name__ == '__main__':
    coordinate_lims = {'latitude': {'min': 36, 'max': 43}, 'longitude': {'min': -76, 'max': -66}}
    ph_user_lims = {'min': 6.5, 'max': 9}
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed_shifted.nc'
    main(coordinate_lims, ph_user_lims, ncfile)
