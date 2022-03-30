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
