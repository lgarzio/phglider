#!/usr/bin/env python

"""
Author: Lori Garzio on 4/2/2021
Last modified: 6/16/2022
Download a user-specified pH glider netCDF dataset from RUCOOL's glider ERDDAP server and save to a local directory
"""

import numpy as np
import os
import functions.common as cf


def flatten(lst):
    return [item for sublist in lst for item in sublist]


def main(deploy, sdir):
    os.makedirs(sdir, exist_ok=True)
    ru_server = 'http://slocum-test.marine.rutgers.edu//erddap'

    ds_vars = cf.get_dataset_variables(ru_server, deploy)

    glider_vars = [['latitude', 'longitude', 'depth', 'water_depth', 'profile_lat', 'profile_lon', 'profile_time',
                   'conductivity', 'salinity', 'density', 'pressure', 'temperature']]

    search_str = ['_qartod_summary_flag', 'instrument_', 'ph_ref_voltage', 'oxygen_', 'chlorophyll_a',
                  '_hysteresis_test']
    for ss in search_str:
        append_vars = [x for x in ds_vars if ss in x]
        if search_str == 'chlorophyll_a':
            if len(append_vars) == 0:
                append_vars = ['sci_flntu_chlor_units']
        glider_vars.append(append_vars)

    glider_vars = flatten(glider_vars)
    glider_vars = list(np.unique(glider_vars))

    kwargs = dict()
    kwargs['variables'] = glider_vars
    ds = cf.get_erddap_dataset(ru_server, deploy, **kwargs)
    fname = f'{deploy}.nc'
    ds.to_netcdf(os.path.join(sdir, fname))


if __name__ == '__main__':
    deployment = 'maracoos_02-20220420T2011-profile-sci-delayed'
    savedir = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/maracoos_02-20220420T2011/delayed'
    main(deployment, savedir)
