#!/usr/bin/env python

"""
Author: Lori Garzio on 4/2/2021
Last modified: 8/10/2021
Download a user-specified netCDF dataset from RUCOOL's glider ERDDAP server and save to a local directory
"""

import os
import functions.common as cf


def main(deploy, sdir):
    os.makedirs(sdir, exist_ok=True)
    ru_server = 'http://slocum-data.marine.rutgers.edu//erddap'
    glider_vars = ['latitude', 'longitude', 'depth', 'conductivity', 'salinity', 'sci_water_pressure',
                   'temperature', 'sbe41n_ph_ref_voltage', 'chlorophyll_a', 'oxygen_concentration', 'water_depth',
                   'profile_lat', 'profile_lon', 'profile_time']

    kwargs = dict()
    kwargs['variables'] = glider_vars
    ds = cf.get_erddap_dataset(ru_server, deploy, **kwargs)
    fname = f'{deploy}.nc'
    ds.to_netcdf(os.path.join(sdir, fname))


if __name__ == '__main__':
    deployment = 'ru30-20210503T1929-profile-sci-delayed'
    savedir = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210503T1929/delayed'
    main(deployment, savedir)
