#!/usr/bin/env python

"""
Author: Lori Garzio on 4/2/2021
Last modified: 4/2/2021
Download a profile-sci-delayed netCDF dataset and save to a local directory
"""

import os
import functions.common as cf


def main(deploy, sdir):
    ru_server = 'http://slocum-data.marine.rutgers.edu//erddap'
    glider_id = '{}-profile-sci-delayed'.format(deploy)
    glider_vars = ['latitude', 'longitude', 'depth', 'conductivity', 'salinity', 'sci_water_pressure',
                   'temperature', 'sbe41n_ph_ref_voltage', 'chlorophyll_a', 'oxygen_concentration', 'water_depth',
                   'profile_lat', 'profile_lon', 'profile_time']
    ds = cf.get_erddap_dataset(ru_server, glider_id, glider_vars)
    fname = '{}-profile-sci-delayed.nc'.format(deployment)
    ds.to_netcdf(os.path.join(sdir, fname))


if __name__ == '__main__':
    deployment = 'ru30-20210226T1647'
    savedir = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed'
    main(deployment, savedir)
