#!/usr/bin/env python

"""
Author: Lori Garzio on 3/2/2021
Last modified: 3/2/2021
Plot realtime glider data
"""

import argparse
import numpy as np
import json
from erddapy import ERDDAP
import functions.phcalc as phcalc


def get_erddap_dataset(server, protocol, file_type, ds_id, var_list=None):
    e = ERDDAP(server=server,
               protocol=protocol,
               response=file_type)
    e.dataset_id = ds_id
    if var_list:
        e.variables = var_list
    ds = e.to_xarray()
    ds = ds.sortby(ds.time)
    return ds


def main(deploy, calfile):
    ru_server = 'http://slocum-data.marine.rutgers.edu//erddap'
    id = '{}-profile-sci-rt'.format(deploy)
    glider_vars = ['latitude', 'longitude', 'depth', 'conductivity', 'salinity', 'sci_water_pressure',
                   'temperature', 'sbe41n_ph_ref_voltage', 'chlorophyll_a', 'oxygen_concentration', 'water_depth']
    ds = get_erddap_dataset(ru_server, 'tabledap', 'nc', id, glider_vars)
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)

    # read cal file
    with open(calfile) as json_file:
        caldata = json.load(json_file)

    temp = ds.temperature
    pressure = ds.sci_water_pressure
    vrs = ds.sbe41n_ph_ref_voltage
    sal = ds.salinity
    oxy = (ds.oxygen_concentration * 32) / 1000  # change oxygen from umol/L to mg/L

    ph = np.array([])
    for i, press in enumerate(pressure):
        f_P = np.polyval([caldata['f6'], caldata['f5'], caldata['f4'], caldata['f3'], caldata['f2'], caldata['f1'], 0],
                         press)
        phfree, phtot = phcalc.phcalc(vrs.values[i], press.values, temp.values[i], sal.values[i], caldata['k0'],
                                      caldata['k2'], f_P)
        ph = np.append(ph, phtot)


if __name__ == '__main__':
    deployment = 'ru30-20210226T1647'
    cal = '/Users/lgarzio/Documents/repo/lgarzio/phglider/calibration/sbe10344_20200306.txt'
    main(deployment, cal)
    # arg_parser = argparse.ArgumentParser(description='Plot real time glider pH data',
    #                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #
    # arg_parser.add_argument('-deployment',
    #                         dest='deploy',
    #                         default='ru30-20210226T1647',
    #                         type=str,
    #                         help='Glider deployment to plot. e.g. glider-yyyymmddTHHMM')
    #
    # arg_parser.add_argument('-calfile',
    #                         dest='calfile',
    #                         default='/Users/lgarzio/Documents/repo/lgarzio/phglider/calibration/sbe10344_20200306.txt',
    #                         type=str,
    #                         help='Calibration file for SBE pH sensor.')
    #
    # parsed_args = arg_parser.parse_args()
    # sys.exit(main(parsed_args))
