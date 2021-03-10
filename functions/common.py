#! /usr/bin/env python3

"""
Author: Lori Garzio on 3/8/2021
Last modified: 3/8/2021
"""
import os
import datetime as dt
import glob


def find_calfile(sn):
    """
    Searches the calibration files located in phglider/calibration and returns the most recent calibration file for a
    sensor serial number based on date in file name (sn-%Y%m%d.txt) where sn = serial number and %Y%m%d = calibration
    date.
    :param sn: sensor serial number (e.g. 'sbe10344')
    :return: full file path to the most recent calibration file
    """
    # dirname = os.path.dirname(os.getcwd())
    dirname = '/home/lgarzio/repo/lgarzio/phglider'  # in server
    caldir = os.path.join(dirname, 'calibration')
    calfiles = sorted(glob.glob(caldir + '/{}*.txt'.format(sn)))
    if len(calfiles) > 1:
        datestrs = [cf.split('/')[-1].split('_')[-1].split('.')[0] for cf in calfiles]
        caldates = [dt.datetime.strptime(d, '%Y%m%d') for d in datestrs]
        max_caldate = max(caldates).strftime('%Y%m%d')
        cfile = caldir + '/{}_{}.txt'.format(sn, max_caldate)
    elif len(calfiles) == 1:
        cfile = calfiles[0]
    else:
        raise ValueError('No cal file found for {}'.format(sn))

    return cfile
