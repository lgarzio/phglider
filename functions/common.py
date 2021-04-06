#! /usr/bin/env python3

"""
Author: Lori Garzio on 3/8/2021
Last modified: 4/2/2021
"""
import datetime as dt
import glob
import pandas as pd
import numpy as np
from erddapy import ERDDAP
from sklearn.linear_model import LinearRegression


def find_calfile(deployment, sn):
    """
    Searches the calibration files located in phglider/calibration and returns the most recent calibration file for a
    sensor deployment and serial number based on date in file name (sn-%Y%m%d.txt) where sn = serial number and
    %Y%m%d = calibration date.
    :param deployment: glider deployment, including deployment date (e.g. ru30-20210226T1647)
    :param sn: sensor serial number (e.g. 'sbe10344')
    :return: full file path to the most recent calibration file
    """
    caldir = '/Users/garzio/Documents/repo/lgarzio/phglider/calibration'
    #caldir = '/home/lgarzio/repo/lgarzio/phglider/calibration'  # in server
    calfiles = sorted(glob.glob(caldir + '/{}*.txt'.format(sn)))  # get all cal files for the serial number
    deploy_date = pd.to_datetime(deployment.split('-')[-1])
    if len(calfiles) > 1:
        datestrs = [cf.split('/')[-1].split('_')[-1].split('.')[0] for cf in calfiles]
        caldates = [dt.datetime.strptime(d, '%Y%m%d') for d in datestrs]

        # calculate the difference between the calibration date and the deployment date
        # only keep differences that are positive (ignore cal dates after the deployment date)
        cals = []
        diff = []
        for cd in caldates:
            dd = deploy_date - cd
            if dd.days > 0:
                cals.append(cd)
                diff.append(dd)
        closest_caldate = cals[np.argmin(diff)].strftime('%Y%m%d')  # find the closest date to the deployment date
        cfile = caldir + '/{}_{}.txt'.format(sn, closest_caldate)

    elif len(calfiles) == 1:
        cfile = calfiles[0]
    else:
        raise ValueError('No cal file found for {}'.format(sn))

    return cfile


def get_erddap_dataset(server, ds_id, var_list=None):
    e = ERDDAP(server=server,
               protocol='tabledap',
               response='nc')
    e.dataset_id = ds_id
    if var_list:
        e.variables = var_list
    ds = e.to_xarray()
    ds = ds.sortby(ds.time)
    return ds


def linear_regression(x, y):
    """
    :param x: independent variable
    :param y: dependent variable
    :return: r-squared, slope, intercept, and the predicted y-values
    """

    x = x.reshape((-1, 1))
    model = LinearRegression().fit(x, y)
    r_squared = model.score(x, y)
    slope = model.coef_
    intercept = model.intercept_

    y_pred = model.predict(x)

    return r_squared, slope[0], intercept, y_pred
