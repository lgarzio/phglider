#! /usr/bin/env python3

"""
Author: Lori Garzio on 3/8/2021
Last modified: 3/30/2022
"""
import datetime as dt
import glob
import json
import pandas as pd
import numpy as np
import PyCO2SYS as pyco2
from erddapy import ERDDAP
from sklearn.linear_model import LinearRegression
import cmocean as cmo


def attribute_mapping(varname):
    """
    Map new variable to old for attribute assignment
    :param varname: variable name
    :return: variable name from original dataset to use for attributes with additional relevant attributes
    """
    variable_mapping = dict(
        pressure_interpolated=dict(
            comment='linear interpolated pressure',
            ancillary_variables='pressure'
        ),
        temperature_interpolated=dict(
            comment='linear interpolated temperature',
            ancillary_variables='temperature'
        ),
        salinity_interpolated=dict(
            comment='linear interpolated salinity',
            ancillary_variables='salinity'
        )
    )

    return variable_mapping[varname]


def assign_ph_attrs(varname):
    """
    Assign attributes for pH variables
    :param varname: variable name
    :return: attributes for pH variables
    """
    ph_attributes = dict(
        f_p=dict(
            observation_type='calculated',
            comment='Polynomial evaluation of sensor pressure response polynomial coefficients (f6-f1) and pressure',
            ancillary_variables='pressure_interpolated, polynomial_coefficients',
            units='1'
        ),
        ph_total=dict(
            observation_type='calculated',
            comment='Calculated from instrument calibration coefficents, interpolated pressure, salinity, temperature, '
                    'and measured reference voltage',
            ancillary_variables='sbe41n_ph_ref_voltage pressure_interpolated temperature_interpolated '
                                'salinity_interpolated f_p k0 k2',
            units='1'
        ),
        ph_total_shifted=dict(
            observation_type='calculated',
            comment='Calculated from instrument calibration coefficents, interpolated pressure, salinity, temperature, '
                    'and measured reference voltage shifted by values defined in sbe41n_ph_ref_voltage_optimal_shift',
            ancillary_variables='sbe41n_ph_ref_voltage_shifted pressure_interpolated temperature_interpolated '
                                'salinity_interpolated f_p k0 k2',
            units='1'
        )
    )

    return ph_attributes[varname]


def calculate_ta(deployment, salinity):
    """
    Calculate Total Alkalinity from salinity using a linear relationship determined from in-situ water sampling data
    taken during glider deployment and recovery. See ../ta_equation/ta_sal_regression.py
    :param deployment: glider deployment, including deployment date (e.g. ru30-20210226T1647)
    :param salinity: data array of salinity data from glider
    :return: data array of calculated Total Alkalinity
    """
    # find TA-salinity equation file calculates from ta_sal_regression.py
    tadir = '/Users/garzio/Documents/repo/lgarzio/phglider/ta_equation'
    #tadir = '/home/lgarzio/repo/lgarzio/phglider/ta_equation'  # in server
    tafiles = sorted(glob.glob(tadir + '/{}_ta_equation.txt'.format(deployment)))
    if len(tafiles) > 1:
        raise ValueError('More than 1 TA-salinity equation file found for deployment: {}'.format(deployment))
    elif len(tafiles) < 1:
        raise ValueError('No TA-salinity equation file found for deployment: {}'.format(deployment))

    ta_equ_file = tafiles[0]
    with open(ta_equ_file) as json_file:
        ta_equ = json.load(json_file)

    ta = ta_equ['m'] * salinity + ta_equ['b']
    ta.name = 'total_alkalinity'
    ta.attrs['units'] = 'umol/kg'
    ta.attrs['long_name'] = 'Total Alkalinity'
    ta.attrs['ancillary_variables'] = 'salinity'
    ta.attrs['observation_type'] = 'calculated'
    ta.attrs['comment'] = 'Calculated from salinity using the linear relationship (TA = {} * salinity + {}) ' \
                          'determined from in-situ water sampling data taken during glider deployment and ' \
                          'recovery in addition to ship-based water samples as described in Saba ' \
                          'et al 2019 https://doi.org/10.3389/fmars.2019.00664'.format(np.round(ta_equ['m'], 2),
                                                                                       np.round(ta_equ['b'], 2))

    return ta


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


def find_configs(deployment):
    # configdir = '/Users/garzio/Documents/repo/lgarzio/phglider/config'
    configdir = '/home/lgarzio/repo/lgarzio/phglider/config'  # in server
    global_attributes = '{}/{}/global_attributes.json'.format(configdir, deployment)
    variable_attrs = '{}/variable_attrs.json'.format(configdir)
    instruments = '{}/{}/instruments.json'.format(configdir, deployment)

    return global_attributes, variable_attrs, instruments


def get_erddap_dataset(server, ds_id, variables=None, constraints=None):
    variables = variables or None
    constraints = constraints or None

    e = ERDDAP(server=server,
               protocol='tabledap',
               response='nc')
    e.dataset_id = ds_id
    if constraints:
        e.constraints = constraints
    if variables:
        e.variables = variables
    ds = e.to_xarray()
    ds = ds.sortby(ds.time)
    return ds


def glider_region(ds):
    try:
        extent = [np.nanmin(ds.longitude.values) - 1.75, np.nanmax(ds.longitude.values) + 1.75,
                  np.nanmin(ds.latitude.values) - 1.5, np.nanmax(ds.latitude.values) + 1.5]
    except AttributeError:
        extent = [np.nanmin(ds.Longitude.values) - 1.75, np.nanmax(ds.Longitude.values) + 1.75,
                  np.nanmin(ds.Latitude.values) - 1.5, np.nanmax(ds.Latitude.values) + 1.5]
    region = dict()
    region['extent'] = extent

    return region


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


def plot_vars():
    plt_vars = {'conductivity': {'cmap': 'jet', 'ttl': 'Conductivity (S m-1)'},
                'temperature': {'cmap': cmo.cm.thermal, 'ttl': 'Temperature ({})'.format(r'$\rm ^oC$')},
                'salinity': {'cmap': cmo.cm.haline, 'ttl': 'Salinity'},
                'density': {'cmap': cmo.cm.dense, 'ttl': 'Density (kg m-3)'},
                'chlorophyll_a': {'cmap': cmo.cm.algae, 'ttl': 'Chlorophyll ({}g/L)'.format(chr(956))},
                'oxygen_concentration_mgL': {'cmap': cmo.cm.oxy, 'ttl': 'Oxygen (mg/L)'},
                'oxygen_concentration': {'cmap': cmo.cm.oxy, 'ttl': 'Oxygen (umol/L)'},
                'oxygen_concentration_shifted': {'cmap': cmo.cm.oxy, 'ttl': 'Oxygen (shifted) (umol/L)'},
                'sbe41n_ph_ref_voltage': {'cmap': cmo.cm.matter, 'ttl': 'pH Reference Voltage'},
                'sbe41n_ph_ref_voltage_shifted': {'cmap': cmo.cm.matter, 'ttl': 'pH Reference Voltage (shifted)'},
                'ph_total': {'cmap': cmo.cm.matter, 'ttl': 'pH'},
                'ph_total_shifted': {'cmap': cmo.cm.matter, 'ttl': 'pH (shifted)'}
                }
    # plt_vars = {'ph_total': {'cmap': cmo.cm.matter, 'ttl': 'pH'},
    #             'ph_total_shifted': {'cmap': cmo.cm.matter, 'ttl': 'pH (shifted)'}
    #             }
    return plt_vars


def plot_vars_dac():
    plt_vars = {'Conductivity': {'cmap': 'jet', 'ttl': 'Conductivity (S m-1)'},
                'Temperature': {'cmap': cmo.cm.thermal, 'ttl': 'Temperature ({})'.format(r'$\rm ^oC$')},
                'Salinity': {'cmap': cmo.cm.haline, 'ttl': 'Salinity'},
                'Chlorophyll': {'cmap': cmo.cm.algae, 'ttl': 'Chlorophyll ({}g/L)'.format(chr(956))},
                'Oxygen_mgL': {'cmap': cmo.cm.oxy, 'ttl': 'Oxygen (mg/L)'},
                'Oxygen_molar': {'cmap': cmo.cm.oxy, 'ttl': 'Oxygen (umol/L)'},
                'pHReferenceVoltage': {'cmap': cmo.cm.matter, 'ttl': 'pH Reference Voltage'},
                'pH': {'cmap': cmo.cm.matter, 'ttl': 'pH'},
                'TotalAlkalinity': {'cmap': cmo.cm.matter, 'ttl': 'Total Alkalinity'},
                'AragoniteSaturationState': {'cmap': cmo.cm.matter, 'ttl': 'Aragonite Saturation State'}
                }
    return plt_vars


def run_co2sys_ta_ph(ta, ph, sal, temp=25, press_dbar=0):
    """
    Runs the PyCO2SYS function using input TA and pH data.
    opt_pH_scale=1 is the default (Total scale), including here for clarity.
    opt_k_carbonic=4 is MCHP73 (Mehrbach et al 1973) refit by DM87 (Dickson & Millero 1987)
    opt_k_bisulfate=1 is the default (Dickson 1990), including here for clarity
    opt_total_borate=1 is the default (Uppstrom 1974), including here for clarity
    opt_k_fluoride=2 is PF87 (Perez & Fraga 1987)
    :param ta: Total Alkalinity array
    :param ph: pH array
    :param sal: salinity array
    :param temp: temperature array, default=25
    :param press_dbar: pressure (in units=dbar) array, default=0
    :return: calculated aragonite saturation, pCO2, revelle factor
    """
    # define input conditions
    par1 = ta  # Total Alkalinity
    par1_type = 1  # parameter 1 type (TA)
    par2 = ph
    par2_type = 3  # parameter 2 type (pH)

    kwargs = dict(salinity=sal,
                  temperature=temp,
                  pressure=press_dbar,
                  opt_pH_scale=1,
                  opt_k_carbonic=4,
                  opt_k_bisulfate=1,
                  opt_total_borate=1,
                  opt_k_fluoride=2)

    results = pyco2.sys(par1, par2, par1_type, par2_type, **kwargs)
    omega_arag = results['saturation_aragonite']  # aragonite saturation state
    pco2 = results['pCO2']  # units = uatm
    revelle = results['revelle_factor']

    return omega_arag, pco2, revelle


def yos(ds):
    """
    Identify the down-up profile pairs (yo), profile times and indices
    """
    profiletimes, idx = np.unique(ds.profile_time.values, return_index=True)
    depth = ds.depth.values

    # calculate segments
    idxs = []
    yos = []
    indices = []
    directions = []
    for i, ii in enumerate(idx):
        if i > 0:
            if ii < i_f:  # it's already been included in a previous profile
                continue
        if i < len(idx) - 1:
            i_f = idx[i + 1]
            # check the next profile time index. if the next profile contains <4 data points, include it in this profile
            try:
                i_f2 = idx[i + 2]
                if i_f2 - i_f < 4:
                    i_f = i_f2
            except IndexError:
                print('')
        else:
            i_f = len(depth)
        while i_f < ii:
            # skip if the profile indices are out of order, it looks like there are just a few profile times inserted
            # in the wrong place for some rt datasets (e.g. ru30-20210503T1929)
            print("Profile time indices aren't in ascending order {} {}".format(ii, ds.time.values[ii]))
            i_f = idx[i + 2]
            i += 1
        y = depth[ii:i_f]
        y = y[~np.isnan(y)]
        if len(y) > 20:
            y1 = np.nanmean(y[0:10])
            y2 = np.nanmean(y[-10:])
        else:
            y1 = np.nanmean(y[0:3])
            y2 = np.nanmean(y[-3:])
        if y1 < y2:
            dd = 'down'
        else:
            dd = 'up'
        if len(directions) > 0:
            if dd not in directions:
                directions.append(dd)
                indices.append(ii)
                yos.append(directions)
                idxs.append(indices)
                indices = []
                directions = []
            else:
                yos.append(directions)
                idxs.append(indices)
                indices = [ii]
                directions = [dd]
        else:
            directions.append(dd)
            indices.append(ii)
            if dd == 'up':
                yos.append(directions)
                idxs.append(indices)
                indices = []
                directions = []

    for i, idx in enumerate(idxs):
        if i < len(idxs) - 1:
            idx.append(idxs[i + 1][0])
        else:
            idx.append(len(ds.time))

    return profiletimes, yos, idxs
