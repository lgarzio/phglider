#!/usr/bin/env python

"""
Author: Lori Garzio on 3/28/2022
Last modified: 4/5/2022
Process delayed-mode pH glider data.
1. Apply QARTOD QC to CTD data (set data flagged as 3/SUSPECT and 4/FAIL to nan).
2. Set profiles flagged as 3/SUSPECT and 4/FAIL from CTD hysteresis tests to nan (conductivity, temperature,
salinity and density).
3. Run location QARTOD test.
4. Convert pH voltages of 0.0 to nan.
5. Interpolate (method=linear) pressure and remove data where pressure <1 dbar.
6. Interpolate (method=linear) temperature and salinity.
7. Remove interpolated data for profiles that failed hysteresis tests.
8. Calculate pH from original and corrected voltages and interpolated CTD data.
9. Calculate Total Alkalinity from interpolated salinity using a linear relationship determined from in-situ water
sampling.
10. Run ioos_qc gross range test on additional variables defined in the gross_range.yml config file (e.g. calculated pH,
chlorophyll-a) and apply test results to data.
11. Run QARTOD spike test on pH and corrected pH, and apply test results to data.
12. TODO Run QARTOD rate of change test on pH and corrected pH (not implemented)
13. Calculate CO2SYS variables using TA, corrected pH, interpolated salinity, interpolated temperature, interpolated
pressure.
14. Convert oxygen concentration to mg/L.
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import json
from ioos_qc.config import Config
from ioos_qc.streams import XarrayStream
from ioos_qc.results import collect_results
import functions.common as cf
import functions.phcalc as phcalc
import functions.qc as qc
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def assign_co2sys_attrs(array, units, long_name):
    """
    Define the CO2SYS variable attributes
    :param array: numpy array containing data values
    :param units: variable units
    :param long_name: variable long name
    """
    vmin = np.round(np.nanmin(array), 2)
    vmax = np.round(np.nanmax(array), 2)

    # Defining gross/flatline QC variable attributes
    attrs = {
        'actual_range': [vmin, vmax],
        'observation_type': 'calculated',
        'units': units,
        'long_name': long_name,
        'comment': 'Calculated using the PyCO2SYS function with inputs of interpolated pressure, '
                   'interpolated temperature, interpolated salinity, Total Alkalinity and pH'
    }

    return attrs


def assign_dissolved_oxygen_attrs(array, original_varname):
    """
    Define the transformed dissolved oxygen variable attributes
    :param array: numpy array containing data values
    :param original_varname: original DO variable name
    """
    vmin = np.round(np.nanmin(array), 2)
    vmax = np.round(np.nanmax(array), 2)

    # Defining gross/flatline QC variable attributes
    attrs = {
        'actual_range': [vmin, vmax],
        'ancillary_variables': f'sci_oxy4_oxygen {original_varname}',
        'observation_type': 'calculated',
        'units': 'mg L-1',
        'long_name': 'Dissolved Oxygen',
        'source_sensor': 'sci_oxy4_oxygen',
        'standard_name': 'mass_concentration_of_chlorophyll_a_in_seawater',
        'comment': f'transformed from umol/L to mg/L from {original_varname}'
    }

    return attrs


def delete_attr(da):
    """
    Delete these local attributes because they are no longer valid with QC'd dataset
    :param da: DataArray of variable
    :return: DataArray of variable with local attributes removed
    """

    for item in ['actual_range']:
        try:
            del da.attrs[item]
        except KeyError:
            continue
    return da


def main(coord_lims, grconfig, fname):
    deploy = '-'.join(fname.split('/')[-1].split('-')[0:2])
    ds = xr.open_dataset(fname)
    ds = ds.drop_vars(names=['profile_id', 'rowSize'])
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)

    # apply QARTOD QC to all variables except pressure
    qcvars = [x for x in list(ds.data_vars) if '_qartod_summary_flag' in x]
    for qv in qcvars:
        if 'pressure' in qv:
            continue
        target_var = list([qv.split('_qartod_summary_flag')[0]])
        if target_var[0] in ['conductivity', 'temperature']:
            target_var.append('salinity')
            target_var.append('density')
        qc_idx = np.where(np.logical_or(ds[qv].values == 3, ds[qv].values == 4))[0]
        if len(qc_idx) > 0:
            for tv in target_var:
                ds[tv][qc_idx] = np.nan

    # apply CTD hysteresis test QC
    qcvars = [x for x in list(ds.data_vars) if '_hysteresis_test' in x]
    for qv in qcvars:
        target_var = list([qv.split('_hysteresis_test')[0]])
        target_var.append('salinity')
        target_var.append('density')
        qc_idx = np.where(np.logical_or(ds[qv].values == 3, ds[qv].values == 4))[0]
        if len(qc_idx) > 0:
            for tv in target_var:
                ds[tv][qc_idx] = np.nan

    # QARTOD Location Test
    ds = qc.location_test(ds, coord_lims)
    qv = 'location_qartod_test'
    qc_idx = np.where(np.logical_or(ds[qv].values == 3, ds[qv].values == 4))[0]
    if len(qc_idx) > 0:
        for varname in list(ds.data_vars):
            ds[varname][qc_idx] = np.nan

    # convert pH voltages of 0.0 to nan
    phvolts = ['sbe41n_ph_ref_voltage', 'sbe41n_ph_ref_voltage_shifted']
    for phv in phvolts:
        ds[phv][ds[phv] == 0.0] = np.nan

    # read cal file
    sn = f'sbe{ds.instrument_ph.serial_number}'
    calfile = cf.find_calfile(deploy, sn)
    with open(calfile) as json_file:
        cc = json.load(json_file)

    # convert to dataframe and drop duplicated timestamps
    df = ds.to_dataframe()

    # interpolate pressure, and drop all data where pressure < 1 dbar
    df['pressure_interpolated'] = df['pressure'].interpolate(method='linear', limit_direction='both')
    df = df[df['pressure_interpolated'] >= 1].copy()

    # drop duplicated timestamps
    df = df[~df.index.duplicated(keep='first')]

    # interpolate pressure again, this time limiting the amount of interpolation
    df['pressure_interpolated'] = df['pressure'].interpolate(method='linear', limit_direction='both', limit=2)

    # interpolate temperature and salinity (needed to calculate pH)
    df['temperature_interpolated'] = df['temperature'].interpolate(method='linear', limit_direction='both', limit=2)
    df['salinity_interpolated'] = df['salinity'].interpolate(method='linear', limit_direction='both', limit=2)

    # if entire profiles were removed due to the CTD hysteresis tests, make sure the interpolated values are nan
    hysteresis_tests = ['conductivity_hysteresis_test', 'temperature_hysteresis_test']
    interp_cols = ['temperature_interpolated', 'salinity_interpolated']
    for ht in hysteresis_tests:
        ht_fail = np.logical_or(df[ht] == 3, df[ht] == 4)
        for col in interp_cols:
            df.loc[ht_fail, col] = np.nan

    # remove data where the pressure QARTOD QC flag value = 4/FAIL
    df = df[df.pressure_qartod_summary_flag != 4]

    # calculate pH and add to dataframe
    df['f_p'] = np.polyval([cc['f6'], cc['f5'], cc['f4'], cc['f3'], cc['f2'], cc['f1'], 0], df.pressure_interpolated)
    phfree, phtot = phcalc.phcalc(df.sbe41n_ph_ref_voltage, df.pressure_interpolated, df.temperature_interpolated,
                                  df.salinity_interpolated, cc['k0'], cc['k2'], df.f_p)
    df['ph_total'] = phtot

    phfree, phtot = phcalc.phcalc(df.sbe41n_ph_ref_voltage_shifted, df.pressure_interpolated,
                                  df.temperature_interpolated, df.salinity_interpolated, cc['k0'], cc['k2'], df.f_p)
    df['ph_total_shifted'] = phtot

    # convert the dataframe to an xarray dataset
    phds = df.to_xarray()

    # assign attributes from the original dataset
    original_time = delete_attr(ds['time'])
    phds['time'].attrs = original_time.attrs
    phds['time'].encoding = original_time.encoding
    for varname in list(phds.data_vars):
        try:
            original_da = delete_attr(ds[varname])
            phds[varname].attrs = original_da.attrs
            phds[varname].encoding = original_da.encoding
        except KeyError:
            try:
                attr_mapping = cf.attribute_mapping(varname)
                original_da = delete_attr(ds[attr_mapping['ancillary_variables']])
                phds[varname].attrs = original_da.attrs
                for k, v in attr_mapping.items():
                    phds[varname].attrs[k] = v
                phds[varname].encoding = original_da.encoding
            except KeyError:
                attr_mapping = cf.assign_ph_attrs(varname)
                phds[varname].attrs = attr_mapping

    # add calibration coefficients to dataset
    name = 'polynomial_coefficients'
    attributes = dict(
        comment='sensor calibration coefficients',
        ancillary_variables='instrument_ph',
        calibration_coefficients=str(cc)
    )
    da = xr.DataArray(phds.instrument_ph.values, coords=phds.instrument_ph.coords,
                      dims=phds.instrument_ph.dims, name=name, attrs=attributes)
    phds[name] = da

    phds = phds.assign_coords(
        {'latitude': phds.latitude, 'longitude': phds.longitude, 'depth': phds.depth})

    ds_global_attrs = ds.attrs
    phds = phds.assign_attrs(ds_global_attrs)

    # calculate Total Alkalinity and add to dataset
    # TA calculated from salinity using a linear relationship determined from in-situ water sampling data taken during
    # glider deployment and recovery. See ../ta_equation/ta_sal_regression.py
    ta = cf.calculate_ta(deploy, phds.salinity_interpolated)
    phds['total_alkalinity'] = ta

    # Run ioos_qc gross range test based on the configuration file
    c = Config(grconfig)
    xs = XarrayStream(phds, time='time', lat='latitude', lon='longitude')
    qc_results = xs.run(c)
    collected_list = collect_results(qc_results, how='list')

    # Parse each gross/flatline QC result
    for cl in collected_list:
        sensor = cl.stream_id
        test = cl.test
        qc_varname = f'{sensor}_{cl.package}_{test}'
        flag_results = cl.results.data

        # Define gross range QC variable attributes
        attrs = qc.assign_qartod_attrs(test, sensor, c.config[sensor]['qartod'][test])
        if not hasattr(phds[sensor], 'ancillary_variables'):
            phds[sensor].attrs['ancillary_variables'] = qc_varname
        else:
            phds[sensor].attrs['ancillary_variables'] = ' '.join((phds[sensor].ancillary_variables, qc_varname))

        # Create gross range data array and add to dataset
        da = xr.DataArray(flag_results.astype('int32'), coords=phds[sensor].coords, dims=phds[sensor].dims,
                          name=qc_varname,
                          attrs=attrs)
        phds[qc_varname] = da

        # Make a copy of the data and apply the QC flags
        original_da = phds[sensor].copy()
        if not hasattr(original_da, 'comment'):
            original_da.attrs['comment'] = 'No QC applied'
        else:
            original_da.attrs['comment'] = '. '.join((original_da.comment, 'No QC applied'))
        new_name = f'{original_da.name}_noqc'
        phds[new_name] = original_da

        qc_idx = np.where(np.logical_or(flag_results == 3, flag_results == 4))[0]
        if len(qc_idx) > 0:
            phds[sensor][qc_idx] = np.nan
        if not hasattr(phds[sensor], 'comment'):
            phds[sensor].attrs['comment'] = f'{len(qc_idx)} {test} flags applied'
        else:
            phds[sensor].attrs['comment'] = '. '.join((phds[sensor].comment, f'{len(qc_idx)} {test} flags applied'))

    # QARTOD Spike Test
    ph_thresholds = {'suspect': 0.05, 'fail': 0.2, 'seconds': 30}  # 0.05 is the detection limit of the sensor
    da = qc.spike_test(phds, 'ph_total', ph_thresholds)
    phds[da.name] = da

    da = qc.spike_test(phds, 'ph_total_shifted', ph_thresholds)
    phds[da.name] = da

    # apply spike test flags
    for qv in ['ph_total_qartod_spike_test', 'ph_total_shifted_qartod_spike_test']:
        target_var = qv.split('_qartod_spike_test')[0]
        qc_idx = np.where(np.logical_or(phds[qv].values == 3, phds[qv].values == 4))[0]
        if len(qc_idx) > 0:
            phds[target_var][qc_idx] = np.nan
        if not hasattr(phds[target_var], 'comment'):
            phds[target_var].attrs['comment'] = f'{len(qc_idx)} spike_test flags applied'
        else:
            phds[target_var].attrs['comment'] = '. '.join((phds[target_var].comment, f'{len(qc_idx)} spike_test flags applied'))

    # Rate of Change Test - TODO

    # run CO2SYS
    omega_arag, pco2, revelle = cf.run_co2sys_ta_ph(ta.values,
                                                    phds.ph_total_shifted.values,
                                                    phds.salinity_interpolated.values,
                                                    phds.temperature_interpolated.values,
                                                    phds.pressure_interpolated.values)

    # add the CO2SYS variables to the dataset
    attrs = assign_co2sys_attrs(omega_arag, '1', 'Aragonite Saturation State')
    da = xr.DataArray(omega_arag.astype('float32'), coords=ta.coords, dims=ta.dims, attrs=attrs,
                      name='saturation_aragonite')
    phds['saturation_aragonite'] = da

    attrs = assign_co2sys_attrs(pco2, 'uatm', 'pCO2')
    da = xr.DataArray(pco2.astype('float32'), coords=ta.coords, dims=ta.dims, attrs=attrs,
                      name='pco2_calculated')
    phds['pco2_calculated'] = da

    attrs = assign_co2sys_attrs(revelle, '1', 'Revelle Factor')
    da = xr.DataArray(revelle.astype('float32'), coords=ta.coords, dims=ta.dims, attrs=attrs, name='revelle_factor')
    phds['revelle_factor'] = da

    # calculate oxygen concentration in mg/L
    for dov in ['oxygen_concentration', 'oxygen_concentration_shifted']:
        name = f'{dov}_mgL'
        data = phds[dov].values * 32 / 1000  # convert from umol/L to mg/L
        attrs = assign_dissolved_oxygen_attrs(data, dov)
        da = xr.DataArray(data.astype('float32'), coords=ta.coords, dims=ta.dims, attrs=attrs, name=name)
        phds[name] = da

    phds.to_netcdf('{}_qc.nc'.format(fname.split('.nc')[0]))


if __name__ == '__main__':
    coordinate_lims = {'latitude': [36, 43], 'longitude': [-76, -66]}
    gross_range_config = '/Users/garzio/Documents/repo/lgarzio/phglider/config/gross_range.yml'
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed.nc'
    main(coordinate_lims, gross_range_config, ncfile)
