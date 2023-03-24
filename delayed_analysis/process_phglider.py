#!/usr/bin/env python

"""
Author: Lori Garzio on 3/28/2022
Last modified: 3/23/2023
Process delayed-mode pH glider data.
1. Apply QARTOD QC flags (downloaded in the files) to CTD and DO data (set data flagged as 3/SUSPECT and 4/FAIL to nan).
2. Set profiles flagged as 3/SUSPECT and 4/FAIL from CTD hysteresis tests to nan (conductivity, temperature,
salinity and density).
3. Run location QARTOD test.
4. Convert pH voltages of 0.0 to nan.
5. Apply CTD thermal lag correction.
6. Interpolate (method=linear) pressure.
7. Interpolate (method=linear) temperature and salinity.
8. Remove interpolated data for profiles that failed hysteresis tests.
9. Calculate pH from original and corrected voltages and interpolated CTD data.
10. Calculate Total Alkalinity from interpolated salinity using a linear relationship determined from in-situ water
sampling.
11. Run ioos_qc gross range test on additional variables defined in the gross_range.yml config file (e.g. calculated pH,
chlorophyll-a) and apply test results to data.
12. Run ioos_qc spike test on pH and corrected pH, and apply test results to data.
13. TODO Run QARTOD rate of change test on pH and corrected pH (not implemented)
14. Calculate CO2SYS variables using TA, corrected pH, interpolated salinity, interpolated temperature, interpolated
pressure.
15. Convert oxygen concentration to mg/L.
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import json
from ioos_qc import qartod
from ioos_qc.config import Config
from ioos_qc.streams import XarrayStream
from ioos_qc.results import collect_results
from ioos_qc.utils import load_config_as_dict as loadconfig
import functions.common as cf
import functions.phcalc as phcalc
import functions.qc as fqc
import thermal_lag.ctd_thermal_lag_correction as ctd_thermal_lag
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
        'ancillary_variables': 'pressure_interpolated temperature_interpolated salinity_interpolated '
                               'ph_total_shifted total_alkalinity',
        'observation_type': 'calculated',
        'units': units,
        'long_name': long_name,
        'comment': 'Calculated using the PyCO2SYS function with inputs of interpolated pressure, '
                   'interpolated temperature, interpolated salinity, Total Alkalinity and corrected pH'
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


def combine_data(var1, var2):
    """
    Combine two datasets. var2 is added only where var1=nan
    :param var1: numpy array
    :param var2: numpy array
    """
    var_combined = var1
    nan_idx = np.isnan(var_combined)
    var_combined[nan_idx] = var2[nan_idx]

    return var_combined


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


def main(coord_lims, configdir, fname, method):
    deploy = '-'.join(fname.split('/')[-1].split('-')[0:2])
    ds = xr.open_dataset(fname)
    ds = ds.drop_vars(names=['profile_id', 'rowSize'])
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)
    savefile = f'{fname.split(".nc")[0]}_qc.nc'

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
    ds = fqc.location_test(ds, coord_lims)
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

    # read file containing attributes for additional variables added to the dataset (pH and thermal lag variables)
    va = '/Users/garzio/Documents/repo/lgarzio/phglider/config/ph_thermal_lag_attrs.json'
    with open(va) as json_file:
        var_attrs = json.load(json_file)

    # convert to dataframe and drop duplicated timestamps
    df = ds.to_dataframe()

    # drop duplicated timestamps
    df = df.loc[~df.index.duplicated(keep='first')]

    # apply CTD thermal lag adjustment
    df = ctd_thermal_lag.apply_thermal_lag(df)

    # rename thermal lag columns
    df.rename(columns={'salt_outside': 'salinity_lag_shifted',
                       'ctemp_outside': 'temperature_lag_shifted',
                       'rho_outside': 'density_lag_shifted'}, inplace=True)

    # generate combined data arrays - use thermal lag corrected data when available, otherwise use raw data
    df['temperature_combined'] = combine_data(np.array(df['temperature_lag_shifted']), np.array(df['temperature']))
    df['salinity_combined'] = combine_data(np.array(df['salinity_lag_shifted']), np.array(df['salinity']))
    df['density_combined'] = combine_data(np.array(df['density_lag_shifted']), np.array(df['density']))

    # interpolate pressure temperature and salinity (needed to calculate pH)
    df['pressure_interpolated'] = df['pressure'].interpolate(method='linear', limit_direction='both', limit=2)
    df['temperature_interpolated'] = df['temperature_combined'].interpolate(method='linear', limit_direction='both', limit=2)
    df['salinity_interpolated'] = df['salinity_combined'].interpolate(method='linear', limit_direction='both', limit=2)

    # if entire profiles were removed due to the CTD hysteresis tests, make sure the interpolated values are nan
    hysteresis_tests = ['conductivity_hysteresis_test', 'temperature_hysteresis_test']
    interp_cols = ['temperature_interpolated', 'salinity_interpolated']
    for ht in hysteresis_tests:
        ht_fail = np.logical_or(df[ht] == 3, df[ht] == 4)
        for col in interp_cols:
            df.loc[ht_fail, col] = np.nan

    # remove data where the pressure QARTOD QC flag value = 4/FAIL
    df = df.loc[df.pressure_qartod_summary_flag != 4].copy()

    # calculate pH and add to dataframe
    df['f_p'] = np.polyval([cc['f6'], cc['f5'], cc['f4'], cc['f3'], cc['f2'], cc['f1'], 0], df.pressure_interpolated)
    phfree, phtot = phcalc.phcalc(df.sbe41n_ph_ref_voltage, df.pressure_interpolated, df.temperature_interpolated,
                                  df.salinity_interpolated, cc['k0'], cc['k2'], df.f_p)
    df['ph_total'] = phtot

    phfree, phtot = phcalc.phcalc(df.sbe41n_ph_ref_voltage_shifted, df.pressure_interpolated,
                                  df.temperature_interpolated, df.salinity_interpolated, cc['k0'], cc['k2'], df.f_p)
    df['ph_total_shifted'] = phtot

    # calculate pH shifted without the thermal lag applied
    df['temperature_interpolated_notl'] = df['temperature'].interpolate(method='linear', limit_direction='both', limit=2)
    df['salinity_interpolated_notl'] = df['salinity'].interpolate(method='linear', limit_direction='both', limit=2)

    phfree, phtot = phcalc.phcalc(df.sbe41n_ph_ref_voltage_shifted, df.pressure_interpolated,
                                  df.temperature_interpolated_notl, df.salinity_interpolated_notl, cc['k0'], cc['k2'], df.f_p)
    df['ph_total_shifted_notl'] = phtot

    # there's a lot of noise in pH at the surface, so set pH values to nan when pressure < 1 dbar
    df.loc[df['pressure_interpolated'] < 1, 'ph_total'] = np.nan
    df.loc[df['pressure_interpolated'] < 1, 'ph_total_shifted'] = np.nan
    df.loc[df['pressure_interpolated'] < 1, 'ph_total_shifted_notl'] = np.nan

    df.set_index('time', inplace=True)

    # convert the dataframe to an xarray dataset
    phds = df.to_xarray()

    # assign attributes from the original dataset
    original_time = delete_attr(ds['time'])
    phds['time'].attrs = original_time.attrs
    phds['time'].encoding = original_time.encoding
    for varname in list(phds.data_vars):
        try:
            # if the variable is from the original dataset, use those attributes and encoding
            original_da = delete_attr(ds[varname])
            phds[varname].attrs = original_da.attrs
            phds[varname].encoding = original_da.encoding
        except KeyError:
            # for new variables that were added to this dataset, assign new attributes (grabbing from the attributes
            # from the original dataset for the new CTD variables
            attr_mapping = var_attrs[varname]
            varname0 = varname.split('_')[0]
            if varname0 in ['conductivity', 'temperature', 'salinity', 'density', 'pressure']:
                # for CTD variables, get the attributes from the original variable to add here
                original_da = delete_attr(ds[varname0])
                phds[varname].attrs = original_da.attrs
                for k, v in attr_mapping.items():
                    phds[varname].attrs[k] = v
                phds[varname].encoding = original_da.encoding
            else:
                # for the pH variables, just use the attributes from the config file
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

    # Run ioos_qc gross range test based on the configuration file
    c = Config(os.path.join(configdir, 'gross_range.yml'))
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
        attrs = fqc.assign_qartod_attrs(test, sensor, c.config[sensor]['qartod'][test])
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

    # run ioos_qc spike test on pH - original and time-shifted
    spike_settings = loadconfig(os.path.join(configdir, 'spike_test.yml'))['spike_settings']

    varnames = ['ph_total', 'ph_total_shifted']
    for vn in varnames:
        data = phds[vn]
        non_nan_ind = np.invert(np.isnan(data))  # identify where not nan
        non_nan_i = np.where(non_nan_ind)[0]  # get locations of non-nans
        tdiff = np.diff(data.time[non_nan_ind]).astype('timedelta64[s]').astype(float)  # get time interval (s) between non-nan points
        tdiff_long = np.where(tdiff > 60 * 5)[0]  # locate time intervals > 5 min
        tdiff_long_i = np.append(non_nan_i[tdiff_long], non_nan_i[tdiff_long + 1])  # original locations of where time interval is long

        # convert original threshold from units/s to units/average-timestep
        spike_settings['suspect_threshold'] = spike_settings['suspect_threshold'] * np.nanmedian(tdiff)
        spike_settings['fail_threshold'] = spike_settings['fail_threshold'] * np.nanmedian(tdiff)

        flag_vals = 2 * np.ones(np.shape(data))
        flag_vals[np.invert(non_nan_ind)] = qartod.QartodFlags.MISSING

        # only run the test if the array has values
        if len(non_nan_i) > 0:
            flag_vals[non_nan_ind] = qartod.spike_test(inp=data.values[non_nan_ind],
                                                       **spike_settings)

            # flag as not evaluated/unknown on either end of long time gap
            flag_vals[tdiff_long_i] = qartod.QartodFlags.UNKNOWN

            qc_varname = f'{vn}_qartod_spike_test'
            attrs = fqc.assign_qartod_attrs('spike_test', vn, spike_settings)
            da = xr.DataArray(flag_vals.astype('int32'), coords=data.coords, dims=data.dims, attrs=attrs,
                              name=qc_varname)
            phds[qc_varname] = da

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

    if not method == 'phonly':
        # calculate Total Alkalinity and add to dataset
        # TA calculated from salinity using a linear relationship determined from in-situ water sampling data taken during
        # glider deployment and recovery. See ../ta_equation/ta_sal_regression.py
        ta = cf.calculate_ta(deploy, phds.salinity_interpolated)
        phds['total_alkalinity'] = ta

        # run CO2SYS
        omega_arag, pco2, revelle = cf.run_co2sys_ta_ph(ta.values,
                                                        phds.ph_total_shifted.values,
                                                        phds.salinity_interpolated.values,
                                                        phds.temperature_interpolated.values,
                                                        phds.pressure_interpolated.values)

        # add the CO2SYS variables to the dataset
        attrs = assign_co2sys_attrs(omega_arag, '1', 'Aragonite Saturation State')
        da = xr.DataArray(omega_arag.astype('float32'), coords=ta.coords, dims=ta.dims, attrs=attrs,
                          name='aragonite_saturation_state')
        phds['aragonite_saturation_state'] = da

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
        da = xr.DataArray(data.astype('float32'), coords=phds[dov].coords, dims=phds[dov].dims, attrs=attrs, name=name)
        phds[name] = da

    phds.to_netcdf(savefile)


if __name__ == '__main__':
    coordinate_lims = {'latitude': [36, 43], 'longitude': [-76, -66]}
    configs = '/Users/garzio/Documents/repo/lgarzio/phglider/config'
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed.nc'
    method = 'final'  # 'final' 'phonly'
    main(coordinate_lims, configs, ncfile, method)
