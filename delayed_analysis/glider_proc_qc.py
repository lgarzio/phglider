#!/usr/bin/env python

"""
Author: Lori Garzio on 4/5/2021
Last modified: 8/10/2021
Modified from MATLAB code written by Liza Wright-Fairbanks, Spring 2020.
Based on IOOS manual for real-time quality control of pH data observations.
Tests to be applied as a post-processing QA/QC measure.
Also calculate TA using salinity data.
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import functions.common as cf
import functions.qc as qc
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def time_test(da, sdir):
    """
    Test 1: IDs data gaps of pH data >1 hour
    :param da: xarray data array dimensioned by time
    :param sdir: full path to save file directory
    """
    # drop nan values
    da = da.dropna(dim='time')
    time_df = pd.DataFrame(da.time.values, columns=['time'])

    # calculate time difference between each timestamp and add to dictionary
    gap_list = []
    time_df['diff'] = time_df['time'].diff()
    idx_gap = time_df['diff'][time_df['diff'] > pd.Timedelta(hours=1)].index.tolist()
    if len(idx_gap) > 0:
        for i in idx_gap:
            t0 = time_df['time'][i - 1]
            tf = time_df['time'][i]
            gaplen_totalhours = ((tf - t0).days / 24) + np.round((tf - t0).seconds/3600, 2)
            gap_list.append([pd.to_datetime(t0).strftime('%Y-%m-%dT%H:%M:%S'),
                             pd.to_datetime(tf).strftime('%Y-%m-%dT%H:%M:%S'), gaplen_totalhours])
        ttest_df = pd.DataFrame(gap_list, columns=['gap_start', 'gap_end', 'gap_length_hours'])
        ttest_df.to_csv(os.path.join(sdir, '{}_timegaps.csv'.format(da.name)), index=False)

    print('\n{} data gaps >1 hour: {}'.format(da.name, len(gap_list)))


def main(coord_lims, phu_lims, fname):
    deploy = '-'.join(fname.split('/')[-1].split('-')[0:2])
    ds = xr.open_dataset(fname)
    savedir = os.path.join(os.path.dirname(fname), 'qc')
    os.makedirs(savedir, exist_ok=True)

    # calculate Total Alkalinity and add to dataset
    # TA calculated from salinity using a linear relationship determined from in-situ water sampling data taken during
    # glider deployment and recovery. See ../ta_equation/ta_sal_regression.py
    ta = cf.calculate_ta(deploy, ds.salinity)
    ds['total_alkalinity'] = ta

    # run CO2SYS
    if np.sum(~np.isnan(ds.ph_total_shifted.values)) > 0:
        ph = ds.ph_total_shifted.values
    else:
        ph = ds.ph_total.values
    omega_arag, pco2, revelle = cf.run_co2sys_ta_ph(ta.values, ph, ds.salinity.values, ds.temperature.values,
                                                    ds.sci_water_pressure_dbar.values)

    # add the variables to the dataset
    da = xr.DataArray(omega_arag, coords=ta.coords, dims=ta.dims, name='saturation_aragonite')
    da.attrs['observation_type'] = 'calculated'
    da.attrs['units'] = '1'
    da.attrs['long_name'] = 'Aragonite Saturation State'
    da.attrs['comment'] = 'Calculated using the PyCO2SYS function with inputs of Total Alkalinity and pH'
    ds['saturation_aragonite'] = da

    da = xr.DataArray(pco2, coords=ta.coords, dims=ta.dims, name='pco2_calculated')
    da.attrs['observation_type'] = 'calculated'
    da.attrs['units'] = 'uatm'
    da.attrs['long_name'] = 'pCO2'
    da.attrs['comment'] = 'Calculated using the PyCO2SYS function with inputs of Total Alkalinity and pH'
    ds['pco2_calculated'] = da

    da = xr.DataArray(revelle, coords=ta.coords, dims=ta.dims, name='revelle_factor')
    da.attrs['observation_type'] = 'calculated'
    da.attrs['units'] = '1'
    da.attrs['long_name'] = 'Revelle Factor'
    da.attrs['comment'] = 'Calculated using the PyCO2SYS function with inputs of Total Alkalinity and pH'
    ds['revelle_factor'] = da

    # calculate oxygen concentration in mg/L
    da = ds.oxygen_concentration * 32 / 1000  # change oxygen from umol/L to mg/L
    da.name = 'oxygen_concentration_mgL'
    da.attrs['observation_type'] = 'calculated'
    da.attrs['units'] = 'mg L-1'
    da.attrs['long_name'] = 'Dissolved Oxygen'
    da.attrs['ancillary_variables'] = 'sci_oxy4_oxygen'
    da.attrs['source_sensor'] = 'sci_oxy4_oxygen'
    da.attrs['standard_name'] = 'mass_concentration_of_chlorophyll_a_in_seawater'
    da.attrs['comment'] = 'transformed from umol/L to mg/L'
    ds['oxygen_concentration_mgL'] = da

    da = ds.oxygen_concentration_shifted * 32 / 1000  # change oxygen from umol/L to mg/L
    da.name = 'oxygen_concentration_mgL_shifted'
    da.attrs['observation_type'] = 'calculated'
    da.attrs['units'] = 'mg L-1'
    da.attrs['long_name'] = 'Dissolved Oxygen'
    da.attrs['ancillary_variables'] = 'sci_oxy4_oxygen'
    da.attrs['source_sensor'] = 'sci_oxy4_oxygen'
    da.attrs['standard_name'] = 'mass_concentration_of_chlorophyll_a_in_seawater'
    da.attrs['comment'] = 'transformed from umol/L to mg/L. {}'.format(ds.oxygen_concentration_shifted.comment)
    ds['oxygen_concentration_mgL_shifted'] = da

    # Test 1: QARTOD Time Test
    time_test(ds['ph_total'], savedir)
    time_test(ds['sbe41n_ph_ref_voltage'], savedir)
    if np.sum(~np.isnan(ds.ph_total_shifted.values)) > 0:
        time_test(ds['ph_total_shifted'], savedir)
        time_test(ds['sbe41n_ph_ref_voltage_shifted'], savedir)
    time_test(ds['oxygen_concentration'], savedir)
    if np.sum(~np.isnan(ds.oxygen_concentration_shifted.values)) > 0:
        time_test(ds['oxygen_concentration_shifted'], savedir)
    time_test(ds.chlorophyll_a, savedir)
    time_test(ds.temperature, savedir)
    time_test(ds.salinity, savedir)

    # Test 3: QARTOD Location Test
    ds = qc.location_test(ds, coord_lims)

    # Test 4: QARTOD Gross Range Test
    ph_sensor_lims = {'min': 0, 'max': 14}
    ds = qc.gross_range(ds, 'ph_total', ph_sensor_lims, phu_lims)
    if np.sum(~np.isnan(ds.ph_total_shifted.values)) > 0:
        ds = qc.gross_range(ds, 'ph_total_shifted', ph_sensor_lims, phu_lims)

    ph_volt_lims = {'min': np.nan, 'max': 0}
    ds = qc.gross_range(ds, 'sbe41n_ph_ref_voltage', ph_volt_lims)
    if np.sum(~np.isnan(ds.sbe41n_ph_ref_voltage_shifted.values)) > 0:
        ds = qc.gross_range(ds, 'sbe41n_ph_ref_voltage_shifted', ph_volt_lims)

    chl_lims = {'min': 0, 'max': 50}
    ds = qc.gross_range(ds, 'chlorophyll_a', chl_lims)

    # Test 6: QARTOD Spike Test
    ph_thresholds = {'low': 0.05, 'high': 0.2}  # 0.05 is the detection limit of the sensor
    ds = qc.spike_test(ds, 'ph_total', ph_thresholds)
    if np.sum(~np.isnan(ds.ph_total_shifted.values)) > 0:
        ds = qc.spike_test(ds, 'ph_total_shifted', ph_thresholds)

    # Test 7: Rate of Change Test - TODO

    ds.to_netcdf('{}_qc.nc'.format(fname.split('.')[0]))


if __name__ == '__main__':
    coordinate_lims = {'latitude': {'min': 36, 'max': 43}, 'longitude': {'min': -76, 'max': -66}}
    ph_user_lims = {'min': 6.5, 'max': 9}
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed_shifted.nc'
    main(coordinate_lims, ph_user_lims, ncfile)
