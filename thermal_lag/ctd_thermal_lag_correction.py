#!/usr/bin/env python

"""
Authors: Daniel Wang and Jack Slater, Virginia Institute of Marine Science
Apply CTD thermal lag corrections
Modified by Lori Garzio 3/1/2023
Last modified by Lori Garzio on 3/30/2023
"""

import pandas as pd
import numpy as np
from datetime import datetime
import gsw
import thermal_lag.thermal_lag_functions as tl

pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console

# ignores divide by 0 errors in later step
np.seterr(divide='ignore', invalid='ignore')


def apply_thermal_lag(glider_sci):
    # Step 1: Prepare deployment data

    temperature_diff_threshold = 0.25

    glider_sci.reset_index(inplace=True)  # reset indices

    # convert source file names from char to string, in preparation for data extraction
    glider_sci.source_file = str(glider_sci.source_file)
    glider_sci.trajectory = str(glider_sci.trajectory)

    # drop all rows without ctd timestamp and duplicate ctd timestamps and rename to sci_data
    sci_data = glider_sci[glider_sci['ctd41cp_timestamp'].notnull()]
    sci_data = sci_data[sci_data['pressure'] > 1]
    sci_data = sci_data.drop_duplicates(subset=['ctd41cp_timestamp'])
    sci_data.reset_index(inplace=True)  # reset indices

    # rename variables
    sci_data['z'] = sci_data['depth'].mul(-1)  # create variable z, negative of depth
    sci_data['ctd_time'] = sci_data['ctd41cp_timestamp'].apply(datetime.timestamp)

    # correction for CTD sensor response time (lag from measurement to recording).
    # this is not used in practice, because pressure sensor lag is assumed to be 0.
    # 0 means no correction, assuming pressure and temperature are recorded correctly and instantly at the CTD time stamp
    P_sensor_lag = 0
    T_sensor_lag = 0

    tm = sci_data['ctd_time']

    sci_data['pressure_lag_shifted'] = tl.correct_sensor_lag(tm, sci_data['pressure'], P_sensor_lag)
    sci_data['z_lag_shifted'] = tl.correct_sensor_lag(tm, sci_data['z'], P_sensor_lag)
    sci_data['temperature_lag_shifted'] = tl.correct_sensor_lag(tm, sci_data['temperature'], T_sensor_lag)

    # correct conductivity for CTD sensor lag time
    Vol = 1.5  # based on diagram from Kim Martini of Sea-Bird.
    Q = 10  # flow rate in ml/s.
    TC_sensor_lag = Vol / Q

    sci_data['conductivity_lag_shifted'] = tl.correct_sensor_lag(tm, sci_data['conductivity'], TC_sensor_lag)

    # smooth data using a moving mean filter
    # this step helps correctly indentify thermocline depth; otherwise result, otherwise the result is too spiky
    # conductivity might not need to be smoothed?
    mov_window = 5
    sci_data['pressure_lag_shifted_smooth'] = sci_data['pressure_lag_shifted'].rolling(mov_window, min_periods=1,
                                                                                       center=True).mean()
    sci_data['z_lag_shifted_smooth'] = sci_data['z_lag_shifted'].rolling(mov_window, min_periods=1, center=True).mean()
    sci_data['temperature_lag_shifted_smooth'] = sci_data['temperature_lag_shifted'].rolling(mov_window, min_periods=1,
                                                                                             center=True).mean()
    sci_data['conductivity_lag_shifted_smooth'] = sci_data['conductivity_lag_shifted'].rolling(mov_window, min_periods=1,
                                                                                               center=True).mean()

    # Thermister reponse correction for temperature data
    # This step is neccesary. assuming tau_T = 0.53 sec.
    # according to Kim Martini's slides:
    # "Response time for this temperature sensor construction can range from 0.1-0.6 seconds.
    # This can vary depending on the pump and profiling speed of the platform."

    # smooth raw temperature data. 3 measurement points corresponds to about 6 seconds.
    # this avoids extremely large dT/dt, dT/dz
    dt = sci_data['ctd_time'].diff()
    dt = dt[1:, ]
    tlsm = sci_data['temperature_lag_shifted_smooth'].diff()
    dT_dt_smooth = np.divide(tlsm[1:, ], dt)
    tau_T = 0.53  # in seconds. nominal value is 0.5 second based on Johnson et al. 2007
    sci_data['temperature_response_corrected_smooth'] = sci_data['temperature_lag_shifted_smooth'].add(
        np.multiply(tau_T, dT_dt_smooth))

    # Step 2 Prepare Profile Data

    sci_data['profile_id'] = sci_data.groupby('profile_time').ngroup()  # identify profiles by id in profile_id

    n_profiles = sci_data['profile_id'].nunique()

    # find profile direction (up or down), find real profiles (pressure range > 5 dbar)
    profile_pressure_range_cutoff = 5  # dbar
    temperature_diff_cutoff = 4  # C
    profile_stats = pd.DataFrame()  # create dataframe for profile statistics

    profile_stats['profile_time'] = sci_data.groupby('profile_id')['profile_time'].agg('first')
    profile_stats['n_profile_values'] = sci_data.groupby('profile_id')['ctd_time'].agg('count')
    profile_stats['pressure_diff'] = sci_data.groupby('profile_id')['pressure'].agg('last') - \
                                     sci_data.groupby('profile_id')['pressure'].agg('first')
    profile_stats['temperature_diff'] = sci_data.groupby('profile_id')['temperature'].agg('max') - \
                                        sci_data.groupby('profile_id')['temperature'].agg('min')

    profile_stats['profile_direction'] = np.select([profile_stats['pressure_diff'] >= profile_pressure_range_cutoff, \
                                                    profile_stats['pressure_diff'] <= -profile_pressure_range_cutoff],
                                                   [1, -1], 0)  # 1 is downcast, -1 is upcast, 0 is null

    profile_stats['stratification_flag'] = np.select([profile_stats['temperature_diff'] >= temperature_diff_cutoff], [1], 0)

    profile_groups = sci_data.groupby('profile_id')

    def find_interface_thickness(group):
        """
        Calculates interface (thermocline) thickness for two layer water column (if stratification_flag == 1) .
        Interface thickness defined as the depth range for middle 70 % temperature values.
        range(min_T + 0.15 * T_diff, max_T - 0.15 * T_diff)
        Args:
            group (pandas groupby object): Profiles grouped by profile_id
        Returns:
            Profile stats dataframe with interface_thickness column
        """
        profile_id = group.iloc[0]['profile_id']

        if profile_stats.loc[profile_id, 'stratification_flag'] == 1:

            temp = group['temperature']
            pressure = group['pressure']
            tempdiff = profile_stats.loc[profile_id, 'temperature_diff']
            mintemp = temp.min() + 0.15 * tempdiff
            maxtemp = temp.max() - 0.15 * tempdiff
            indices = temp[(temp >= mintemp) & (temp <= maxtemp)].index
            interface_measurements_count = len(indices)
            profile_stats.loc[profile_id, 'interface_measurements_count'] = interface_measurements_count
            max_pressure = pressure[indices].max()
            min_pressure = pressure[indices].min()
            interface_thickness = np.abs(max_pressure - min_pressure)
            if np.isnan(interface_thickness):
                interface_thickness = 0.1

        else:
            interface_thickness = np.nan

        profile_stats.loc[profile_id, 'interface_thickness'] = interface_thickness

    profile_groups.apply(find_interface_thickness)
    profile_groups = sci_data.groupby("profile_id")

    def find_gradient_per_profile(group):
        """
        Finds first and second order numpy gradient of dT/dz for each profile (calculate dT/dz)
        Args:
            group (pandas groupby object): Profiles grouped by profile_id
        Returns:
            Profile dataframe with gradient columns
        """
        if len(group) > 1:
            gradient1 = np.gradient(group['temperature_response_corrected_smooth'], group['z_lag_shifted_smooth'])
            group['dT_dz_smooth'] = gradient1
            gradient2 = np.gradient(group['dT_dz_smooth'], group['z_lag_shifted_smooth'])
            group['d2T_dz2_smooth'] = gradient2
        else:
            group['dT_dz_smooth'] = np.nan
            group['d2T_dz2_smooth'] = np.nan
        return group

    sci_data = profile_groups.apply(find_gradient_per_profile)
    profile_groups = sci_data.groupby("profile_id")

    def find_thermocline_z_p(group):
        """
        Calculates thermocline depth and pressure using the point of maximum dT/dz
        Args:
            group (pandas groupby object): Profiles grouped by profile_id
        Returns:
            Profile stats dataframe with thermocline_z and thermocline_pressure columns
        """
        profile_id = group.iloc[0]['profile_id']
        ind_zrange = (group['z'] < -4) & (group['z'] > (2 + group['z'].min()))
        if group['z'][ind_zrange].empty:
            profile_stats.loc[profile_id, 'thermocline_z'] = np.nan
            profile_stats.loc[profile_id, 'thermocline_pressure'] = np.nan
        else:
            ind1 = group['dT_dz_smooth'].abs() == group['dT_dz_smooth'][ind_zrange].abs().max()
            profile_stats.loc[profile_id, 'thermocline_z'] = group['z_lag_shifted_smooth'][ind1].mean()
            profile_stats.loc[profile_id, 'thermocline_pressure'] = group['pressure_lag_shifted_smooth'][ind1].mean()

    profile_groups.apply(find_thermocline_z_p)
    profile_groups = sci_data.groupby("profile_id")

    # Identify which thermal lag method to apply for each profile
    # 0: no correction
    # 1: correction in T/S (or normalized T/S) space
    # 2: correction in in Pressure (depth) - Salinity space, adjusted to thermocline depth

    def assign_TS_flag(group):
        """
        Assign thermal lag flag as 1 (TS) or 0 (no correction) based on conditions
        Args:
            group (pandas groupby object): Profiles grouped by profile_id
        Returns:
            Profile stats dataframe with thermal_lag_flag column
        """
        profile_id = group.iloc[0]['profile_id']
        profile_stats.loc[profile_id, 'thermal_lag_flag'] = 0

        cond3 = group['pressure'].max() >= 10  # is the profile at least 10 dbar in the water column
        cond4 = np.abs(profile_stats.loc[profile_id, 'pressure_diff']) >= 10  # did the glider traverse at least 10 dbar
        # was the temperature difference throughout the water column past a threshold?
        cond5 = np.abs(profile_stats.loc[profile_id, 'temperature_diff']) >= temperature_diff_threshold

        if profile_id == 0:  # first profile
            # compare the first profile to the next profile to make sure one is down and one is up
            cond1 = profile_stats.loc[profile_id, 'profile_direction'] * profile_stats.loc[
                (profile_id + 1), 'profile_direction'] == -1
            profile_stats.loc[profile_id, 'thermal_lag_flag'] = np.select([cond1 & cond3 & cond4 & cond5], [1], 0)

        elif profile_id < (n_profiles - 1):
            cond1 = profile_stats.loc[profile_id, 'profile_direction'] * profile_stats.loc[
                (profile_id + 1), 'profile_direction'] == -1
            cond2 = profile_stats.loc[profile_id, 'profile_direction'] * profile_stats.loc[
                (profile_id - 1), 'profile_direction'] == -1
            profile_stats.loc[profile_id, 'thermal_lag_flag'] = np.select([(cond1 or cond2) & cond3 & cond4 & cond5], [1], 0)

        elif profile_id == (n_profiles - 1):  # last profile
            cond2 = profile_stats.loc[profile_id, 'profile_direction'] * profile_stats.loc[
                (profile_id - 1), 'profile_direction'] == -1
            profile_stats.loc[profile_id, 'thermal_lag_flag'] = np.select([cond2 & cond3 & cond4 & cond5], [1], 0)

    profile_groups.apply(assign_TS_flag)
    profile_groups = sci_data.groupby("profile_id")

    def assign_SP_flag(group):
        """
        Assign thermal lag flag as 2 (SP) or remain as current based on conditions
        Args:
            group (pandas groupby object): Profiles grouped by profile_id
        Returns:
            Profile stats dataframe with updated thermal_lag_flag column
        """
        profile_id = group.iloc[0]['profile_id']
        current = profile_stats.loc[profile_id, 'thermal_lag_flag']
        cond6 = profile_stats.loc[profile_id, 'interface_thickness'] < 8
        cond7 = ~np.isnan(profile_stats.loc[profile_id, 'interface_thickness'])
        cond8 = profile_stats.loc[profile_id, 'interface_measurements_count'] <= 32
        cond9 = group['pressure'].max() - profile_stats.loc[profile_id, 'thermocline_pressure'] >= 2
        cond10 = profile_stats.loc[profile_id, 'thermocline_pressure'] - group['pressure'].min() >= 2
        cond11 = profile_stats.loc[profile_id, 'thermal_lag_flag'] == 1

        profile_stats.loc[profile_id, 'thermal_lag_flag'] = np.select([cond6 & cond7 & cond8 & cond9 & cond10 & cond11],
                                                                      [2], current)

    profile_groups.apply(assign_SP_flag)
    profile_groups = sci_data.groupby("profile_id")

    # Step 3

    def run_thermal_lag_params(group):
        """
        Performs thermal lag correction using optimization function, then calculates final corrected profile values
        High-level processing steps:
            1. Identify which thermal lag correction to use (TS or SP)
            2. Identify whether it is a downcast or and upcast, and choose the corresponding upcast or downcast to pair
            with. The pairing profile should be either the previous profile or the next profile adjacent to the current
            profile.
            3. Perform pair-wise thermal lag correction for the current profile.
            4. Calculate derivative variables (corrected salinity, density, etc.) for the current profile
        Args:
            group (pandas groupby object): Profiles grouped by profile_id
        Returns:
            sci_data_cor dataframe with new columns of corrected values
        """

        def compare_profile_ranges(pressure1, pressure2, max_threshold=.25, diff_threshold=.25):
            # check maximum pressure
            press1_max = np.nanmax(pressure1)
            press2_max = np.nanmax(pressure2)
            min_maxpress = min(press1_max, press2_max)
            max_maxpress = max(press1_max, press2_max)
            criteria1 = (max_maxpress - (max_maxpress * max_threshold)) <= min_maxpress

            # check the pressure range
            press1_diff = abs(np.nanmin(pressure1) - press1_max)
            press2_diff = abs(np.nanmin(pressure2) - press2_max)
            min_diff = min(press1_diff, press2_diff)
            max_diff = max(press1_diff, press2_diff)
            criteria2 = (max_diff - (max_diff * diff_threshold)) <= min_diff

            return (criteria1 & criteria2)

        profile_id = group.iloc[0]['profile_id']
        if profile_stats.loc[profile_id, 'thermal_lag_flag'] != 0:
            try:
                # first profile
                if (profile_id == 0) & (profile_stats.loc[(profile_id + 1), 'thermal_lag_flag'] != 0):
                    pair_group = profile_groups.get_group(profile_id + 1)

                    # check to make sure the profiles don't span different depth ranges
                    pressure_range_test = compare_profile_ranges(group.pressure, pair_group.pressure)
                    if not pressure_range_test:
                        raise Exception("Profiles failed the comparison test: no valid profile to correct with")

                # last profile
                elif (profile_id == n_profiles - 1) & (profile_stats.loc[(profile_id - 1), 'thermal_lag_flag'] != 0):
                    pair_group = profile_groups.get_group(profile_id - 1)

                    # check to make sure the profiles don't span different depth ranges
                    pressure_range_test = compare_profile_ranges(group.pressure, pair_group.pressure)
                    if not pressure_range_test:
                        raise Exception("Profiles failed the comparison test: no valid profile to correct with")

                else:
                    belowid = profile_id - 1
                    aboveid = profile_id + 1

                    below = np.abs(
                        profile_stats.loc[profile_id, 'profile_time'] - profile_stats.loc[belowid, 'profile_time'])
                    above = np.abs(
                        profile_stats.loc[profile_id, 'profile_time'] - profile_stats.loc[aboveid, 'profile_time'])

                    group_below = profile_groups.get_group(belowid)
                    group_above = profile_groups.get_group(aboveid)

                    pressure_range_test_below = compare_profile_ranges(group.pressure, group_below.pressure)
                    pressure_range_test_above = compare_profile_ranges(group.pressure, group_above.pressure)

                    # find the closest profile in time that meets all of the criteria to correct for thermal lag and
                    # has a similar depth range as the profile being tested
                    if (below < above) & (profile_stats.loc[belowid, 'thermal_lag_flag'] != 0) & pressure_range_test_below:
                        pair_group = profile_groups.get_group(belowid)
                    elif (below > above) & (profile_stats.loc[aboveid, 'thermal_lag_flag'] != 0) & pressure_range_test_above:
                        pair_group = profile_groups.get_group(aboveid)

                    # if the closest profile doesn't meet all of the criteria to correct for thermal lag and is not a
                    # comparable depth range, check the other profile to 1. make sure it isn't more than twice as far
                    # away in time as the closer profile, 2. make sure it meets the criteria to correct for thermal lag,
                    # and 3. make sure it spans a comparable depth range
                    elif (below < above * 2) & (profile_stats.loc[belowid, 'thermal_lag_flag'] != 0) & pressure_range_test_below:
                        pair_group = profile_groups.get_group(belowid)
                    elif (below * 2 > above) & (profile_stats.loc[aboveid, 'thermal_lag_flag'] != 0) & pressure_range_test_above:
                        pair_group = profile_groups.get_group(aboveid)
                    else:
                        raise Exception("No valid profile to correct with")

                profile_id2 = pair_group.iloc[0]['profile_id']

                time1 = np.array(group['ctd_time'])
                temp1 = np.array(group['temperature'])
                cond1 = np.array(group['conductivity'])
                pres1 = np.array(group['pressure'])
                thermocline_pres1 = profile_stats.loc[profile_id, 'thermocline_pressure']
                time2 = np.array(pair_group['ctd_time'])
                temp2 = np.array(pair_group['temperature'])
                cond2 = np.array(pair_group['conductivity'])
                pres2 = np.array(pair_group['pressure'])
                thermocline_pres2 = profile_stats.loc[profile_id2, 'thermocline_pressure']
                lat1 = np.array(group['latitude'])
                lon1 = np.array(group['longitude'])

                # make sure none of the arrays are all nans (e.g. data removed in quality control)
                data_condition1 = (np.sum(~np.isnan(temp1)) > 1) & (np.sum(~np.isnan(cond1)) > 1) & (np.sum(~np.isnan(pres1)) > 1)
                data_condition2 = (np.sum(~np.isnan(temp2)) > 1) & (np.sum(~np.isnan(cond2)) > 1) & (np.sum(~np.isnan(pres2)) > 1)
                if data_condition1 & data_condition2:
                    # calculate the optimal alpha and tau for the profile pair by correcting for thermal lag and
                    # minimizing the area between the T-S curves
                    if profile_stats.loc[profile_id, 'thermal_lag_flag'] == 1:
                        # thermal lag correction in T/S (or normalized T/S) space
                        params = tl.find_thermal_lag_params_ts(time1, cond1, temp1, pres1, time2, cond2, temp2, pres2)

                    elif profile_stats.loc[profile_id, 'thermal_lag_flag'] == 2:
                        # thermal lag correction in Pressure (depth) - Salinity space, adjusted to thermocline depth
                        params = tl.find_thermal_lag_params_sp(time1, cond1, temp1, pres1, thermocline_pres1, time2, cond2,
                                                               temp2, pres2, thermocline_pres2)

                    else:
                        print('Invalid thermal lag flag')

                    [temp_inside1, cond_outside1] = tl.correct_thermal_lag(time1, cond1, temp1, params.x)
                    [temp_inside2, cond_outside2] = tl.correct_thermal_lag(time2, cond2, temp2, params.x)

                    # check if the corrected temperatures are an order of magnitude different than the original
                    temp_condition1a = np.nanmax(temp_inside1) > np.nanmax(temp1) * 10
                    temp_condition1b = abs(np.nanmin(temp_inside1)) > abs(np.nanmin(temp1)) * 10
                    temp_condition2 = (np.nanmax(temp_inside2) > np.nanmax(temp2) * 10)
                    temp_condition2b = abs(np.nanmin(temp_inside2)) > abs(np.nanmin(temp2)) * 10

                    # if any of these conditions is met, do no apply correction
                    if temp_condition1a | temp_condition1b | temp_condition2 | temp_condition2b:
                        raise Exception("Correction is invalid, not applying correction")

                    salt_cor1 = gsw.SP_from_C(np.multiply(cond_outside1, 10), temp1, pres1)  # corrected Practical Salinity

                    saltA_outside1 = gsw.SA_from_SP(salt_cor1, pres1, lon1, lat1)  # corrected Absolute Salinity

                    ctemp_outside1 = gsw.CT_from_t(saltA_outside1, temp1, pres1)  # corrected conservative temperature

                    ptemp_outside1 = gsw.pt_from_CT(saltA_outside1, ctemp_outside1)  # corrected potential temperature

                    rho_outside1 = gsw.rho(saltA_outside1, ctemp_outside1, pres1)  # in-situ density

                    sigma0_outside1 = gsw.sigma0(saltA_outside1, ctemp_outside1)  # potential density anomaly

                    profile_stats.loc[profile_id, 'alpha'] = params.x[0]
                    profile_stats.loc[profile_id, 'tau'] = params.x[1]

                    group['salt_outside'] = salt_cor1  # corrected Practical Salinity
                    group['saltA_outside'] = saltA_outside1  # corrected Absolute Salinity
                    group['ctemp_outside'] = ctemp_outside1  # corrected conservative temperature
                    group['ptemp_outside'] = ptemp_outside1  # corrected potential temperature
                    group['rho_outside'] = rho_outside1  # in-situ density
                    group['sigma0_outside'] = sigma0_outside1  # potential density anomaly
                else:
                    print('Not all required data available to calculate thermal lag')
            except Exception as e:
                print(profile_id)
                print(e)
        return group

    sci_data_cor = profile_groups.apply(run_thermal_lag_params)
    sci_data_cor.reset_index(drop=True, inplace=True)  # reset indices
    keep_cols = ['time', 'conductivity_lag_shifted', 'salt_outside', 'ctemp_outside', 'rho_outside']

    # merge the original glider_sci dataframe with the corrected CTD dataframe
    sci_data_cor = sci_data_cor[keep_cols]
    merged = glider_sci.merge(sci_data_cor, how='outer', on='time')

    return merged
