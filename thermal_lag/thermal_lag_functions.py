#!/usr/bin/env python

"""
Authors: Daniel Wang and Jack Slater, Virginia Institute of Marine Science
CTD thermal lag correction functions
"""

import numpy as np
from scipy.interpolate import interp1d
import gsw
import shapely.geometry as spg
import scipy.optimize as scop


def correct_sensor_lag(timestamp, raw, params, flow=0):
    """
    Authors: Daniel Wang and Jack Slater, Virginia Institute of Marine Science
    Python translation of SOCIB Glider Toolbox correctSensorLag function
    """
    if flow == 0:
        constant_flow = True
    else:
        constant_flow = False

    if constant_flow:
        # Positive time check to filter out bad initial lines on Slocum data.
        valid = (timestamp > 0) & ~(np.isnan(raw))
        # Time lag parameter is a constant scalar.
        tau = params
    else:
        # Positive time check to filter out bad initial lines on Slocum data.
        valid = (timestamp > 0) & ~(np.isnan(raw) | np.isnan(flow))
        # Compute dynamic time lag parameter inversely proportional to flow speed.
        tau_offset = params[0]
        tau_slope = params[1]
        tau = tau_offset + np.divide(tau_slope, flow[valid])

    cor = np.full(np.shape(raw), np.nan)

    # get the indices of the non-nan data
    timestamp_valid = timestamp[valid]
    raw_valid = raw[valid]
    timestamp_unique = timestamp_valid.unique()

    if len(timestamp) > 1:
    # del = [0; diff(raw_valid) ./ diff(timestamp_valid)];
    # cor(valid) = raw_val + tau .* del;
        f = interp1d(timestamp_unique, raw_valid, 'linear', fill_value='extrapolate')
        # cor = f(timestamp_valid.add(tau))

        # insert the interpolated data back into the original data,
        # preserving nans and the length of the original dataset
        cor[valid] = f(timestamp_valid.add(tau))

    return cor


def correct_thermal_lag(timestamp, cond_inside, temp_outside, params, flow_speed=0):
    """
    Authors: Daniel Wang and Jack Slater, Virginia Institute of Marine Science
    Python translation of SOCIB Glider Toolbox correctThermalLag function
    Correct CTD conductivity and temperature sequence from thermal lag effects.
    Note: This function is a recoding of the function by Tomeu Garau
    :param timestamp: sample timestamp
    :param cond_inside: conductivity inside CTD cell
    :param temp_outside: temperature outside CTD cell
    :param params: [error magnitude (alpha), error time constant (tau)]
    :returns
        temp_inside: temperature inside CTD cell
        cond_outside: conductivity outside CTD cell

    References:
        Garau, B.; Ruiz, S.; G. Zhang, W.; Pascual, A.; Heslop, E.; Kerfoot, J.; and Tintoré, J.; 2011:
        Thermal Lag Correction on Slocum CTD Glider Data.
        Journal of Atmospheric and Oceanic Technology, vol. 28, pages 1065-1071.

        Morison, J.; Andersen, R.; Larson, N.; D'Asaro, E.; and Boyd, T.; 1994:
        The Correction for Thermal-Lag Effects in Sea-Bird CTD Data.
        Journal of Atmospheric and Oceanic Technology, vol. 11, pages 1151-1164.

        Lueck, R. G.; and Picklo, J. J.; 1990:
        Thermal Inertia of Conductivity Cells: Observations with a Sea-Bird cell.
        Journal of Atmospheric and Oceanic Technology, vol. 7, pages 756–768.

        Lueck, R. G.; 1990:
        Thermal Inertia of Conductivity Cells: Theory.
        Journal of Atmospheric and Oceanic Technology, vol. 7, pages 741–755.

    Original Authors:
        Joan Pau Beltran  <joanpau.beltran@socib.cat>

        Copyright (C) 2013-2016
        ICTS SOCIB - Servei d'observacio i prediccio costaner de les Illes Balears
        <http://www.socib.es>

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.
    """

    if flow_speed == 0:
        constant_flow = True
    else:
        constant_flow = False

    if constant_flow:
        valid = (timestamp > 0) & ~(np.isnan(cond_inside)) & ~(np.isnan(temp_outside))
        time_val = timestamp[valid]
        temp_val = temp_outside[valid]
        cond_val = cond_inside[valid]

        alpha = params[0]
        tau = params[1]
    else:
        valid = (timestamp > 0) & ~(np.isnan(cond_inside)) & ~(np.isnan(temp_outside) & ~(np.isnan(flow_speed)))
        time_val = timestamp[valid]
        temp_val = temp_outside[valid]
        cond_val = cond_inside[valid]
        flow_val = flow_speed[valid]

        alpha_offset = params[1]
        alpha_slope = params[2]
        tau_offset = params[3]
        tau_slope = params[4]

        # Compute dynamic thermal error and error time parameters for variable flow
        # speed. The formula is given in the references above (Morison 1994).

        alpha = alpha_offset + np.divide(alpha_slope, flow_val[1:-2])
        tau = tau_offset + np.divide(tau_slope, np.sqrt(flow_val[1:-2]))

    # Compute the coefficients of the correction formula.
    # Definitions in references use the Nyquist frequency (half the sampling
    # frequency). This might be wrong in the original implementation by Tomeu
    # Garau, where the sampling frequency was used.
    # These are three equivalent formulas for coefficients.

    dtime = np.diff(time_val)

    # sampling_freq = 1 ./ dtime;
    # nyquist_freq = 0.5 * sampling_freq;
    # coefa = alpha .* (4 * nyquist_freq .* tau) ./ (1 + 4 * nyquist_freq .* tau);
    # coefb = 1 - 2  * (4 * nyquist_freq .* tau) ./ (1 + 4 * nyquist_freq .* tau);
    # coefa = 2 * alpha ./ (2 + dtime .* beta); % from SBE Data Processing.
    # coefb = 1 - 2 .* coefa ./ alpha;

    coefa = 2 * np.divide(alpha, (2 + np.divide(dtime, tau)))  # same, but using tau instead of beta.
    coefb = 1 - np.divide(4, (2 + np.divide(dtime, tau)))

    # Haixing. Follow formula on page93 of SeaBird manual-Seassoft_DataProcessing_7.26.8.pdf;
    # John Kerfoot 2019 poster uses this formula as well.
    dcond_dtemp = np.multiply(0.1, (1 + np.multiply(0.006, (temp_val - 20))))

    # Compute auxiliary vector of consecutive temperature differences.

    dtemp = np.diff(temp_val)

    cond_correction = np.zeros_like(time_val, dtype=float)
    temp_correction = np.zeros_like(time_val, dtype=float)

    for n in range(len(time_val) - 1):
        cond_correction[n + 1] = -coefb[n] * cond_correction[n] + coefa[n] * dcond_dtemp[n] * dtemp[n]
        temp_correction[n + 1] = -coefb[n] * temp_correction[n] + coefa[n] * dtemp[n]

    temp_inside = np.empty_like(timestamp) * np.nan
    cond_outside = np.empty_like(timestamp) * np.nan

    # Apply corrections to valid values in original sequences,
    # preserving invalid values in the output.

    temp_inside[valid] = temp_val - temp_correction
    cond_outside[valid] = cond_val + cond_correction

    return temp_inside, cond_outside


def find_thermal_lag_params_sp(time1, cond1, temp1, pres1, thermocline_pres1, time2, cond2, temp2, pres2,
                               thermocline_pres2, flow1=0, flow2=0):
    """
    Authors: Daniel Wang and Jack Slater, Virginia Institute of Marine Science
    Python translation of find thermal lag params SP
    Thermal lag correction in Pressure (depth) - Salinity space, adjusted to thermocline depth
    """

    if (flow1 & flow2) == 0:
        constant_flow = True
    else:
        constant_flow = False

    def optimobjArea(params):
        """
        Compute area enclosed by profiles in TS diagram after thermal lag correction with supplied alpha and tau
        :param params: [alpha, tau]
        """

        # correct conductivity for thermal lag (calculate conductivity outside of the cell - aligned with thermistor)
        if constant_flow:
            cond_cor1 = correct_thermal_lag(time1, cond1, temp1, params)[1]
            cond_cor2 = correct_thermal_lag(time2, cond2, temp2, params)[1]
        else:
            cond_cor1 = correct_thermal_lag(time1, cond1, temp1, flow1, params)[1]
            cond_cor2 = correct_thermal_lag(time2, cond2, temp2, flow1, params)[1]

        # calculate practical salinity outside of conductivity cell (aligned with thermistor)
        salt_cor1 = gsw.SP_from_C(np.multiply(cond_cor1, 10), temp1, pres1)
        salt_cor2 = gsw.SP_from_C(np.multiply(cond_cor2, 10), temp2, pres2)

        # calculate corrected density
        dens_cor1 = gsw.rho(salt_cor1, temp1, pres1)
        dens_cor2 = gsw.rho(salt_cor2, temp2, pres2)

        #dens_min = np.maximum(np.amin(dens_cor1), np.amin(dens_cor2))
        #dens_max = np.minimum(np.amax(dens_cor1), np.amax(dens_cor2))
        dens_min = np.minimum(np.nanmin(dens_cor1), np.nanmin(dens_cor2))
        dens_max = np.maximum(np.nanmax(dens_cor1), np.nanmax(dens_cor2))

        dens_mask1 = (dens_min <= dens_cor1) & (dens_cor1 <= dens_max)
        dens_mask2 = (dens_min <= dens_cor2) & (dens_cor2 <= dens_max)

        min_idx1 = np.argwhere(dens_mask1)[0][0]
        min_idx2 = np.argwhere(dens_mask2)[0][0]
        max_idx1 = np.argwhere(dens_mask1)[-1][0]
        max_idx2 = np.argwhere(dens_mask2)[-1][0]

        # rescaled Salt-Pressure area
        # pressure is adjusted to thermocline depth
        # outside cond cell
        salt_max = np.maximum(np.nanmax(salt_cor1), np.nanmax(salt_cor2))
        salt_min = np.minimum(np.nanmin(salt_cor1), np.nanmin(salt_cor2))

        pressure_max = np.maximum(np.nanmax(pres1 - thermocline_pres1), np.nanmax(pres2 - thermocline_pres2))
        pressure_min = np.minimum(np.nanmin(pres1 - thermocline_pres1), np.nanmin(pres2 - thermocline_pres2))

        salt1_goodid = salt_cor1[min_idx1:max_idx1 + 1]
        salt2_goodid = salt_cor2[min_idx2:max_idx2 + 1]

        nrmlsalt1 = (salt1_goodid - salt_min) / (salt_max - salt_min)
        nrmlsalt2 = (salt2_goodid - salt_min) / (salt_max - salt_min)

        pres1_goodid = pres1[min_idx1:max_idx1 + 1]
        pres2_goodid = pres2[min_idx2:max_idx2 + 1]

        nrmlpres1 = ((pres1_goodid - thermocline_pres1) - pressure_min) / (pressure_max - pressure_min)
        nrmlpres2 = ((pres2_goodid - thermocline_pres2) - pressure_min) / (pressure_max - pressure_min)

        salt = np.append(nrmlsalt1, nrmlsalt2)
        pres = np.append(nrmlpres1, nrmlpres2)

        # remove nans
        nonan = np.where(~np.isnan(salt))[0]
        salt = salt[nonan]
        pres = pres[nonan]

        nonan = np.where(~np.isnan(pres))[0]
        salt = salt[nonan]
        pres = pres[nonan]

        points = np.concatenate([salt[:, None], pres[:, None]], axis=1)

        outline = spg.Polygon(points)
        area = outline.area

        return area

    params_guess = [0.0677, 11.1431]  # [alpha, tau] for GPCTD
    params = scop.minimize(optimobjArea, params_guess, method='SLSQP', tol=1e-4)
    return params


def find_thermal_lag_params_ts(time1, cond1, temp1, pres1, time2, cond2, temp2, pres2, flow1=0, flow2=0):
    """
    Authors: Daniel Wang and Jack Slater, Virginia Institute of Marine Science
    Python translation of find thermal lag params TS
    Thermal lag correction in T/S (or normalized T/S) space. Not adjusted for thermocline depth.
    """

    if (flow1 & flow2) == 0:
        constant_flow = True
    else:
        constant_flow = False

    def optimobjArea(params):
        """
        Compute area enclosed by profiles in TS diagram after thermal lag correction with supplied alpha and tau
        :param params: [alpha, tau]
        """

        # correct conductivity for thermal lag (calculate conductivity outside of the cell - aligned with thermistor)
        if constant_flow:
            cond_cor1 = correct_thermal_lag(time1, cond1, temp1, params)[1]
            cond_cor2 = correct_thermal_lag(time2, cond2, temp2, params)[1]
        else:
            cond_cor1 = correct_thermal_lag(time1, cond1, temp1, flow1, params)[1]
            cond_cor2 = correct_thermal_lag(time2, cond2, temp2, flow1, params)[1]

        # calculate practical salinity outside of conductivity cell (aligned with thermistor)
        salt_cor1 = gsw.SP_from_C(np.multiply(cond_cor1, 10), temp1, pres1)
        salt_cor2 = gsw.SP_from_C(np.multiply(cond_cor2, 10), temp2, pres2)

        # calculate corrected density
        dens_cor1 = gsw.rho(salt_cor1, temp1, pres1)
        dens_cor2 = gsw.rho(salt_cor2, temp2, pres2)

        #dens_min = np.maximum(np.amin(dens_cor1), np.amin(dens_cor2))
        #dens_max = np.minimum(np.amax(dens_cor1), np.amax(dens_cor2))
        dens_min = np.minimum(np.nanmin(dens_cor1), np.nanmin(dens_cor2))
        dens_max = np.maximum(np.nanmax(dens_cor1), np.nanmax(dens_cor2))

        dens_mask1 = (dens_min <= dens_cor1) & (dens_cor1 <= dens_max)
        dens_mask2 = (dens_min <= dens_cor2) & (dens_cor2 <= dens_max)

        min_idx1 = np.argwhere(dens_mask1)[0][0]
        min_idx2 = np.argwhere(dens_mask2)[0][0]
        max_idx1 = np.argwhere(dens_mask1)[-1][0]
        max_idx2 = np.argwhere(dens_mask2)[-1][0]

        salt_max = np.maximum(np.nanmax(salt_cor1), np.nanmax(salt_cor2))
        salt_min = np.minimum(np.nanmin(salt_cor1), np.nanmin(salt_cor2))

        temp_max = np.maximum(np.nanmax(temp1), np.nanmax(temp2))
        temp_min = np.minimum(np.nanmin(temp1), np.nanmin(temp2))

        salt1_goodid = salt_cor1[min_idx1:max_idx1 + 1]
        salt2_goodid = salt_cor2[min_idx2:max_idx2 + 1]

        nrmlsalt1 = (salt1_goodid - salt_min) / (salt_max - salt_min)
        nrmlsalt2 = (salt2_goodid - salt_min) / (salt_max - salt_min)

        temp1_goodid = temp1[min_idx1:max_idx1 + 1]
        temp2_goodid = temp2[min_idx2:max_idx2 + 1]

        nrmltemp1 = (temp1_goodid - temp_min) / (temp_max - temp_min)
        nrmltemp2 = (temp2_goodid - temp_min) / (temp_max - temp_min)

        salt = np.append(nrmlsalt1, nrmlsalt2)
        temp = np.append(nrmltemp1, nrmltemp2)

        # remove nans, 0.0 and 1.0
        nonan = (~np.isnan(salt)) & (salt > 0.0) & (salt < 1.0)
        salt = salt[nonan]
        temp = temp[nonan]

        nonan = (~np.isnan(temp)) & (temp > 0.0) & (temp < 1.0)
        salt = salt[nonan]
        temp = temp[nonan]

        points = np.concatenate([salt[:, None], temp[:, None]], axis=1)

        outline = spg.Polygon(points)
        area = outline.area

        return area

    params_guess = [0.0677, 11.1431]  # [alpha, tau] for GPCTD
    params = scop.minimize(optimobjArea, params_guess, method='SLSQP', tol=1e-4)
    return params