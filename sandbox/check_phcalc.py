#!/usr/bin/env python

"""
Author: Lori Garzio on 4/18/2023
Last modified: 4/28/2023
Quick check of the pH calculation
"""

import json
import numpy as np
from gsw import SP_from_C
import functions.phcalc as phcalc


def calculate_practical_salinity(conductivity, temperature, pressure):
    """Calculates practical salinity given glider conductivity, temperature,
    and pressure using Gibbs gsw SP_from_C function.
    Parameters:
        timestamp, conductivity (S/m), temperature (C), and pressure (bar).
    Returns:
        salinity (psu PSS-78).
    """

    # Convert S/m to mS/cm
    ms_conductivity = conductivity * 10

    return SP_from_C(
        ms_conductivity,
        temperature,
        pressure
    )

# pH in the tank = 7.99 with these values
volt = -0.81229
pres_dbar = 0.12
cond = 4.3229
temp = 21.3003
sal = calculate_practical_salinity(cond, temp, pres_dbar/10)

calfile = '/Users/garzio/Documents/repo/lgarzio/phglider/calibration/sbe12110_20221021.txt'
with open(calfile) as json_file:
    cc = json.load(json_file)

try:
    f_p = np.polyval([cc['f12'], cc['f11'], cc['f10'], cc['f9'], cc['f8'], cc['f7'], cc['f6'], cc['f5'], cc['f4'],
                      cc['f3'], cc['f2'], cc['f1'], 0], pres_dbar)
    k2 = [cc['k2f3'], cc['k2f2'], cc['k2f1'], cc['k2f0']]
except KeyError:
    f_p = np.polyval([cc['f6'], cc['f5'], cc['f4'], cc['f3'], cc['f2'], cc['f1'], 0], pres_dbar)
    k2 = cc['k2']
phfree, phtot = phcalc.phcalc(volt, pres_dbar, temp, sal, cc['k0'], k2, f_p)

print(phtot)
