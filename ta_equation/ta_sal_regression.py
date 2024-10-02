#!/usr/bin/env python

"""
Author: Lori Garzio on 4/6/2021
Last modified: 10/2/2024
Modified from MATLAB code written by Liza Wright-Fairbanks.
Calculate the linear relationship between Total Alkalinity and salinity from in-situ water sampling data taken during
glider deployment and recovery (for nearshore samples = lower salinity region) and historical vessel-based water
sampling data (for the offshore = higher salinity region).
In-situ TA data are in units of umol/kg.
Vessel-based data comes from discrete sample TA analysis during vessel-based oceanographic cruises. Data were taken from
the CODAP-NA dataset (documented here: https://essd.copernicus.org/articles/13/2777/2021/), and more recent datasets
from ECOMON and ECOA cruises were downloaded from the NCEI OCADs data portal
(https://www.ncei.noaa.gov/products/ocean-carbon-acidification-data-system). Data for this analysis were restricted
to a defined bounding box spanning only the New York Bight and limited to 200m depth (max glider depth)
"""
import os
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import functions.common as cf


def main(deployment_season, deployment_region, sDir):
    vessel_data = f'/Users/garzio/Documents/repo/lgarzio/phglider/water_sampling/vessel_based_TA_salinity_{deployment_region}.csv'
    glider_ws_data = f'/Users/garzio/Documents/repo/lgarzio/phglider/water_sampling/ph_water_sampling_{deployment_region}.csv'
    
    # grab vessel data for specified season
    vessel_df = pd.read_csv(vessel_data)
    vessel_df = vessel_df[vessel_df.season == deployment_season]
    vessel_sub = vessel_df[['salinity', 'total_alkalinity']]  # TA units = umol/kg

    # grab water sampling done during glider ops for specified season
    glider_df = pd.read_csv(glider_ws_data)
    glider_df = glider_df[glider_df.season == deployment_season]
    glider_df.dropna(subset='depth', inplace=True)
    glider_sub = glider_df[['salinity', 'total_alkalinity']]  # TA units = umol/kg

    # combine dataframes
    df = vessel_sub.append(glider_sub)
    df.dropna(inplace=True)

    r_sq, m, b, y_predict = cf.linear_regression(np.array(df.salinity), np.array(df.total_alkalinity))
    equation = 'y = {}x + {}'.format(np.round(m, 2), np.round(b, 2))

    source = np.unique(glider_df['source']).tolist()
    if len(source) == 1:
        glider_label = '{}'.format(source[0].upper())
    elif len(source) == 2:
        glider_label = '{}/{}'.format(source[0].upper(), source[1].upper())
    elif len(source) == 3:
        glider_label = '{}/{}/{}'.format(source[0].upper(), source[1].upper(), source[2].upper())

    # plot
    fig, ax = plt.subplots()
    vlab = f'Vessel (n={len(vessel_sub.salinity)})'
    glab = f'{glider_label} (n={len(glider_sub.salinity)})'
    ax.plot(vessel_sub.salinity, vessel_sub.total_alkalinity, 'r.', ms=5, label=vlab)
    ax.plot(glider_sub.salinity, glider_sub.total_alkalinity, 'b.', ms=5, label=glab)
    ax.plot(df.salinity, y_predict, 'k-', label=equation)
    #ax.text(np.nanmin(df.salinity), np.nanmax(df.total_alkalinity) * 0.9, 'R2: {}'.format(np.round(r_sq, 2)), fontsize=8)
    ax.set_ylim([1900, 2400])
    ax.set_xlim([28, 37])
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Total Alkalinity')
    plt.legend(fontsize=8)
    plt.title(f'{deployment_season} (R2={np.round(r_sq, 2)})')
    sname = os.path.join(sDir, f'{deployment_season}_ta_sal_regression_{deployment_region}.png')
    plt.savefig(sname, dpi=200)
    plt.close()

    values = {'m': np.round(m, 2),
              'b': np.round(b, 2)
              }

    fname = f'{deployment_season}_ta_equation_{deployment_region}.txt'
    with open(fname, 'w') as outfile:
        json.dump(values, outfile)


if __name__ == '__main__':
    season = 'SON'
    region = 'NYB'
    save_dir = '/Users/garzio/Documents/repo/lgarzio/phglider/water_sampling/regression_plots'
    main(season, region, save_dir)
