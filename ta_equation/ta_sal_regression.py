#!/usr/bin/env python

"""
Author: Lori Garzio on 4/6/2021
Last modified: 4/5/2022
Modified from MATLAB code written by Liza Wright-Fairbanks.
Calculate the linear relationship between Total Alkalinity and salinity from in-situ water sampling data taken during
glider deployment and recovery (for nearshore samples = lower salinity region) and historical ECOA water sampling data
(for the offshore = higher salinity region).
In-situ TA data are in units of umol/kg.
ECOA discrete data comes from discrete sample TA analysis during the 2015 East Coast Ocean Acidification cruise.
Samples were collected along 3 lines on the Mid-Atlantic shelf. Only ECOA samples taken at depths shallower than 200 m
are included (max glider depth)
"""
import os
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import functions.common as cf


def main(deploy, deployment_season, deployment_region, sDir):
    ecoa_data = '/Users/garzio/Documents/repo/lgarzio/phglider/water_sampling/ECOA.csv'
    glider_ws_data = f'/Users/garzio/Documents/repo/lgarzio/phglider/water_sampling/ph_water_sampling_{deployment_region}.csv'
    
    ecoa_df = pd.read_csv(ecoa_data)
    glider_df = pd.read_csv(glider_ws_data)
    glider_df = glider_df[glider_df.season == deployment_season]
    glider_sub = glider_df[['discrete_sal', 'discrete_ta']]  # TA units = umol/kg
    df = ecoa_df.append(glider_sub)
    df.dropna(inplace=True)

    # discrete_sal = np.append(np.array(ecoa_df.discrete_sal), np.array(glider_df.discrete_sal))
    # discrete_ta = np.append(np.array(ecoa_df.discrete_ta), np.array(glider_df.discrete_ta))  # units = umol/kg

    r_sq, m, b, y_predict = cf.linear_regression(np.array(df.discrete_sal), np.array(df.discrete_ta))
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
    ax.plot(ecoa_df.discrete_sal, ecoa_df.discrete_ta, 'r.', ms=5, label='ECOA')
    ax.plot(glider_df.discrete_sal, glider_df.discrete_ta, 'b.', ms=5, label=glider_label)
    ax.plot(df.discrete_sal, y_predict, 'k-', label=equation)
    ax.text(np.nanmin(df.discrete_sal), np.nanmax(df.discrete_ta) * 0.97, 'R2: {}'.format(np.round(r_sq, 2)), fontsize=8)
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Total Alkalinity')
    plt.legend(fontsize=8)
    plt.title(deploy)
    sname = os.path.join(sDir, '{}_ta_sal_regression.png'.format(deploy))
    plt.savefig(sname, dpi=200)
    plt.close()

    values = {'m': m,
              'b': b
              }

    fname = '{}_ta_equation.txt'.format(deploy)
    with open(fname, 'w') as outfile:
        json.dump(values, outfile)


if __name__ == '__main__':
    deployment = 'ru30-20210503T1929'
    season = 'spring'
    region = 'mab'
    save_dir = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210503T1929/delayed/'
    main(deployment, season, region, save_dir)
