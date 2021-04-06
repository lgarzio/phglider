#!/usr/bin/env python

"""
Author: Lori Garzio on 4/6/2021
Last modified: 4/6/2021
Modified from MATLAB code written by Liza Wright-Fairbanks.
Calculate the linear relationship between TA and salinity from in-situ water sampling data taken during glider
deployment and recovery (for nearshore samples = lower salinity region) and historical ECOA water sampling data (for
the offshore = higher salinity region).
ECOA discrete data comes from discrete sample TA analysis during the 2015 East Coast Ocean Acidification cruise.
Samples were collected along 3 lines on the Mid-Atlantic shelf. Only samples taken at depths shallower than 200 m are
included (max glider depth)
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import functions.common as cf


def main(deploy, ecoa_discrete, glider_discrete, sDir):
    ecoa_df = pd.read_csv(ecoa_discrete)
    glider_df = pd.read_csv(glider_discrete)

    discrete_sal = np.append(np.array(ecoa_df.discrete_sal), np.array(glider_df.discrete_sal))
    discrete_ta = np.append(np.array(ecoa_df.discrete_ta), np.array(glider_df.discrete_ta))

    r_sq, m, b, y_predict = cf.linear_regression(discrete_sal, discrete_ta)
    equation = 'y = {}x + {}'.format(np.round(m, 2), np.round(b, 2))

    # plot
    fig, ax = plt.subplots()
    ax.plot(ecoa_df.discrete_sal, ecoa_df.discrete_ta, 'r.', ms=5, label='ECOA')
    ax.plot(glider_df.discrete_sal, glider_df.discrete_ta, 'b.', ms=5, label='RU')
    ax.plot(discrete_sal, y_predict, 'k-', label=equation)
    ax.text(np.nanmin(discrete_sal), np.nanmax(discrete_ta) * 0.97, 'R2: {}'.format(np.round(r_sq, 2)), fontsize=8)
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Total Alkalinity')
    plt.legend(fontsize=8)
    plt.title(deploy)
    sname = os.path.join(sDir, '{}_ta_sal_regression.png'.format(deploy))
    plt.savefig(sname, dpi=200)
    plt.close()


if __name__ == '__main__':
    deployment = 'ru30-20210226T1647'
    ecoa_data = '/Users/garzio/Documents/repo/lgarzio/phglider/water_sampling/ECOA.csv'
    glider_discrete_data = '/Users/garzio/Documents/repo/lgarzio/phglider/water_sampling/ru30-20210226T1647.csv'
    save_dir = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/plots'
    main(deployment, ecoa_data, glider_discrete_data, save_dir)
