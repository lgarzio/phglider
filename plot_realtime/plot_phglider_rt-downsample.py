#!/usr/bin/env python

"""
Author: Lori Garzio on 3/2/2021
Last modified: 11/11/2021
Plot realtime pH glider data variables: seawater temperature, salinity, chlorophyll, dissolved oxygen, pH reference
voltage, and pH (not corrected for time lag). Data are downsampled to line up mismatched timestamps between CTD and
other instruments
"""

import argparse
import sys
import numpy as np
import pandas as pd
import datetime as dt
import json
import cmocean as cmo
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import functions.phcalc as phcalc
import functions.common as cf
plt.rcParams.update({'font.size': 16})


# def main(deploy, sensor_sn, sfilename):
def main(args):
    ru_server = 'http://slocum-data.marine.rutgers.edu//erddap'
    deploy = args.deployment
    print('\nPlotting {}'.format(deploy))
    glider_id = '{}-profile-sci-rt'.format(deploy)
    glider_vars = ['latitude', 'longitude', 'depth', 'conductivity', 'salinity', 'pressure',
                   'temperature', 'sbe41n_ph_ref_voltage', 'oxygen_concentration', 'water_depth']
    if 'um_242' in deploy:
        chlvar = 'sci_flntu_chlor_units'
    else:
        chlvar = 'chlorophyll_a'
    glider_vars.append(chlvar)

    gargs = dict()
    gargs['variables'] = glider_vars
    ds = cf.get_erddap_dataset(ru_server, glider_id, **gargs)
    ds = ds.drop_vars(names=['profile_id', 'rowSize'])
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)

    # read cal file
    # calfile = cf.find_calfile(deploy, sensor_sn)
    calfile = cf.find_calfile(deploy, args.sensor_sn)
    with open(calfile) as json_file:
        cc = json.load(json_file)

    for v in ['conductivity', 'temperature', 'sbe41n_ph_ref_voltage', 'oxygen_concentration', chlvar]:
        ds[v][ds[v] == 0.0] = np.nan  # convert zeros to nan

    # resample the dataframe to 20s to line up possible mismatched timestamps between CTD and other instruments
    df = ds.to_dataframe()
    df = df.resample('20s').mean()
    col_subset = ['conductivity', 'salinity', 'pressure', 'temperature', 'sbe41n_ph_ref_voltage',
                  'oxygen_concentration', chlvar]

    # convert oxygen from umol/L to mg/L
    df['oxygen_concentration'] = (df['oxygen_concentration'] * 32) / 1000

    # drop timestamps with duplicated and missing data and any rows where pressure <1 dbar
    df = df.dropna(axis=0, how='all', subset=col_subset)
    df.drop_duplicates(inplace=True)
    df = df[df.pressure > 1]

    # calculate pH and add to dataframe
    df['f_p'] = np.polyval([cc['f6'], cc['f5'], cc['f4'], cc['f3'], cc['f2'], cc['f1'], 0], df.pressure)
    phfree, phtot = phcalc.phcalc(df.sbe41n_ph_ref_voltage, df.pressure, df.temperature, df.salinity,
                                  cc['k0'], cc['k2'], df.f_p)
    df['ph_total'] = phtot

    # plot
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(20, 10), sharex=True, sharey=True)
    plt_vars = {'temp': {'var': df.temperature, 'axes': ax1, 'cmap': cmo.cm.thermal, 'ttl': 'Temperature ({})'.format(r'$\rm ^oC$')},
                'salinity': {'var': df.salinity, 'axes': ax4, 'cmap': cmo.cm.haline, 'ttl': 'Salinity'},
                'chl': {'var': df[chlvar], 'axes': ax2, 'cmap': cmo.cm.algae, 'ttl': 'Chlorophyll ({}g/L)'.format(chr(956))},
                'oxy': {'var': df.oxygen_concentration, 'axes': ax5, 'cmap': cmo.cm.oxy, 'ttl': 'Oxygen (mg/L)'},
                'ph_volt': {'var': df.sbe41n_ph_ref_voltage, 'axes': ax3, 'cmap': cmo.cm.matter, 'ttl': 'pH Reference Voltage'},
                'ph': {'var': df.ph_total, 'axes': ax6, 'cmap': cmo.cm.matter, 'ttl': 'pH (uncorrected)'}
                }
    for pv, info in plt_vars.items():
        xc = info['axes'].scatter(df.index, df.depth, c=info['var'], cmap=info['cmap'], s=10, edgecolor='None')

        info['axes'].set_title(info['ttl'], fontsize=18)

        # format colorbar
        divider = make_axes_locatable(info['axes'])
        cax = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
        fig.add_axes(cax)
        cb = plt.colorbar(xc, cax=cax)

        # format x-axis
        t_start = pd.to_datetime(df.index[0])
        t_end = pd.to_datetime(df.index[-1])
        if t_end - t_start < dt.timedelta(days=5):
            t0 = t_start
            tf = t_end
        else:
            delta = dt.timedelta(hours=6)
            t0 = t_start + delta
            tf = t_end - delta
        xdates = pd.date_range(t0, tf, periods=6)  # create 6 time bins
        info['axes'].set_xticks(xdates)
        xfmt = mdates.DateFormatter('%d-%b')
        info['axes'].xaxis.set_major_formatter(xfmt)
        info['axes'].xaxis.set_tick_params(labelsize=15)

        # add grid
        info['axes'].grid(ls='--', lw=.5)

    # invert the last iteration through the axes
    info['axes'].invert_yaxis()

    plt_vars['temp']['axes'].set_ylabel('Depth (m)')
    plt_vars['salinity']['axes'].set_ylabel('Depth (m)')

    splitter = deploy.split('-')
    glider = splitter[0]
    deploy_date = dt.datetime.strptime(splitter[1].split('T')[0], '%Y%m%d')
    main_ttl = '{}: deployed {} (updated {} EST)'.format(glider, deploy_date.strftime('%d-%b-%Y'),
                                                     dt.datetime.now().strftime('%d-%b-%Y %H:%M'))
    fig.suptitle(main_ttl, fontsize=22, y=.96)

    plt.subplots_adjust(left=0.08, right=0.91, bottom=0.1, top=0.88)
    #plt.savefig(sfilename, dpi=300)
    plt.savefig(args.sfilename, dpi=300)
    plt.close()


if __name__ == '__main__':
    # deployment = 'ru30-20211020T1316'
    # ph_sn = 'sbeC18'
    # savefile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20211020T1316/rt_plotting/phoxy_live.png'
    # deployment = 'um_242-20211105T1601'
    # ph_sn = 'sbe10490'
    # savefile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/um_242-20211105T1601/um_242-20211105T1601_live-test.png'
    # main(deployment, ph_sn, savefile)
    arg_parser = argparse.ArgumentParser(description='Plot real time glider pH data',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('-deploy',
                            dest='deployment',
                            type=str,
                            help='Glider deployment to plot (glider-yyyymmddTHHMM), e.g. ru30-20210226T1647')

    arg_parser.add_argument('-sn',
                            dest='sensor_sn',
                            type=str,
                            help='pH sensor serial number, e.g. sbe10344')

    arg_parser.add_argument('-s', '--save_file',
                            dest='sfilename',
                            type=str,
                            help='Full file path to save directory and save filename')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
