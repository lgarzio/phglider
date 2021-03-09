#!/usr/bin/env python

"""
Author: Lori Garzio on 3/2/2021
Last modified: 3/8/2021
Plot realtime glider data
"""

import argparse
import numpy as np
import datetime as dt
import json
from erddapy import ERDDAP
import cmocean as cmo
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import functions.phcalc as phcalc
import functions.common as cf
plt.rcParams.update({'font.size': 16})


def get_erddap_dataset(server, protocol, file_type, ds_id, var_list=None):
    e = ERDDAP(server=server,
               protocol=protocol,
               response=file_type)
    e.dataset_id = ds_id
    if var_list:
        e.variables = var_list
    ds = e.to_xarray()
    ds = ds.sortby(ds.time)
    return ds


def main(deploy, sensor_sn, sfilename):
    ru_server = 'http://slocum-data.marine.rutgers.edu//erddap'
    id = '{}-profile-sci-rt'.format(deploy)
    glider_vars = ['latitude', 'longitude', 'depth', 'conductivity', 'salinity', 'sci_water_pressure',
                   'temperature', 'sbe41n_ph_ref_voltage', 'chlorophyll_a', 'oxygen_concentration', 'water_depth']
    ds = get_erddap_dataset(ru_server, 'tabledap', 'nc', id, glider_vars)
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)

    # read cal file
    calfile = cf.find_calfile(sensor_sn)
    with open(calfile) as json_file:
        cc = json.load(json_file)

    tm = ds.time.values
    depth = ds.depth.values
    pressure_dbar = ds.sci_water_pressure.values * 10
    temp = ds.temperature.values
    sal = ds.salinity.values
    vrs = ds.sbe41n_ph_ref_voltage.values
    vrs[vrs == 0.0] = np.nan  # convert voltages of zero to nan
    oxy = (ds.oxygen_concentration.values * 32) / 1000  # change oxygen from umol/L to mg/L
    chl = ds.chlorophyll_a.values

    ph = np.array([])
    for i, press in enumerate(pressure_dbar):
        f_P = np.polyval([cc['f6'], cc['f5'], cc['f4'], cc['f3'], cc['f2'], cc['f1'], 0], press)
        phfree, phtot = phcalc.phcalc(vrs[i], press, temp[i], sal[i], cc['k0'], cc['k2'], f_P)
        if phtot > 14:
            phtot = np.nan
        ph = np.append(ph, phtot)

    # plot
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(20, 10), sharex=True, sharey=True)
    plt_vars = {'temp': {'var': temp, 'axes': ax1, 'cmap': cmo.cm.thermal, 'ttl': 'Temperature ({})'.format(r'$\rm ^oC$')},
                'salinity': {'var': sal, 'axes': ax4, 'cmap': cmo.cm.haline, 'ttl': 'Salinity'},
                'chl': {'var': chl, 'axes': ax2, 'cmap': cmo.cm.algae, 'ttl': 'Chlorophyll ({}g/L)'.format(chr(956))},
                'oxy': {'var': oxy, 'axes': ax5, 'cmap': cmo.cm.oxy, 'ttl': 'Oxygen (mg/L)'},
                'ph_volt': {'var': vrs, 'axes': ax3, 'cmap': cmo.cm.matter, 'ttl': 'pH Reference Voltage'},
                'ph': {'var': ph, 'axes': ax6, 'cmap': cmo.cm.matter, 'ttl': 'pH (uncorrected)'}
                }
    for pv, info in plt_vars.items():
        xc = info['axes'].scatter(tm, depth, c=info['var'], cmap=info['cmap'], s=10, edgecolor='None')

        info['axes'].set_title(info['ttl'], fontsize=18)

        # format colorbar
        divider = make_axes_locatable(info['axes'])
        cax = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
        fig.add_axes(cax)
        cb = plt.colorbar(xc, cax=cax)

        # format x-axis
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
    main_ttl = '{}: deployed {} (updated {})'.format(glider, deploy_date.strftime('%d-%b-%Y'),
                                                     dt.datetime.now().strftime('%d-%b-%Y %H:%M'))
    fig.suptitle(main_ttl, fontsize=22, y=.96)

    plt.subplots_adjust(left=0.08, right=0.91, bottom=0.1, top=0.88)
    plt.savefig(sfilename, dpi=300)
    plt.close()


if __name__ == '__main__':
    deployment = 'ru30-20210226T1647'
    #cal = '/Users/lgarzio/Documents/repo/lgarzio/phglider/calibration/sbe10344_20200306.txt'
    ph_sn = 'sbe10344'
    savefile = '/Users/lgarzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/phoxy_live-test.png'
    main(deployment, ph_sn, savefile)
    # arg_parser = argparse.ArgumentParser(description='Plot real time glider pH data',
    #                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #
    # arg_parser.add_argument('-deployment',
    #                         dest='deploy',
    #                         default='ru30-20210226T1647',
    #                         type=str,
    #                         help='Glider deployment to plot. e.g. glider-yyyymmddTHHMM')
    #
    # arg_parser.add_argument('-calfile',
    #                         dest='calfile',
    #                         default='/Users/lgarzio/Documents/repo/lgarzio/phglider/calibration/sbe10344_20200306.txt',
    #                         type=str,
    #                         help='Calibration file for SBE pH sensor.')
    #
    # parsed_args = arg_parser.parse_args()
    # sys.exit(main(parsed_args))
