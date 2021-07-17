#!/usr/bin/env python

"""
Author: Lori Garzio on 3/2/2021
Last modified: 7/17/2021
Plot realtime pH glider data variables: seawater temperature, salinity, chlorophyll, dissolved oxygen, pH reference
voltage, and pH (not corrected for time lag)
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


#def main(deploy, sensor_sn, sfilename):
def main(args):
    ru_server = 'http://slocum-data.marine.rutgers.edu//erddap'
    deploy = args.deployment
    print('\nPlotting {}'.format(deploy))
    glider_id = '{}-profile-sci-rt'.format(deploy)
    glider_vars = ['latitude', 'longitude', 'depth', 'conductivity', 'salinity', 'sci_water_pressure',
                   'temperature', 'sbe41n_ph_ref_voltage', 'oxygen_concentration', 'water_depth']
    if 'um_242' in deploy:
        glider_vars.append('sci_flntu_chlor_units')
    else:
        glider_vars.append('chlorophyll_a')

    gargs = dict()
    gargs['variables'] = glider_vars
    ds = cf.get_erddap_dataset(ru_server, glider_id, **gargs)
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)

    # read cal file
    # calfile = cf.find_calfile(deploy, sensor_sn)
    calfile = cf.find_calfile(deploy, args.sensor_sn)
    with open(calfile) as json_file:
        cc = json.load(json_file)

    tm = ds.time.values
    depth = ds.depth.values
    pressure_dbar = ds.sci_water_pressure.values * 10
    temp = ds.temperature.values
    sal = ds.salinity.values
    vrs = ds.sbe41n_ph_ref_voltage.values
    vrs[vrs == 0.0] = np.nan  # convert zero values to nan
    oxy = (ds.oxygen_concentration.values * 32) / 1000  # change oxygen from umol/L to mg/L
    oxy[oxy == 0.0] = np.nan  # convert zero values to nan
    chl = ds.chlorophyll_a.values
    chl[chl == 0.0] = np.nan  # convert zero values to nan

    # calculate pH
    ph = np.array([])
    for i, press in enumerate(pressure_dbar):
        f_p = np.polyval([cc['f6'], cc['f5'], cc['f4'], cc['f3'], cc['f2'], cc['f1'], 0], press)
        phfree, phtot = phcalc.phcalc(vrs[i], press, temp[i], sal[i], cc['k0'], cc['k2'], f_p)
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
        t_start = pd.to_datetime(tm[0])
        t_end = pd.to_datetime(tm[-1])
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
    main_ttl = '{}: deployed {} (updated {})'.format(glider, deploy_date.strftime('%d-%b-%Y'),
                                                     dt.datetime.now().strftime('%d-%b-%Y %H:%M'))
    fig.suptitle(main_ttl, fontsize=22, y=.96)

    plt.subplots_adjust(left=0.08, right=0.91, bottom=0.1, top=0.88)
    #plt.savefig(sfilename, dpi=300)
    plt.savefig(args.sfilename, dpi=300)
    plt.close()


if __name__ == '__main__':
    # deployment = 'ru30-20210503T1929'
    # ph_sn = 'sbe10344'
    # savefile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210503T1929/rt_plotting/phoxy_live.png'
    # deployment = 'sbu01-20210226T1902'
    # ph_sn = 'sbe10528'
    # savefile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/sbu01-20210226T1902/phoxy_live_sbu01-test.png'
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
