#!/usr/bin/env python

"""
Author: Lori Garzio on 5/13/2021
Last modified: 5/18/2021
Plot realtime pH glider data from yesterday for the purposes of QC: seawater conductivity, temperature, salinity,
chlorophyll, dissolved oxygen, pH reference voltage, and pH (not corrected for time lag)
"""

import argparse
import sys
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import json
import os
import cmocean as cmo
import matplotlib.pyplot as plt
import functions.phcalc as phcalc
import functions.common as cf
import functions.plotting as pf
plt.rcParams.update({'font.size': 12})


def main(args):
#def main(deploy, sensor_sn, savedir):
    ru_server = 'http://slocum-data.marine.rutgers.edu//erddap'
    deploy = args.deployment
    savedir = args.save_dir
    glider = deploy.split('-')[0]
    glider_id = '{}-profile-sci-rt'.format(deploy)

    # get yesterday's date
    today = dt.date.today()
    #today = dt.date.today() - dt.timedelta(days=10)  # for debugging
    t0 = today - dt.timedelta(days=1)
    t0_savestr = t0.strftime('%Y%m%d')
    t0str = t0.strftime('%Y-%m-%dT00:00:00Z')
    t1str = today.strftime('%Y-%m-%dT00:00:00Z')

    constraints = dict({'time>=': t0str, 'time<=': t1str})
    glider_vars = ['latitude', 'longitude', 'depth', 'conductivity', 'salinity', 'sci_water_pressure',
                   'temperature', 'sbe41n_ph_ref_voltage', 'chlorophyll_a', 'oxygen_concentration', 'water_depth',
                   'profile_time']

    gargs = dict()
    gargs['variables'] = glider_vars
    gargs['constraints'] = constraints
    ds = cf.get_erddap_dataset(ru_server, glider_id, **gargs)
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)

    # drop everything where all key values are nan
    ds = ds.where((~np.isnan(ds.sci_water_pressure) | ~np.isnan(ds.temperature) | ~np.isnan(ds.sbe41n_ph_ref_voltage)
                   | ~np.isnan(ds.chlorophyll_a) | ~np.isnan(ds.oxygen_concentration)), drop=True)

    # read cal file
    #calfile = cf.find_calfile(deploy, sensor_sn)
    calfile = cf.find_calfile(deploy, args.sensor_sn)
    with open(calfile) as json_file:
        cc = json.load(json_file)

    ds['pressure_dbar'] = ds.sci_water_pressure * 10
    for v in ['sbe41n_ph_ref_voltage', 'oxygen_concentration', 'chlorophyll_a']:
        ds[v][ds[v] == 0.0] = np.nan  # convert zeros to nan

    ds['oxygen_concentration_mgL'] = ds.oxygen_concentration * 32 / 1000  # convert oxygen from umol/L to mg/L

    # calculate pH
    ph = np.array([])
    for i, press in enumerate(ds.pressure_dbar.values):
        f_p = np.polyval([cc['f6'], cc['f5'], cc['f4'], cc['f3'], cc['f2'], cc['f1'], 0], press)
        vrs = ds.sbe41n_ph_ref_voltage.values[i]
        temp = ds.temperature.values[i]
        sal = ds.salinity.values[i]
        phfree, phtot = phcalc.phcalc(vrs, press, temp, sal, cc['k0'], cc['k2'], f_p)
        if phtot > 14:
            phtot = np.nan
        ph = np.append(ph, phtot)

    # add pH to the dataset
    da = xr.DataArray(ph, coords=ds.sbe41n_ph_ref_voltage.coords, dims=ds.sbe41n_ph_ref_voltage.dims, name='ph')
    ds['ph'] = da

    # plot
    plt_vars = {'conductivity': {'cmap': cmo.cm.thermal, 'ttl': 'Conductivity (S m-1)'},
                'temperature': {'cmap': cmo.cm.thermal, 'ttl': 'Temperature ({})'.format(r'$\rm ^oC$')},
                'salinity': {'cmap': cmo.cm.haline, 'ttl': 'Salinity'},
                'chlorophyll_a': {'cmap': cmo.cm.algae, 'ttl': 'Chlorophyll ({}g/L)'.format(chr(956))},
                'oxygen_concentration_mgL': {'cmap': cmo.cm.oxy, 'ttl': 'Oxygen (mg/L)'},
                'sbe41n_ph_ref_voltage': {'cmap': cmo.cm.matter, 'ttl': 'pH Reference Voltage'},
                'ph': {'cmap': cmo.cm.matter, 'ttl': 'pH (uncorrected)'}
                }

    savedir_dt = os.path.join(savedir, deploy, t0_savestr)
    os.makedirs(savedir_dt, exist_ok=True)

    # define the down-up profile pairs (yos), profile times and indices for profile plots
    pft, pf_directions, idxs = cf.yos(ds)
    n = len(pf_directions)
    if (n & 1) == 1:  # if the number of casts is odd
        n = n + 1

    yo_idx0 = idxs[0:int(n / 2)]
    yo_idx1 = idxs[int(n / 2):]
    yos = dict()
    yos['0'] = dict()
    yos['0']['idx'] = yo_idx0
    yos['0']['ds_sel'] = [yo_idx0[0][0], yo_idx1[0][0]]
    yos['1'] = dict()
    yos['1']['idx'] = yo_idx1
    yos['1']['ds_sel'] = [yo_idx1[0][0], len(ds.time)]

    for pv, info in plt_vars.items():
        dst0 = pd.to_datetime(ds.time.values[0]).strftime('%Y-%m-%dT%H:%M')
        dst1 = pd.to_datetime(ds.time.values[-1]).strftime('%Y-%m-%dT%H:%M')

        # plot xsection - entire day
        fig, ax = plt.subplots(figsize=(12, 8))
        figttl_xsection = f'{glider} {info["ttl"].split(" (")[0]}\n{dst0} to {dst1}'
        xargs = dict()
        xargs['clabel'] = info['ttl']
        xargs['title'] = figttl_xsection
        xargs['date_fmt'] = '%m-%d\n%H:%M'
        xargs['grid'] = True
        pf.xsection(fig, ax, ds.time.values, ds.depth.values, ds[pv].values, **xargs)

        sfilename = f'{glider}_xsection_{pv}_{t0_savestr}.png'
        sfile = os.path.join(savedir_dt, sfilename)
        plt.savefig(sfile, dpi=300)
        plt.close()

        # plot profiles
        for jj, si in yos.items():
            colors = plt.cm.rainbow(np.linspace(0, 1, len(si['idx'])))
            fig, ax = plt.subplots(figsize=(8, 10))
            da = ds[pv]
            da_sel = da[si['ds_sel'][0]:si['ds_sel'][1]]
            da_t0 = pd.to_datetime(da_sel.time.values[0]).strftime('%Y-%m-%dT%H:%M')
            da_t1 = pd.to_datetime(da_sel.time.values[-1]).strftime('%Y-%m-%dT%H:%M')
            da_t0save = pd.to_datetime(da_sel.time.values[0]).strftime('%Y%m%dT%H%M')

            figttl_profile = f'{glider} {info["ttl"].split(" (")[0]}\n{da_t0} to {da_t1}'
            pargs = dict()
            pargs['xlabel'] = info['ttl']
            pargs['title'] = figttl_profile
            pargs['grid'] = True
            pf.profile_yos(fig, ax, da, si['idx'], colors, **pargs)

            sfilename = f'{glider}_profile_{pv}_{da_t0save}.png'
            sfile = os.path.join(savedir_dt, sfilename)
            plt.savefig(sfile, dpi=300)
            plt.close()


if __name__ == '__main__':
    # deployment = 'ru30-20210503T1929'
    # ph_sn = 'sbe10344'
    # save_dir = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021'
    # main(deployment, ph_sn, save_dir)
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

    arg_parser.add_argument('-s', '--save_dir',
                            dest='save_dir',
                            type=str,
                            help='Full file path to save directory')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
