#!/usr/bin/env python

"""
Author: Lori Garzio on 5/13/2021
Last modified: 8/10/2021
Plot delayed-mode pH glider data variables for the purposes of QC: seawater conductivity, temperature, salinity,
chlorophyll, dissolved oxygen, pH reference voltage, and pH
"""

import numpy as np
import pandas as pd
import xarray as xr
import os
import matplotlib.pyplot as plt
import functions.common as cf
import functions.plotting as pf
plt.rcParams.update({'font.size': 12})


def main(fname, num_yos, dtrange, plot_xsection=None, plot_profiles=None, plot_qcvars=None):
    plot_xsection = plot_xsection or True
    plot_profiles = plot_profiles or True
    plot_qcvars = plot_qcvars or False

    ds = xr.open_dataset(fname)
    ds = ds.sortby(ds.time)
    if dtrange:
        ds = ds.sel(time=slice(dtrange[0], dtrange[1]))
    deploy = '-'.join(fname.split('/')[-1].split('-')[0:2])
    glider = deploy.split('-')[0]
    ds_t0savestr = pd.to_datetime(ds.time.values[0]).strftime('%Y%m%dT%H%M')
    ds_t1savestr = pd.to_datetime(ds.time.values[-1]).strftime('%Y%m%dT%H%M')

    # drop everything where all key values are nan
    ds = ds.where((~np.isnan(ds.temperature) | ~np.isnan(ds.chlorophyll_a) | ~np.isnan(ds.sbe41n_ph_ref_voltage)
                   | ~np.isnan(ds.sbe41n_ph_ref_voltage_shifted) | ~np.isnan(ds.sbe41n_ph_ref_voltage_shifted)
                   | ~np.isnan(ds.oxygen_concentration) | ~np.isnan(ds.oxygen_concentration_shifted)
                   | ~np.isnan(ds.ph_total) | ~np.isnan(ds.ph_total_shifted)), drop=True)

    for v in ['sbe41n_ph_ref_voltage', 'oxygen_concentration', 'chlorophyll_a', 'sbe41n_ph_ref_voltage_shifted',
              'oxygen_concentration_shifted', 'ph_total', 'ph_total_shifted']:
        ds[v][ds[v] == 0.0] = np.nan  # convert zeros to nan

    # plot
    plt_vars = cf.plot_vars()

    sdir_xsection = os.path.join(os.path.dirname(fname), 'qc_figs_ph', f'{ds_t0savestr}_{ds_t1savestr}', 'xsections')
    sdir_profile = os.path.join(os.path.dirname(fname), 'qc_figs_ph', f'{ds_t0savestr}_{ds_t1savestr}', 'profiles')
    os.makedirs(sdir_xsection, exist_ok=True)
    os.makedirs(sdir_profile, exist_ok=True)

    for pv, info in plt_vars.items():
        try:
            ds[pv]
            sdir_xs_var = os.path.join(sdir_xsection, pv)
            sdir_pf_var = os.path.join(sdir_profile, pv)
            os.makedirs(sdir_xs_var, exist_ok=True)
            os.makedirs(sdir_pf_var, exist_ok=True)
        except KeyError:
            continue

    # define the down-up profile pairs (yos), profile times and indices for profile plots
    pft, pf_directions, idxs = cf.yos(ds)
    n = len(pf_directions)

    # divide the yos into user-defined groups for plotting
    divs = np.arange(num_yos, n, num_yos)
    divs = np.append(divs, n)  # include the last endpoint

    # for each group, subset the dataset
    for div in divs:
        div0 = div - num_yos
        div_idxs = idxs[div0:div]
        div_tm = ds.time.values[div_idxs[0][0]: div_idxs[-1][-1]]
        t0 = div_tm[0]
        t1 = div_tm[-1]
        tds = ds.sel(time=slice(t0, t1))
        t0str = pd.to_datetime(t0).strftime('%Y-%m-%dT%H:%M')
        t1str = pd.to_datetime(t1).strftime('%Y-%m-%dT%H:%M')
        t0_savestr = pd.to_datetime(t0).strftime('%Y%m%dT%H%M')
        for pv, info in plt_vars.items():
            try:
                variable = tds[pv]
            except KeyError:
                continue
            # plot xsection
            if len(variable) > 1:
                if plot_xsection:
                    fig, ax = plt.subplots(figsize=(12, 8))
                    figttl_xsection = f'{glider} {info["ttl"].split(" (")[0]}\n{t0str} to {t1str}'
                    xargs = dict()
                    xargs['clabel'] = info['ttl']
                    xargs['title'] = figttl_xsection
                    xargs['date_fmt'] = '%m-%d\n%H:%M'
                    xargs['grid'] = True
                    pf.xsection(fig, ax, tds.time.values, tds.depth.values, variable.values, **xargs)

                    sfilename = f'{glider}_xsection_{pv}_{t0_savestr}.png'
                    sfile = os.path.join(sdir_xsection, pv, sfilename)
                    plt.savefig(sfile, dpi=300)
                    plt.close()

                if plot_profiles:
                    # plot profiles
                    colors = plt.cm.rainbow(np.linspace(0, 1, len(div_idxs)))
                    fig, ax = plt.subplots(figsize=(8, 10))
                    figttl_profile = f'{glider} {info["ttl"].split(" (")[0]}\n{t0str} to {t1str}'
                    pargs = dict()
                    pargs['xlabel'] = info['ttl']
                    pargs['title'] = figttl_profile
                    pargs['grid'] = True
                    pf.profile_yos(fig, ax, ds[pv], div_idxs, colors, **pargs)
                    sfilename = f'{glider}_profile_{pv}_{t0_savestr}.png'
                    if plot_qcvars:
                        if 'ph_total' in pv:
                            ph_qcvars = ['qartod_{}_gross_range_flag'.format(pv), 'qartod_{}_spike_flag'.format(pv)]
                            for pqc in ph_qcvars:
                                ax.plot()
                                tds[pv][np.where(tds[pqc] == 1)] = np.nan  # turn "good" values to nan
                                if np.sum(~np.isnan(tds[pv].values)) > 0:  # if there are any values not "good", plot them red
                                    xmask = ~np.isnan(tds[pv])
                                    ax.scatter(tds[pv][xmask], tds.depth.values[xmask], marker='x', color='k', s=60, label=pqc)
                                    ax.legend()
                                    if sfilename.startswith('qc-'):
                                        continue
                                    else:
                                        sfilename = f'qc-{sfilename}'

                    sfile = os.path.join(sdir_profile, pv, sfilename)
                    plt.savefig(sfile, dpi=300)
                    plt.close()


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210716T1804/delayed/ru30-20210716T1804-profile-sci-delayed_shifted.nc'
    number_yos = 6  # number of yos (down-up profile pairs) to plot in one image
    date_range = None   # [dt.datetime(2021, 3, 16), dt.datetime(2021, 3, 18)]  # None
    kwargs = {}
    # kwargs['plot_xsection'] = False
    # kwargs['plot_profiles'] = False
    # kwargs['plot_qcvars'] = True
    main(ncfile, number_yos, date_range, **kwargs)
