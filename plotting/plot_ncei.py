#!/usr/bin/env python

"""
Author: Lori Garzio on 8/10/2021
Last modified: 8/10/2021
Plot data variables from the dataset that will be sent to the IOOS Glider DAC
"""

import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import functions.common as cf
import functions.plotting as pf
plt.rcParams.update({'font.size': 13})


def main(fname):
    ds = xr.open_dataset(fname)
    save_dir = os.path.join(os.path.dirname(fname), 'ncei_plots')
    os.makedirs(save_dir, exist_ok=True)
    deploy = ds.attrs['title']

    t0str = pd.to_datetime(np.nanmin(ds.time)).strftime('%Y-%m-%dT%H:%M')
    t1str = pd.to_datetime(np.nanmax(ds.time)).strftime('%Y-%m-%dT%H:%M')

    plt_vars = cf.plot_vars_ncei()
    for pv, info in plt_vars.items():
        try:
            variable = ds[pv]
        except KeyError:
            continue

        # plot xsection
        if len(variable) > 1:
            fig, ax = plt.subplots(figsize=(12, 6))
            plt.subplots_adjust(left=0.1)
            figttl_xsection = f'{deploy} {info["ttl"].split(" (")[0]}\n{t0str} to {t1str}'
            xargs = dict()
            xargs['clabel'] = info['ttl']
            xargs['title'] = figttl_xsection
            xargs['date_fmt'] = '%m-%d'
            xargs['grid'] = True
            xargs['ylabel'] = ds.pressure_interpolated.units
            xargs['cmap'] = info['cmap']
            pf.xsection(fig, ax, ds.time.values, ds.pressure_interpolated.values, variable.values, **xargs)

            sfilename = f'{deploy}_xsection_{pv}.png'
            sfile = os.path.join(save_dir, sfilename)
            plt.savefig(sfile, dpi=300)
            plt.close()


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ncei/ru30-20210226T1647-delayed.nc'
    main(ncfile)
