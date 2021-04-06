#!/usr/bin/env python

"""
Author: Lori Garzio on 4/3/2021
Last modified: 4/6/2021
Plot raw and time-shifted glider pH and oxygen data at user-specified intervals to evaluate the applied time shift.
"""

import os
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams.update({'font.size': 16})
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def plot_scatter(figure, axis, x, y, color, plt_ttl, cmin=None, cmax=None):
    if cmin:
        xc = axis.scatter(x, y, c=color, cmap=cm.jet, s=10, edgecolor='None', vmin=cmin, vmax=cmax)
    else:
        xc = axis.scatter(x, y, c=color, cmap=cm.jet, s=10, edgecolor='None')

    # format colorbar
    divider = make_axes_locatable(axis)
    cax = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
    figure.add_axes(cax)
    cb = plt.colorbar(xc, cax=cax)

    # add y-axis label
    axis.set_ylabel('Depth (m)')

    # format x-axis
    xfmt = mdates.DateFormatter('%m-%d\n%H:%M')
    axis.xaxis.set_major_formatter(xfmt)
    axis.xaxis.set_tick_params(labelsize=15)

    axis.set_title(plt_ttl, fontsize=18)


def main(fname, plt_interval):
    ds = xr.open_dataset(fname)
    plotting_vars = [['ph_total', 'ph_total_shifted'],
                     ['oxygen_concentration', 'oxygen_concentration_shifted']]

    for pv in plotting_vars:
        for p in pv:
            if 'ph' in p:
                da = ds[p]
                da = da.where(da < 14, np.nan)  # convert pH values > 14 to nan
                ds[p] = da

    sdir = os.path.join(os.path.dirname(fname), 'timeshift_figs')
    os.makedirs(sdir, exist_ok=True)

    daterange = pd.date_range(ds.time.values[0], ds.time.values[-1], freq=plt_interval)
    for dt_idx, dr in enumerate(daterange):
        if dt_idx > 0:
            ds_sel = ds.sel(time=slice(daterange[dt_idx - 1], dr))
            for plt_vars in plotting_vars:
                fig, axs = plt.subplots(2, figsize=(12, 12), sharex=True, sharey=True)
                for i, ax in enumerate(fig.axes):
                    if i == 0:
                        plot_scatter(fig, ax, ds_sel.time.values, ds_sel.depth.values, ds_sel[plt_vars[0]].values,
                                     plt_vars[0])
                    else:
                        plot_scatter(fig, ax, ds_sel.time.values, ds_sel.depth.values, ds_sel[plt_vars[1]].values,
                                     plt_vars[1])

                ax.invert_yaxis()

                sname = os.path.join(sdir, '{}_{}.png'.format(plt_vars[0], daterange[dt_idx - 1].strftime('%Y%m%dT%H%M')))
                plt.savefig(sname, dpi=300)
                plt.close()


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed_shifted.nc'
    plot_interval = '3H'  # 3 hours
    main(ncfile, plot_interval)
