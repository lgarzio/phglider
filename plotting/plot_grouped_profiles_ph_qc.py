"""
Author: lgarzio on 4/5/2022
Last modified: lgarzio on 4/5/2022
Plot groups of profiles for corrected pH, values flagged by QC variables are highlighted.
"""

import numpy as np
import pandas as pd
import xarray as xr
import os
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})


def define_markers(qc_varname):
    markers = dict(climatology=dict(m='v', s=60, alpha=1),
                   hysteresis=dict(m='s', s=40, alpha=1),
                   flat_line=dict(m='^', s=60, alpha=1),
                   gross_range=dict(m='D', s=40, alpha=1),
                   rate_of_change=dict(m='X', s=80, alpha=1),
                   spike=dict(m='*', s=100, alpha=1),
                   summary=dict(m='o', s=100, alpha=.2)
                   )
    mkey = [key for key in markers.keys() if key in qc_varname][0]
    return markers[mkey]


def flatten(lst):
    return [item for sublist in lst for item in sublist]


def main(ncf, nprof):
    fname = ncf.split('/')[-1].split('.nc')[0]
    deploy = f'{fname.split("-")[0]}-{fname.split("-")[1]}'

    ds = xr.open_dataset(ncf)

    savedir = os.path.join(os.path.dirname(ncf), 'plots', f'profiles_group{nprof}')
    os.makedirs(savedir, exist_ok=True)

    profiletimes = np.unique(ds.profile_time.values)

    plot_sections = np.arange(0, len(profiletimes), nprof)
    plot_sections = np.append(plot_sections, len(profiletimes))

    varlist = ['ph_total_shifted_noqc']

    flag_defs = dict(not_evaluated=dict(value=2, color='cyan'),
                     suspect=dict(value=3, color='orange'),
                     fail=dict(value=4, color='red'))

    for ps_idx, ps in enumerate(plot_sections):
        if ps_idx > 0:
            if ps_idx == 1:
                ii = 0
            else:
                ii = plot_sections[ps_idx - 1] + 1
            ptimes = profiletimes[ii:ps]
            try:
                ptimes_idx = np.where(np.logical_and(ds.profile_time >= ptimes[0], ds.profile_time <= ptimes[-1]))[0]
            except IndexError:
                continue
            time0 = np.nanmin(ds.time.values[ptimes_idx])
            time1 = np.nanmax(ds.time.values[ptimes_idx])
            dss = ds.sel(time=slice(time0, time1))
            t0str = pd.to_datetime(np.nanmin(dss.profile_time.values)).strftime('%Y-%m-%dT%H:%M')
            t1str = pd.to_datetime(np.nanmax(dss.profile_time.values)).strftime('%Y-%m-%dT%H:%M')
            t0save = pd.to_datetime(np.nanmin(dss.profile_time.values)).strftime('%Y%m%dT%H%M')
            t1save = pd.to_datetime(np.nanmax(dss.profile_time.values)).strftime('%Y%m%dT%H%M')
            for cv in varlist:
                save_filename = f'{cv}_qc_{t0save}-{t1save}.png'

                data = dss[cv]

                # in some cases, ERDDAP doesn't set the metadata/fill values, so get rid of any possible fill values
                data[data > 10000] = np.nan

                pressure_interp = dss.pressure_interpolated
                fig, ax = plt.subplots(figsize=(8, 10))

                # iterate through each profile and plot the profile lines
                for pt in ptimes:
                    pt_idx = np.where(dss.profile_time.values == pt)[0]
                    non_nans = np.where(np.invert(np.isnan(data[pt_idx])))[0]
                    ax.plot(data[pt_idx][non_nans], pressure_interp[pt_idx][non_nans], color='gray')  # plot lines

                # add points
                ax.scatter(data, pressure_interp, color='gray', s=20, zorder=5)

                # find the qc variables
                qc_vars = ['ph_total_shifted_qartod_spike_test', 'ph_total_shifted_qartod_gross_range_test']

                for qi, qv in enumerate(qc_vars):
                    try:
                        flag_vals = dss[qv].values
                    except KeyError:
                        continue
                    for fd, info in flag_defs.items():
                        qc_idx = np.where(flag_vals == info['value'])[0]
                        if len(qc_idx) > 0:
                            m_defs = define_markers(qv)
                            ax.scatter(data[qc_idx], pressure_interp[qc_idx], color=info['color'], s=m_defs['s'],
                                       marker=m_defs['m'], edgecolor='k', alpha=m_defs['alpha'],
                                       label=f'{qv}-{fd}', zorder=10)

                # add legend if necessary
                handles, labels = plt.gca().get_legend_handles_labels()
                by_label = dict(zip(labels, handles))
                if len(handles) > 0:
                    ax.legend(by_label.values(), by_label.keys(), loc='best')

                ax.invert_yaxis()
                ax.set_ylabel('Pressure (dbar)')

                ax.ticklabel_format(useOffset=False)  # don't use scientific notation for ticks

                try:
                    units = data.units
                except AttributeError:
                    units = 'no_attributes'

                ax.set_xlabel(f'{data.name} ({units})')
                ttl = f'{deploy} {t0str} to {t1str}'
                ax.set_title(ttl)

                sfile = os.path.join(savedir, save_filename)
                plt.savefig(sfile, dpi=300)
                plt.close()


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210226T1647/delayed/ru30-20210226T1647-profile-sci-delayed_qc.nc'
    profile_group_n = 10
    main(ncfile, profile_group_n)
