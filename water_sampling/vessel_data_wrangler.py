#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2021
Last modified: 10/2/2024
Grab vessel-based salinity and TA data from CODAP-NA and additional ECOMON and ECOA datasets within a defined
bounding box spanning only the New York Bight and limited to 200m depth (max glider depth). QC flags were applied (when
provided). Export as NetCDF and csv to use for TA-salinity regressions for estimating glider-based Total Alkalinity.
CODAP-NA dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
Additional cruise datasets were downloaded from the NCEI OCADs data portal (https://www.ncei.noaa.gov/products/ocean-carbon-acidification-data-system)
"""

import os
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
from collections import OrderedDict
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console


def main(lon_bounds, lat_bounds, codap_file, ecomon_files, ecomon_files2, ecoa_files, ecoa_files2, savedir):
    # initialize dictionary to append salinity and TA data from cruises (CODAP and additional datasets)
    data = {
        "coords": {
            "time": {"dims": "time", "data": np.array([], dtype='datetime64[ns]')}
        },
        "attrs": {
            "comment": "Synthesis of salinity and TA data from vessel-based measurements "
                       "that were spatially limited to the U.S. Northeast Shelf.",
            "data_sources": "CODAP-NA dataset (https://essd.copernicus.org/articles/13/2777/2021/) and subsequent "
                            "ECOMON and ECOA cruises not included in the CODAP-NA dataset. Additional cruise data were "
                            "downloaded from NCEI's Ocean Carbon and Acidification Data Portal "
                            "(https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system-portal/)"
        },
        "dims": "time",
        "data_vars": {
            "data_source": {
                "dims": "time",
                "data": np.array([], dtype='<U32'),
                "attrs": {
                    "units": "1",
                    "comment": "Source of data"
                }
            },
            "cruise": {
                "dims": "time",
                "data": np.array([], dtype='<U32'),
                "attrs": {
                    "units": "1",
                    "comment": "Cruise ID from original data source"
                }
            },
            "obs_type": {
                "dims": "time",
                "data": np.array([], dtype='<U32'),
                "attrs": {
                    "units": "1",
                    "comment": "Observation type"
                }
            },
            "accession": {
                "dims": "time",
                "data": np.array([], dtype='int32'),
                "attrs": {
                    "units": "1",
                    "comment": "NCEI Accession number"
                }
            },
            "lat": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "degrees_north",
                    "long_name": "Latitude",
                    "comment": "Latitude from original data source"
                }
            },
            "lon": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "degrees_east",
                    "long_name": "Longitude",
                    "comment": "Longitude from original data source"
                }
            },
            "depth": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "m",
                    "long_name": "Depth",
                    "description": "Depth at which sample was collected",
                    "comment": "From the original data source"
                }
            },
            "salinity": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "Salinity"
                }
            },
            "total_alkalinity": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "umol/kg",
                    "long_name": "Total Alkalinity"
                }
            },
        },
    }

    # get CODAP data
    ds = xr.open_dataset(codap_file)
    idx = []
    codap_vars = dict(Day_UTC=np.array([]),
                      Month_UTC=np.array([]),
                      Year_UTC=np.array([]),
                      Cruise_ID=np.array([]),
                      Accession=np.array([]),
                      Observation_type=np.array([]),
                      Profile_number=np.array([]),
                      Latitude=np.array([]),
                      Longitude=np.array([]),
                      CTDPRES=np.array([]),
                      Depth=np.array([]),
                      Depth_bottom=np.array([]),
                      CTDTEMP_ITS90=np.array([]),
                      CTDTEMP_flag=np.array([]),
                      recommended_Salinity_PSS78=np.array([]),
                      recommended_Salinity_flag=np.array([]),
                      pH_TS_insitu_calculated=np.array([]),
                      pH_TS_insitu_measured=np.array([]),
                      pH_flag=np.array([]),
                      TALK=np.array([]),
                      TALK_flag=np.array([]),
                      Aragonite=np.array([]))

    # make sure the data are within the defined extent
    for i, lon in enumerate(ds.Longitude.values):
        if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(lon, ds.Latitude.values[i])):
            idx.append(i)
            for key in codap_vars.keys():
                if key in ['Cruise_ID', 'Observation_type']:
                    cid = ds[key].values[:, i]
                    cid = [x.decode('UTF-8') for x in cid]
                    codap_vars[key] = np.append(codap_vars[key], ''.join(cid).strip())
                else:
                    codap_vars[key] = np.append(codap_vars[key], ds[key].values[i])

    df = pd.DataFrame(codap_vars)

    # generate timestamp
    df['year'] = df['Year_UTC'].apply(int)
    df['month'] = df['Month_UTC'].apply(int)
    df['day'] = df['Day_UTC'].apply(int)
    df['time'] = pd.to_datetime(df[['year', 'month', 'day']])

    # remove questionable (3) and missing (9) TALK
    df = df[df.TALK_flag != 3]
    df = df[df.TALK_flag != 9]

    # remove questionable (3) and missing (9) salinity
    df = df[df.recommended_Salinity_flag != 3]
    df = df[df.recommended_Salinity_flag != 9]

    # remove data >200m depth
    df = df[df.Depth <= 200]

    # add data to dictionary
    data['coords']['time']['data'] = np.array(df.time)
    data['data_vars']['data_source']['data'] = np.repeat('CODAP_NA_v2021', len(df.time))
    data['data_vars']['cruise']['data'] = np.array(df.Cruise_ID, dtype=data['data_vars']['cruise']['data'].dtype)
    data['data_vars']['obs_type']['data'] = np.array(df.Observation_type, dtype=data['data_vars']['obs_type']['data'].dtype)
    data['data_vars']['accession']['data'] = np.array(df.Accession, dtype=data['data_vars']['accession']['data'].dtype)
    data['data_vars']['depth']['data'] = np.array(df.Depth, dtype=data['data_vars']['depth']['data'].dtype)
    data['data_vars']['lat']['data'] = np.array(df.Latitude, dtype=data['data_vars']['lat']['data'].dtype)
    data['data_vars']['lon']['data'] = np.array(df.Longitude, dtype=data['data_vars']['lon']['data'].dtype)
    data['data_vars']['salinity']['data'] = np.array(df.recommended_Salinity_PSS78, dtype=data['data_vars']['salinity']['data'].dtype)
    data['data_vars']['total_alkalinity']['data'] = np.array(df.TALK, dtype=data['data_vars']['total_alkalinity']['data'].dtype)

    # additional datasets that aren't included in CODAP
    accession_mapping = {'HB1902': 209045,
                         'GU1902': 209156,
                         'GU1905': 210238,
                         'GU2102': 248269,
                         'PC2104': 249432,
                         'PC2106': 249517,
                         'PC2205': 283758,
                         'HB2204': 276023,
                         'ECOA3': 283329}
    accession_underway_mapping = {'ECOA-1': 157389,
                                  'ECOA-2': 215462,
                                  'ECOA3': 295751,
                                  'HB2204': 288996
                                  }
    for ef in ecomon_files:
        df = pd.read_csv(ef)

        # remove questionable (3) bad (4) and missing (5 or 9 or -999) TA
        df = df[df.TA_Flag != 3]
        df = df[df.TA_Flag != 4]
        df = df[df.TA_Flag != 5]
        df = df[df.TA_Flag != 9]
        df = df[df.TA_Flag != -999]

        # remove missing salinity
        try:
            df = df[df.CTDSAL_PSS78 != -999]
        except AttributeError:
            df = df[df['CTDSAL (PSS-78)'] != -999]

        # remove data >200m depth
        try:
            df = df[df.Depth_meters <= 200]
        except AttributeError:
            df = df[df['Depth_sampling (M)'] <= 200]

        # format date
        try:
            df['year'] = df['Year_UTC'].apply(int)
            df['month'] = df['Month_UTC'].apply(int)
            df['day'] = df['Day_UTC'].apply(int)
            df['time'] = pd.to_datetime(df[['year', 'month', 'day']])
        except KeyError:
            df['time'] = pd.to_datetime(df['Date_UTC'])

        # make sure the data are within the defined extent
        df['in_region'] = ''
        for i, row in df.iterrows():
            try:
                lon = row.Longitude_Dec_Deg
                lat = row.Latitude_Dec_Deg
            except AttributeError:
                lon = row.Longitude
                lat = row.Latitude
            if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(lon, lat)):
                df.loc[i, 'in_region'] = 'yes'
            else:
                df.loc[i, 'in_region'] = 'no'

        # drop data if it's not in the region specified
        df = df[df.in_region == 'yes']

        # add data to dictionary
        if len(df) > 0:
            data['coords']['time']['data'] = np.append(data['coords']['time']['data'], np.array(df.time))
            data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'], np.repeat('ECOMON-NCEI', len(df.time)))
            data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'], np.array(df.Cruise_ID))
            data['data_vars']['obs_type']['data'] = np.append(data['data_vars']['obs_type']['data'], np.array(df.Observation_Type))
            cruise_accession = accession_mapping[np.unique(df.Cruise_ID)[0]]
            data['data_vars']['accession']['data'] = np.append(data['data_vars']['accession']['data'], np.repeat(cruise_accession, len(df.time)))
            try:
                data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'], np.array(df.Depth_meters))
                data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], np.array(df.Latitude_Dec_Deg))
                data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], np.array(df.Longitude_Dec_Deg))
                data['data_vars']['salinity']['data'] = np.append(data['data_vars']['salinity']['data'], np.array(df.CTDSAL_PSS78))
                data['data_vars']['total_alkalinity']['data'] = np.append(data['data_vars']['total_alkalinity']['data'], np.array(df['TA_umol/kg']))
            except AttributeError:  # column names changed
                data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'], np.array(df['Depth_sampling (M)']))
                data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], np.array(df.Latitude))
                data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], np.array(df.Longitude))
                data['data_vars']['salinity']['data'] = np.append(data['data_vars']['salinity']['data'], np.array(df['CTDSAL (PSS-78)']))
                data['data_vars']['total_alkalinity']['data'] = np.append(data['data_vars']['total_alkalinity']['data'], np.array(df['TA (umol/kg)']))

    # ECOMON underway flow through data
    for ef in ecomon_files2:
        df = pd.read_csv(ef)
        df['time'] = df['DATE_UTC (yyyymmdd)'].map(lambda t: pd.to_datetime(str(t)))

        # make sure the data are within the defined extent
        df['in_region'] = ''
        for i, row in df.iterrows():
            try:
                lon = row.LONGITUDE
                lat = row.LATITUDE
            except AttributeError:
                lon = row.Longitude
                lat = row.Latitude
            if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(lon, lat)):
                df.loc[i, 'in_region'] = 'yes'
            else:
                df.loc[i, 'in_region'] = 'no'

        # drop data if it's not in the region specified
        df = df[df.in_region == 'yes']

        # add data to dictionary
        if len(df) > 0:
            data['coords']['time']['data'] = np.append(data['coords']['time']['data'], np.array(df.time))
            data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'],
                                                                 np.repeat('ECOMON-NCEI', len(df.time)))
            data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'],
                                                            np.array(df.CRUISE_ID))
            data['data_vars']['obs_type']['data'] = np.append(data['data_vars']['obs_type']['data'],
                                                              np.repeat('underway', len(df.time)))
            cruise_accession = accession_underway_mapping[np.unique(df.CRUISE_ID)[0]]
            data['data_vars']['accession']['data'] = np.append(data['data_vars']['accession']['data'],
                                                               np.repeat(cruise_accession, len(df.time)))

            data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'],
                                                           np.repeat(5, len(df.time)))
            data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'],
                                                         np.array(df.LATITUDE))
            data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'],
                                                         np.array(df.LONGITUDE))
            data['data_vars']['salinity']['data'] = np.append(data['data_vars']['salinity']['data'],
                                                              np.array(df.SSS))
            data['data_vars']['total_alkalinity']['data'] = np.append(data['data_vars']['total_alkalinity']['data'],
                                                                      np.array(df.TA))

    # additional ECOA datasets that aren't included in CODAP
    for ef in ecoa_files:
        df = pd.read_csv(ef)

        # remove questionable (3) bad (4) and missing (9) TA
        df = df[df.TA_flag != 3]
        df = df[df.TA_flag != 4]
        df = df[df.TA_flag != -999]

        # remove missing salinity
        df = df[df.CTDSAL_PSS78_flag != -999]

        # remove data >200m depth
        df = df[df.Depth <= 200]

        # format date
        df['year'] = df['Year_UTC'].apply(int)
        df['month'] = df['Month_UTC'].apply(int)
        df['day'] = df['Day_UTC'].apply(int)
        df['time'] = pd.to_datetime(df[['year', 'month', 'day']])

        # make sure the data are within the defined extent
        df['in_region'] = ''
        for i, row in df.iterrows():
            if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(row.Longitude, row.Latitude)):
                df.loc[i, 'in_region'] = 'yes'
            else:
                df.loc[i, 'in_region'] = 'no'

        # drop data if it's not in the region specified
        df = df[df.in_region == 'yes']

        # add data to dictionary
        data['coords']['time']['data'] = np.append(data['coords']['time']['data'], np.array(df.time))
        data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'], np.repeat('ECOA-NCEI', len(df.time)))
        data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'], np.array(df.Cruise_ID))
        data['data_vars']['obs_type']['data'] = np.append(data['data_vars']['obs_type']['data'], np.repeat('Niskin', len(df.time)))
        cruise_accession = accession_mapping[np.unique(df.Cruise_ID)[0]]
        data['data_vars']['accession']['data'] = np.append(data['data_vars']['accession']['data'], np.repeat(cruise_accession, len(df.time)))
        data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'], np.array(df.Depth))
        data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], np.array(df.Latitude))
        data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], np.array(df.Longitude))
        data['data_vars']['salinity']['data'] = np.append(data['data_vars']['salinity']['data'], np.array(df.CTDSAL_PSS78))
        data['data_vars']['total_alkalinity']['data'] = np.append(data['data_vars']['total_alkalinity']['data'], np.array(df['TA']))

    # ECOA underway flow through data
    for ef in ecoa_files2:
        df = pd.read_csv(ef)

        # remove questionable (3) bad (4) and missing (9) TA
        try:
            df = df[df.TA_flag != 3]
            df = df[df.TA_flag != 4]
            df = df[df.TA_flag != -999]
            df = df[df.TA_flag != 9]
        except AttributeError:
            print('no TA flags')

        # remove missing salinity
        df = df[df.SSS != -999]

        # format date
        try:
            df['time'] = pd.to_datetime(df[['Year', 'Month', 'Day']])
        except KeyError:
            df['time'] = df['Date_UTC'].map(lambda t: pd.to_datetime(str(t)))

        # make sure the data are within the defined extent
        df['in_region'] = ''
        for i, row in df.iterrows():
            if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(row.Longitude, row.Latitude)):
                df.loc[i, 'in_region'] = 'yes'
            else:
                df.loc[i, 'in_region'] = 'no'

        # drop data if it's not in the region specified
        df = df[df.in_region == 'yes']

        # add data to dictionary
        data['coords']['time']['data'] = np.append(data['coords']['time']['data'], np.array(df.time))
        data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'],
                                                             np.repeat('ECOA-NCEI', len(df.time)))
        data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'],
                                                        np.array(df.CRUISE_ID))
        data['data_vars']['obs_type']['data'] = np.append(data['data_vars']['obs_type']['data'],
                                                          np.repeat('underway', len(df.time)))
        cruise_accession = accession_underway_mapping[np.unique(df.CRUISE_ID)[0]]
        data['data_vars']['accession']['data'] = np.append(data['data_vars']['accession']['data'],
                                                           np.repeat(cruise_accession, len(df.time)))
        data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'], np.repeat(5, len(df.time)))
        data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], np.array(df.Latitude))
        data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], np.array(df.Longitude))
        data['data_vars']['salinity']['data'] = np.append(data['data_vars']['salinity']['data'],
                                                          np.array(df.SSS))
        data['data_vars']['total_alkalinity']['data'] = np.append(data['data_vars']['total_alkalinity']['data'],
                                                                  np.array(df.TA_umol_kg))

    # save as netcdf
    outds = xr.Dataset.from_dict(data)

    # add season
    outds['season'] = outds['time.season']

    # add created time to global attrs
    datetime_format = '%Y-%m-%dT%H:%M:%SZ'
    savedate = dt.datetime.utcnow().strftime('%Y%m%d')
    created = dt.datetime.utcnow().strftime(datetime_format)  # creation time Timestamp
    time_start = pd.to_datetime(np.nanmin(outds.time.values)).strftime(datetime_format)
    time_end = pd.to_datetime(np.nanmax(outds.time.values)).strftime(datetime_format)

    global_attributes = OrderedDict([
        ('date_created', created),
        ('date_modified', created),
        ('time_coverage_start', time_start),
        ('time_coverage_end', time_end),
        ('creator_email', 'lgarzio@marine.rutgers.edu'),
        ('creator_name', 'Lori Garzio'),
        ('creator_url', 'rucool.marine.rutgers.edu'),
        ('institution', 'Rutgers University'),
        ('contributor_name', 'Grace Saba,Lori Garzio'),
        ('contributor_role', 'Principal Investigator,Data Management')
    ])

    global_attributes.update(outds.attrs)

    outds = outds.assign_attrs(global_attributes)
    outds = outds.sortby(outds.time)

    # Add compression to all variables
    encoding = {}
    for k in outds.data_vars:
        encoding[k] = {'zlib': True, 'complevel': 1}

    encoding['time'] = dict(units='seconds since 1970-01-01 00:00:00', calendar='gregorian', zlib=False,
                            _FillValue=False, dtype=np.double)

    save_file = os.path.join(savedir, f'vessel_based_TA_salinity_NYB_{savedate}.nc')
    outds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')

    # save as csv
    save_file = os.path.join(savedir, f'vessel_based_TA_salinity_NYB_{savedate}.csv')
    df = outds.to_pandas()
    df.reset_index(inplace=True)
    df.to_csv(save_file, index=False)


if __name__ == '__main__':
    lons = [-75, -73.25, -72.5, -71.5, -71.65, -74]  # NYB
    lats = [38.75, 37.5, 37.5, 39.25, 41.2, 40.5]
    codap = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/CODAP_NA_v2021.nc'
    ecomon = ['/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/EcoMon/2019/Accession_0209045-discrete-profiles/33HH20190522-HB1902_data.csv',
              '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/EcoMon/2019/Accession_0209156-discrete-profiles/33GG20190815-GU1902_data.csv',
              '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/EcoMon/2019/Accession_0210238-discrete-profiles/33GG20191015-GU1905_data.csv',
              '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/EcoMon/2021/Accession_0248269-discrete-profiles/33GG20210514-GU2102_data.csv',
              '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/EcoMon/2021/Accession_0249432-discrete-profiles/334B20210805-PC2104_Data.csv',
              '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/EcoMon/2021/Accession_0249517-discrete-profiles/334B20211015-PC2106_Data.csv',
              '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/EcoMon/2022/Accession_0276023-discrete-profiles/33HH20220531_HB2204_Data.csv',
              '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/EcoMon/2022/Accession_0283758-discrete-profiles/334B20221101_PC2205_Data-mod.csv']
    ecomon_underway = ['/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/EcoMon/2022/Accession_0288996-surface-underway/33HH20220531-mod.csv']
    ecoa = ['/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/ECOA/ECOA-3/Accession_0283329-discrete/ECOA_3_CTD_MasterDataSheet_09_26_2023_Accession_0283329-mod.csv']
    ecoa_underway = ['/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/ECOA/ECOA-1/Accession_0157389-underway/Discrete_Underway_Data_12082016_Accession_0157389-mod.csv',
                     '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/ECOA/ECOA-2/Accession_0215462-surface-underway/ECOA_2_2018_Data_Accession_0215462-mod.csv',
                     '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data/ECOA/ECOA-3/Accession_0295751-underway/ECOA_3_Underway_MasterDataSheet_NCEI_20240731._Accession_0295751-mod.csv']
    save_directory = '/Users/garzio/Documents/rucool/Saba/TA_salinity_regressions/vessel-based-data'
    main(lons, lats, codap, ecomon, ecomon_underway, ecoa, ecoa_underway, save_directory)
