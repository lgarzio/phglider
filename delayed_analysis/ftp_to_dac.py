#!/usr/bin/env python

"""
Author: John Kerfoot
Modified by Lori Garzio 8/11/2022
"""

import logging
import sys
import argparse
import os
import glob
import yaml
import ftplib


def main(args):
#def main(dataset_id, ncdir, config_file):
    """FTP a list of NetCDF files to the DAC ftp server for the specified dataset_id"""

    dataset_id = args.dataset_id
    ncdir = args.ncdir
    config_file = args.config_file

    log_level = getattr(logging, args.loglevel.upper())
    #log_level = getattr(logging, 'INFO')
    log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
    logging.basicConfig(format=log_format, level=log_level)

    nc_files = sorted(glob.glob(os.path.join(ncdir, '*.nc')))

    if config_file:
        if not os.path.isfile(config_file):
            logging.error('Invalid configuration file specified: {:}'.format(config_file))
            return 1
        with open(config_file, 'r') as fid:
            cfg_params = yaml.safe_load(fid)
            user = cfg_params.get('username', None)
            pw = cfg_params.get('password', None)
            url = cfg_params.get('url', None)

    if not url:
        logging.error('No FTP site specified')
        return 1
    if not user:
        logging.error('No username specified')
        return 1
    if not pw:
        logging.error('No password specified')
        return 1

    ftp = ftplib.FTP(url, user, pw, timeout=300)

    with ftp:
        # Make sure the dataset directory exists
        try:
            resp = ftp.cwd(dataset_id)
        except ftplib.all_errors as e:
            logging.error('Dataset {:} does not exist on the remote server'.format(dataset_id))
            return 1

        # get the list of files already in that directory
        dac_nc_list = ftp.nlst('*.nc')
        logging.info(f'Found {len(dac_nc_list)} files already in the DAC for {dataset_id}')

        # compare the list of files in the DAC to list of local nc files
        upload_files = list(set([os.path.basename(x) for x in nc_files]) - set(dac_nc_list))

        if len(upload_files) > 0:
            # Upload each file that doesn't already exist in the DAC
            logging.info('Uploading {:} files for {:}'.format(len(upload_files), dataset_id))
            file_completed = []
            for nc_filename in sorted(upload_files):
                filepath = os.path.join(ncdir, nc_filename)
                try:
                    with open(filepath, 'br') as fid:
                        resp = ftp.storbinary('STOR {:s}'.format(nc_filename), fid)
                        file_completed.append(filepath)
                except ftplib.all_errors as e:
                    logging.error(e)
        else:
            logging.info('No new files to upload for {:}'.format(dataset_id))


if __name__ == '__main__':
    # dsid = 'ru30-20220906T1523-delayed'
    # nc_directory = '/Users/garzio/Documents/rucool/Saba/gliderdata/2022/ru30-20220906T1523/delayed/ngdac'
    # c = '/Users/garzio/Documents/rucool/gliders/data_archiving/example_code/dacftp.yml'
    # main(dsid, nc_directory, c)
    arg_parser = argparse.ArgumentParser(description=main.__doc__,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('dataset_id',
                            help='Registered DAC deployment dataset id',
                            type=str)

    arg_parser.add_argument('ncdir',
                            help='Filepath to directory containing NetCDF files to upload',
                            type=str)

    arg_parser.add_argument('-c', '--config_file',
                            help='YAML configuration file specifying FTP url and credentials')

    arg_parser.add_argument('-l', '--loglevel',
                            help='Verbosity level',
                            type=str,
                            choices=['debug', 'info', 'warning', 'error', 'critical'],
                            default='info')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
