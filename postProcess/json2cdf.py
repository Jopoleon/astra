#!/usr/bin/env python

import os, logging, argparse
from scipy.io import netcdf_file
import numpy as np
import json

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')

logger = logging.getLogger('json2cdf')
if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)

logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.INFO)

awd = os.getenv('AWD')
if awd is None:
    awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

dtyp = '>f8' #np.float64


def json_concat(expequ):

    loc = '%s/ncdf_out/%s' %(awd, expequ)

    cdf_out = '%s.CDF' %loc

    j_json = 1
    f_json = '%s-%d.json' %(loc, j_json)

# Don't store a new CDF if the current one is newer than 1st json
    if os.path.isfile(cdf_out):
        if os.stat(f_json).st_mtime < os.stat(cdf_out).st_mtime:
            logger.warning('json files older then %s' %cdf_out)
            logger.warning('No 2D NetCDF file written')
            return

    f_meta = '%s/astra_variables.json' %awd
    with open(f_meta, 'r') as fmeta:
        meta_d = json.load(fmeta)

    logger.info('Started json2cdf')
    ds_astra = {}
    ds_equil = {}
    while True:
        f_json = '%s-%d.json' %(loc, j_json)
        if not os.path.isfile(f_json):
            if j_json == 1:
                return
            else:
                break
        if j_json > 1:
            if os.stat(f_json).st_mtime < os.stat(f_json_prev).st_mtime: #Newer
                break
        logger.debug(f_json)
        with open(f_json, 'r') as fjson:
            json_d = json.load(fjson)

        for lbl, dic in json_d.items():
            for key, dat in dic.items():
                if j_json == 1: # First data, metadata from astra_variables.json
                    ds_astra[key] = {}
                    meta = meta_d[lbl][key]
                    if 'units' in meta:
                        ds_astra[key]['units'] = meta['units']
                    else:
                        ds_astra[key]['units'] = ''
                    if 'desc' in meta:
                        ds_astra[key]['long_name'] = meta['desc']
                    else:
                        ds_astra[key]['long_name'] = ''
                    ds_astra[key]['data'] = [dat]
                else: # append data
                    ds_astra[key]['data'].append(dat)

        f_json_prev = f_json
        j_json += 1

    nt = j_json - 1
    nx   = len(json_d['profiles']['XRHO'])
    n_eq = len(json_d['equil_profiles']['rho_tor_norm'])
    n_th = len(json_d['equil_coord']['teta2d'])
    nR   = len(json_d['equil_rect']['r2d'])
    nZ   = len(json_d['equil_rect']['z2d'])
    logger.debug('nt=%d, nrho=%d, nr_eq=%d, nthe_eq=%d' %(nt, nx, n_eq, n_th))

    for key, val in ds_astra.items():
        val['data'] = np.array(val['data'], dtype=dtyp)
        darr = val['data']
        if darr.shape == (nt, ):
            val['dimensions'] = ['TIME']
        elif darr.shape == (nt, nx) and key[:6] != 'equil_':
            val['dimensions'] = ['TIME', 'XRHO']
        elif darr.shape == (nt, n_eq) and key[:6] == 'equil_':
            val['dimensions'] = ['TIME', 'RHO_SURF']
        elif darr.shape == (nt, n_eq, n_th):
            val['dimensions'] = ['TIME', 'RHO_SURF', 'THETA']
        elif darr.shape == (nt, nR, nZ):
            val['dimensions'] = ['TIME', 'R', 'Z']

    f = netcdf_file(cdf_out, 'w', mmap=False)

    f.createDimension('TIME', nt)
    f.createDimension('XRHO', nx)
    f.createDimension('RHO_SURF', n_eq)
    f.createDimension('THETA', n_th)
    f.createDimension('R', nR)
    f.createDimension('Z', nZ)

    rho = f.createVariable('XRHO', dtyp, ('XRHO', ))
    rho.data  = np.array(json_d['profiles']['XRHO'], dtype=dtyp)
    rho.units = ''
    rho.long_name = meta_d['profiles']['XRHO']['desc']

    time = f.createVariable('TIME', dtyp, ('TIME', ))
    time.data = ds_astra['TIME']['data'].astype(dtyp)
    time.units = 's'
    time.long_name = 'Time'

    rho_surf = f.createVariable('RHO_SURF', dtyp, ('RHO_SURF', ))
    rho_surf.data = np.array(ds_astra['rho_tor_norm']['data'], dtype=dtyp)
    rho_surf.units = '-'
    rho_surf.long_name = meta_d['equil_profiles']['rho_tor_norm']['desc']

    theta = f.createVariable('THETA', dtyp, ('THETA', ))
    theta.data = np.array(ds_astra['teta2d']['data'], dtype=dtyp)
    theta.units = 'rad'
    theta.long_name = meta_d['equil_coord']['teta2d']['desc']

    rgrid = f.createVariable('R', dtyp, ('R', ))
    rgrid.data = np.array(ds_astra['r2d']['data'], dtype=dtyp)
    rgrid.units = 'm'
    rgrid.long_name = meta_d['equil_rect']['r2d']['desc']

    zgrid = f.createVariable('Z', dtyp, ('Z', ))
    zgrid.data = np.array(ds_astra['z2d']['data'], dtype=dtyp)
    zgrid.units = 'm'
    zgrid.long_name = meta_d['equil_rect']['z2d']['desc']

    for key, val in ds_astra.items():
        if key not in ('TIME', 'XRHO', 'rho_tor_norm', 'teta2d', 'r2d', 'z2d'):
            tmp = f.createVariable(key, dtyp, val['dimensions'])
            tmp[:] = val['data']
            tmp.units = val['units']
            tmp.long_name = val['long_name']

    f.close()
    logger.info('Stored %s' %cdf_out)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Concatenate json files to 1 NetCDF along time')
    parser.add_argument('-e', '--expequ', help='<exp><equ>', required=True)
    args = parser.parse_args()
    json_concat(args.expequ)
