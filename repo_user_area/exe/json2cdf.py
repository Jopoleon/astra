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

def json_concat(expequ):
    
    loc = '%s/ncdf_out/%s' %(awd, expequ)

    cdf_out = '%s-js.CDF' %loc

    j_json = 1
    f_json = '%s-%d.json' %(loc, j_json)

# Don't store anything if 1st json is older than tmp/astra.nml
    f_log2 = '%s/tmp/%s.nml' %(awd, expequ)

    if os.stat(f_json).st_mtime < os.stat(f_log2).st_mtime:
        logger.warning('json files older then %s' %f_log2)
        logger.warning('No 2D NetCDF file written')
        return

    ds = {}
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
            var_d = json.load(fjson)
        nx = len(var_d['NEX']['data'])
#        n_eq = len(var_d['RHO_SURF']['data'])
#        n_th = len(var_d['THETA']['data'])

        for key, val in var_d.items():
            ds[key] = []
            dat = val['data']
            if key not in ('XRHO', 'THETA', 'RHO_SURF'):
                if j_json == 1:
                    if type(dat) == type([]): # list
                        if len(dat) > 0:
                            ds[key] = dat
                        else: # array, all zeros
                            ds[key] = nx*[0.]
                    else: # scalar
                        ds[key] = [dat]
#                    logger.debug('#1 %s', key)
                else:
                    if type(dat) == type([]): # list
                        if len(dat) > 0:
                            ds[key] += dat
                        else: # array, all zeros
                            ds[key] += nx*[0.]
                    else: # scalar
                        ds[key].append(dat)
                        
        f_json_prev = f_json
        j_json += 1
    nt = j_json
    logger.info('nt = %d, nrho = %d, n_json = %d' %(nt, nx, j_json-1))
    
    dtyp = '>f8' #np.float64
    for key, val in ds.items():
#        print(key, val)
        if len(val) == nt:
            val = np.array(val, dtype=dtyp)
            var_d[key]['dimensions'] = ['TIME']
        if len(val) == nt*nx:
            val = val.reshape((nt, nx))
            var_d[key]['dimensions'] = ['TIME', 'XRHO']

    f = netcdf_file(cdf_out, 'w', mmap=False)

    f.createDimension('TIME', nt)
    f.createDimension('XRHO', nx)
#    f.createDimension('RHO_SURF', n_eq)
#    f.createDimension('THETA', n_th)

    rho = f.createVariable('XRHO', dtyp, ('XRHO', ))
    rho.data = var_d['XRHO']['data'].astype(dtyp)
    rho.units = var_d['XRHO']['units']
    rho.long_name = var_d['XRHO']['long_name']

    time = f.createVariable('TIME', dtyp, ('TIME', ))
    time.data = ds['TIME'].astype(dtyp)
    time.units = 's'
    time.long_name = 'Time'

#    rho_surf = f.createVariable('RHO_SURF', dtyp, ('RHO_SURF', ))
#    rho_surf.data = var_d['RHO_SURF']['data'].astype(dtyp)
#    rho_surf.units = '-'
#    rho_surf.long_name = var_d['RHO_SURF']['long_name']

#    theta = f.createVariable('THETA', dtyp, ('THETA', ))
#    theta.data = var_d['THETA']['data'].astype(dtyp)
#    theta.units = 'rad'
#    theta.long_name = var_d['THETA']['long_name']

    for key, val in ds.items():
        if key != 'TIME':
            if 'TIME' in var_d[key]['dimensions']:
                dims = var_d[key]['dimensions']
            else:
                dims = ('TIME', ) + var_d[key]['dimensions']
            tmp = f.createVariable(key, dtyp, dims)
            tmp[:] = val
            tmp.units = var_d[key]['units']
            tmp.long_name = var_d[key]['long_name']

    f.close()
    logger.info('Stored %s' %cdf_out)


if __name__ == '__main__':

    exp = 'aug34954'
    equ = 'test'
    json_concat(exp+equ)
