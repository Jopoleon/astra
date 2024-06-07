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
        n_eq = len(var_d['RHO_SURF']['data'])
        n_th = len(var_d['THETA']['data'])

        for key, val in var_d.items():
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
    nt = j_json - 1
    logger.info('nt = %d, nrho = %d, n_json = %d' %(nt, nx, j_json-1))
    
    dtyp = '>f8' #np.float64
    for key, val in ds.items():
        ds[key] = np.array(val, dtype=dtyp)
        if key in ('phi', 'pprime', 'ffprim'):
            ds[key] = ds[key].reshape((nt, n_eq))
            var_d[key]['dimensions'] = ['TIME', 'RHO_SURF']
        elif key in ('r2d', 'z2d', 'psi2d'):
            ds[key] = np.transpose(ds[key].reshape((nt, n_th, n_eq)), (0, 2, 1) )
            var_d[key]['dimensions'] = ['TIME', 'RHO_SURF', 'THETA']
        elif len(val) == nt:
            var_d[key]['dimensions'] = ['TIME']
        elif len(val) == nt*nx:
            ds[key] = ds[key].reshape((nt, nx))
            var_d[key]['dimensions'] = ['TIME', 'XRHO']

    f = netcdf_file(cdf_out, 'w', mmap=False)

    f.createDimension('TIME', nt)
    f.createDimension('XRHO', nx)
    f.createDimension('RHO_SURF', n_eq)
    f.createDimension('THETA', n_th)

    rho = f.createVariable('XRHO', dtyp, ('XRHO', ))
    rho.data  = np.array(var_d['XRHO']['data'], dtype=dtyp)
    rho.units = var_d['XRHO']['units']
    rho.long_name = var_d['XRHO']['long_name']

    time = f.createVariable('TIME', dtyp, ('TIME', ))
    time.data = ds['TIME'].astype(dtyp)
    time.units = 's'
    time.long_name = 'Time'

    rho_surf = f.createVariable('RHO_SURF', dtyp, ('RHO_SURF', ))
    rho_surf.data = np.array(var_d['RHO_SURF']['data'], dtype=dtyp)
    rho_surf.units = '-'
    rho_surf.long_name = var_d['RHO_SURF']['long_name']

    theta = f.createVariable('THETA', dtyp, ('THETA', ))
    theta.data = np.array(var_d['THETA']['data'], dtype=dtyp)
    theta.units = 'rad'
    theta.long_name = var_d['THETA']['long_name']

    for key, val in ds.items():
        if key != 'TIME':
            tmp = f.createVariable(key, dtyp, var_d[key]['dimensions'])
            tmp[:] = val
            tmp.units = var_d[key]['units']
            tmp.long_name = var_d[key]['long_name']

    f.close()
    logger.info('Stored %s' %cdf_out)


if __name__ == '__main__':

    import matplotlib.pylab as plt

    exp = 'aug34954'
    equ = 'test'
    json_concat(exp+equ)

    f_cdf = '/shares/departments/AUG/users/git/a8/ncdf_out/aug34954test-js.CDF'
    cv = netcdf_file(f_cdf, 'r', mmap=False).variables

    plt.figure(1)
    plt.plot(cv['XRHO'].data, cv['TI'][-1, :])
    plt.figure(2)
    plt.plot(cv['r2d'][-1, -1, :], cv['z2d'][-1, -1, :])
    plt.figure(3)
    plt.plot(cv['TIME'].data, cv['IPL'].data)
    plt.show()
