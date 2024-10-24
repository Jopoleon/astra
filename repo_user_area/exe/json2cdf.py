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

    cdf_out = '%s.CDF' %loc

    j_json = 1
    f_json = '%s-%d.json' %(loc, j_json)

# Don't store anything if 1st json is older than tmp/astra.nml
    f_log2 = '%s/tmp/%s.nml' %(awd, expequ)

    if os.stat(f_json).st_mtime < os.stat(f_log2).st_mtime:
        logger.warning('json files older then %s' %f_log2)
        logger.warning('No 2D NetCDF file written')
        return

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
            astra_d = json_d['astra']
            equil_d = json_d['equil']
        nx   = astra_d['NEX']['ndim']
        n_eq = equil_d['rho_tor']['ndim']
        n_th = equil_d['teta2d']['ndim']
        nR   = equil_d['r2d']['ndim']
        nZ   = equil_d['z2d']['ndim']

        for key, val in astra_d.items():
            dat = val['data']
            if key not in ('XRHO', ):
                if j_json == 1:
                    if type(dat) == type([]): # list
                        if len(dat) > 0:
                            ds_astra[key] = dat
                        else: # array, all zeros
                            ds_astra[key] = val['ndim']*[0.]
                    else: # scalar
                        ds_astra[key] = [dat]
                else:
                    if type(dat) == type([]): # list
                        if len(dat) > 0:
                            ds_astra[key] += dat
                        else: # array, all zeros
                            ds_astra[key] += val['ndim']*[0.]
                    else: # scalar
                        ds_astra[key].append(dat)

        for key, val in equil_d.items():
            dat = val['data']
            if key not in ('rho_tor', 'teta2d'):
                if j_json == 1:
                    if type(dat) == type([]): # list
                        if len(dat) > 0:
                            ds_equil[key] = dat
                        else: # array, all zeros
                            ds_equil[key] = val['ndim']*[0.]
                    else: # scalar
                        ds_equil[key] = [dat]
                else:
                    if type(dat) == type([]): # list
                        if len(dat) > 0:
                            ds_equil[key] += dat
                        else: # array, all zeros
                            ds_equil[key] += val['ndim']*[0.]
                    else: # scalar
                        ds_equil[key].append(dat)
                        
        f_json_prev = f_json
        j_json += 1

    nt = j_json - 1
    logger.debug('nt=%d, nrho=%d, nr_eq=%d, nthe_eq=%d' %(nt, nx, n_eq, n_th))
    
    dtyp = '>f8' #np.float64
    for key, val in ds_astra.items():
        ds_astra[key] = np.array(val, dtype=dtyp)
        if len(val) == nt:
            astra_d[key]['dimensions'] = ['TIME']
        elif len(val) == nt*nx:
            ds_astra[key] = ds_astra[key].reshape((nt, nx))
            astra_d[key]['dimensions'] = ['TIME', 'XRHO']

    for key, val in ds_equil.items():
        ds_equil[key] = np.array(val, dtype=dtyp)
        if len(val) == nt:
            equil_d[key]['dimensions'] = ['TIME']
        elif len(val) == nt*n_eq:
            ds_equil[key] = ds_equil[key].reshape((nt, n_eq))
            equil_d[key]['dimensions'] = ['TIME', 'RHO_SURF']
        elif len(val) == nt*n_th*n_eq:
            ds_equil[key] = np.transpose(ds_equil[key].reshape((nt, n_th, n_eq)), (0, 2, 1) )
            equil_d[key]['dimensions'] = ['TIME', 'RHO_SURF', 'THETA']
        elif len(val) == nt*nR*nZ:
            ds_equil[key] = np.transpose(ds_equil[key].reshape((nt, nZ, nR)), (0, 2, 1) )
            equil_d[key]['dimensions'] = ['TIME', 'R', 'Z']

    f = netcdf_file(cdf_out, 'w', mmap=False)

    f.createDimension('TIME', nt)
    f.createDimension('XRHO', nx)
    f.createDimension('RHO_SURF', n_eq)
    f.createDimension('THETA', n_th)
    f.createDimension('R', nR)
    f.createDimension('Z', nZ)

    rho = f.createVariable('XRHO', dtyp, ('XRHO', ))
    rho.data  = np.array(astra_d['XRHO']['data'], dtype=dtyp)
    rho.units = astra_d['XRHO']['units']
    rho.long_name = astra_d['XRHO']['long_name']

    time = f.createVariable('TIME', dtyp, ('TIME', ))
    time.data = ds_astra['TIME'].astype(dtyp)
    time.units = 's'
    time.long_name = 'Time'

    rho_surf = f.createVariable('RHO_SURF', dtyp, ('RHO_SURF', ))
    rho_surf.data = np.array(equil_d['rho_tor']['data'], dtype=dtyp)
    rho_surf.units = '-'
    rho_surf.long_name = equil_d['rho_tor']['long_name']

    theta = f.createVariable('THETA', dtyp, ('THETA', ))
    theta.data = np.array(equil_d['teta2d']['data'], dtype=dtyp)
    theta.units = 'rad'
    theta.long_name = equil_d['teta2d']['long_name']

    rgrid = f.createVariable('R', dtyp, ('R', ))
    rgrid.data = np.array(equil_d['r2d']['data'], dtype=dtyp)
    rgrid.units = 'm'
    rgrid.long_name = equil_d['r2d']['long_name']

    zgrid = f.createVariable('Z', dtyp, ('Z', ))
    zgrid.data = np.array(equil_d['z2d']['data'], dtype=dtyp)
    zgrid.units = 'm'
    zgrid.long_name = equil_d['z2d']['long_name']

    for key, val in ds_astra.items():
        if key != 'TIME':
            tmp = f.createVariable(key, dtyp, astra_d[key]['dimensions'])
            tmp[:] = val
            tmp.units = astra_d[key]['units']
            tmp.long_name = astra_d[key]['long_name']

    for key, val in ds_equil.items():
        if key != 'TIME':
            if 'dimensions' in equil_d[key].keys():
                tmp = f.createVariable(key, dtyp, equil_d[key]['dimensions'])
                tmp[:] = val
                tmp.units = equil_d[key]['units']
                tmp.long_name = equil_d[key]['long_name']

    f.close()
    logger.info('Stored %s' %cdf_out)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='astra.nml writer')
    parser.add_argument('-e', '--expequ', help='<exp><equ>', required=True)
    args = parser.parse_args()
    json_concat(args.expequ)
