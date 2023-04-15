#!/usr/bin/env python

import os, logging, argparse
from scipy.io import netcdf
import numpy as np

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')

logger = logging.getLogger('nc_concat')

if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)

#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

awd = os.getenv('AWD')
if awd is None:
    awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def nc_concat(expequ):

    loc = '%s/.res/ncdf/%s' %(awd, expequ)

    cdf_out = '%s.CDF' %loc
    ds = {}

    j_cdf = 1
    f_cdf = '%s%d.cdf' %(loc, j_cdf)
# Don't store anything if 1st cdf is older than tmp/astra.nml
    f_log2 = '%s/tmp/%s.nml' %(awd, expequ)

    if os.stat(f_cdf).st_mtime < os.stat(f_log2).st_mtime:
        logger.warning('CDF files older then %s' %f_log2)
        logger.warning('No 2D NetCDF file written')
        os.system('rm %s' %f_log2)
        return

# run as a second nc_concat
#    os.system('rm %s' %f_log2)

    while True:
        f_cdf = '%s%d.cdf' %(loc, j_cdf)
        if not os.path.isfile(f_cdf):
            if j_cdf == 1:
                return
            else:
                break
        if j_cdf > 1:
            if os.stat(f_cdf).st_mtime < os.stat(f_cdf_prev).st_mtime: #Newer
                break
        cv = netcdf.netcdf_file(f_cdf, 'r', mmap=False).variables
        for key, val in cv.items():
            if j_cdf == 1:
                ds[key] = val.data
            else:
                ds[key] = np.append(ds[key], val.data)
        f_cdf_prev = f_cdf
        j_cdf += 1

    nt = len(ds['TIME'])
    logger.info('nt = %d, n_cdf = %d' %(nt, j_cdf-1))
    nx = len(cv['XRHO'].data)
    n_eq = len(cv['RHO_SURF'].data)
    n_th = len(cv['THETA'].data)

    for key, val in ds.items():
        if key not in ('XRHO', 'THETA', 'RHO_SURF'):
            if cv[key].dimensions == ('XRHO', ):
                ds[key] = ds[key].reshape((nt, nx))
            if cv[key].dimensions == ('THETA', 'RHO_SURF'):
                ds[key] = ds[key].reshape((nt, n_th, n_eq))

    f = netcdf.netcdf_file(cdf_out, 'w', mmap=False)

    f.createDimension('TIME', nt)
    f.createDimension('XRHO', nx)
    f.createDimension('RHO_SURF', n_eq)
    f.createDimension('THETA', n_th)

    rho = f.createVariable('XRHO', np.float64, ('XRHO', ))
    rho.data = cv['XRHO'].data
    rho.units = cv['XRHO'].units
    rho.long_name = cv['XRHO'].long_name

    print('nc_concat:TIME', ds['TIME'])
    time = f.createVariable('TIME', np.float64, ('TIME', ))
    time.data = np.float64(ds['TIME'])
    time.units = 's'
    time.long_name = 'Time'

    rho_surf = f.createVariable('RHO_SURF', np.float64, ('RHO_SURF', ))
    rho_surf.data = cv['RHO_SURF'].data
    rho_surf.units = '-'
    rho_surf.long_name = cv['RHO_SURF'].long_name

    theta = f.createVariable('THETA', np.float64, ('THETA', ))
    theta.data = cv['THETA'].data
    theta.units = 'rad'
    theta.long_name = cv['THETA'].long_name

    for key, val in ds.items():
        if key not in ('XRHO', 'TIME', 'RHO_SURF', 'THETA'):
            if 'TIME' in cv[key].dimensions:
                dims = cv[key].dimensions
            else:
                dims = ('TIME', ) + cv[key].dimensions
            tmp = f.createVariable(key, np.float64, dims)
            tmp.data = val.astype(np.float64)
            tmp.units = cv[key].units
            tmp.long_name = cv[key].long_name

    f.close()
    logger.info('Stored %s' %cdf_out)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='astra.nml writer')
    parser.add_argument('-e', '--expequ', help='<exp><equ>', required=True)
    args = parser.parse_args()

    nc_concat(args.expequ)
