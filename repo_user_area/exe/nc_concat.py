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
    n_eq, n_the = cv['r2d'].shape

    for key, val in ds.items():
        if len(val) == nt*nx:
            ds[key] = ds[key].reshape((nt, nx))
        if len(val) == nt*n_eq*n_the:
            ds[key] = ds[key].reshape((nt, n_eq, n_the))

    f = netcdf.netcdf_file(cdf_out, 'w', mmap=False)

    f.createDimension('TIME', nt)
    f.createDimension('XRHO', nx)
    f.createDimension('RHO_SURF', n_eq)
    f.createDimension('THETA', n_the)

    rho = f.createVariable('XRHO', np.float32, ('XRHO', ))
    rho.data = cv['XRHO'].data
    rho.units = cv['XRHO'].units
    rho.long_name = cv['XRHO'].long_name

    time = f.createVariable('TIME', np.float32, ('TIME', ))
    time.data = ds['TIME']
    time.units = 's'
    time.long_name = 'Time'

    for key, val in ds.items():
        if key not in ('XRHO', 'TIME', 'RHO_SURF', 'THETA'):
            dims = ('TIME', ) + cv[key].dimensions
            tmp = f.createVariable(key, np.float32, dims)

            tmp[:] = val
            tmp.units = cv[key].units
            tmp.long_name = cv[key].long_name

    f.close()
    logger.info('Stored %s' %cdf_out)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='astra.nml writer')
    parser.add_argument('-e', '--expequ', help='<exp><equ>', required=True)
    args = parser.parse_args()

    nc_concat(args.expequ)
