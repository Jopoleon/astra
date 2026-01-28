import os, logging, argparse
from scipy.io import netcdf_file
import numpy as np

awd = os.path.dirname(os.path.realpath(__file__))
fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('regression')

def compare(fcdf, fcdf_ref, tolerance=1.e-7):

    ncnew = netcdf_file(fcdf    , 'r', mmap=False).variables
    ncref = netcdf_file(fcdf_ref, 'r', mmap=False).variables

    time = ncnew['TIME'].data
    nt = len(time)

    print(fcdf)
# Check keys
    for key in ncref.keys():
        if key not in ncnew.keys():
            logger.error('Key %s is missing', key)
    for key, arr1 in ncnew.items():
        logger.info(key)
        if key not in ncref.keys():
            logger.error('Key %s not in reference nc-file %s', key, fcdf_ref)
        else:
            arr2 = ncref[key].data
            if arr1.shape != arr2.shape:
                logger.error('Shape mismatch')
            else:
                for jt in range(nt):
                    arr_new = np.atleast_1d(arr1[jt])
                    arr_ref = np.atleast_1d(arr2[jt])
                    nlen = len(arr_new)
                    norm_new = np.sum(np.abs(arr_new))
                    norm_ref = np.sum(np.abs(arr_ref))
                    if norm_new == 0:
                        if norm_ref == 0:
                            diff = 0
                        else:
                            diff = np.linalg.norm(arr_new - arr_ref)/norm_ref
                    else:
                        diff = np.linalg.norm(arr_new - arr_ref)/norm_new
                    if diff > tolerance:
                        logger.error('%s: discrepancy %12.4e at time=%8.4f', key, diff, time[jt])


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Regression test')
    
    parser.add_argument('-m', '--equ', help='Model file', required=False, default='fluxes')
    parser.add_argument('-v', '--exp', help='Exp file', required=False, default='aug34954')
    args = parser.parse_args()

    reg_log = '%s/regressions/%s%s.log' %(awd, args.exp, args.equ)
    if len(logger.handlers) == 0:
        hnd = logging.FileHandler(reg_log, mode='w')
        hnd.setFormatter(fmt)
        logger.addHandler(hnd)
        logger.setLevel(logging.ERROR)

    fcdf = '%s%s.CDF' %(args.exp, args.equ)
    fnew = 'ncdf_out/%s'    %fcdf
    fref = 'regressions/%s' %fcdf
    compare(fnew, fref)
    print('Log report in %s' %reg_log)
