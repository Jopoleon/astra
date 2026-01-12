import os, logging, argparse
from scipy.io import netcdf_file
import numpy as np
import matplotlib.pylab as plt

awd = os.path.dirname(os.path.realpath(__file__))

fname = 'aug34954fluxes'
fcdf_ref = '%s/regressions/%s.CDF' %(awd, fname)
fcdf_new = '%s/ncdf_out/%s.CDF'    %(awd, fname)

ncref = netcdf_file(fcdf_ref, 'r', mmap=False).variables
ncnew = netcdf_file(fcdf_new, 'r', mmap=False).variables

rho = ncref['XRHO'].data

t_ref = ncref['TIME'][-1]
var = 'pprime'
yref = ncref[var][-1]
ynew = ncnew[var][-1]

plt.plot(rho, yref)
plt.plot(rho, ynew)
plt.xlabel(r'$\rho_{tor}$')
plt.ylabel(var)
plt.figtext(0.5, 0.85, 'time = %6.3f s' %t_ref, ha='center')
plt.show()
