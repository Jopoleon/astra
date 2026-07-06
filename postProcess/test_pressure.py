import os, datetime, logging, argparse
from scipy.io import netcdf_file
import matplotlib.pylab as plt

keV_to_eV    = 1e3
e19m3_to_m3  = 1e19
eV_m3_to_Pa = 1.602*1e-19
keV_m3_to_Pa = 1.602*1e-16
keV_e19m3_to_Pa = e19m3_to_m3*keV_m3_to_Pa

fcdf = '../ncdf_out/aug34954flux_feqis.CDF'
cv = netcdf_file(fcdf, 'r', mmap=False).variables

jt = -1
ne = e19m3_to_m3*cv['NE'][jt, :]
ni = e19m3_to_m3*cv['NI'][jt, :]
nd = e19m3_to_m3*cv['NDEUT'][jt, :]
te = keV_to_eV*cv['TE'][jt, :]
ti = keV_to_eV*cv['TI'][jt, :]

rho = cv['XRHO'].data

pres_e = eV_m3_to_Pa*ne*te
pres_i_th = eV_m3_to_Pa*ni*ti
pres_i_fast = 0.5*keV_e19m3_to_Pa*cv['PBPER'][jt, :] + 0.5*keV_e19m3_to_Pa*cv['PBLON'][jt, :]
pres_equil = cv['pressure'][jt, :]
pres_tot = pres_e + pres_i_th + pres_i_fast
print(pres_equil)
print(pres_tot)
print(pres_i_fast)

plt.plot(rho, pres_e     , 'g-', label='e')
plt.plot(rho, pres_i_th  , 'k-', label='i,th')
plt.plot(rho, pres_i_fast, 'm-', label='i,fast')
plt.plot(rho, pres_tot   , 'r-', label='e+i_th+0.5*fast')
plt.plot(rho, pres_equil , 'b--', label='equil')
plt.legend()
plt.show()
