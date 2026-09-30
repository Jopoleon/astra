#!/usr/bin/env python3
import numpy as np
from scipy.io import netcdf_file

VMEC_FILE = "dat/wout_VMECoutput.nc"

# Open VMEC file
with netcdf_file(VMEC_FILE, 'r') as f:
    vmec = f.variables

# Load poloidal and toroidal mode numbers
xm = vmec["xm"][:]  # poloidal mode numbers
xn = vmec["xn"][:]  # toroidal mode numbers
print(len(vmec["phi"][:]))

# Find unique modes
unique_m = np.unique(xm)
unique_n = np.unique(xn)

print("Total poloidal modes (unique m):", len(xm))
print("Poloidal mode numbers:", unique_m)
print("Total toroidal modes (unique n):", len(xn))
print("Toroidal mode numbers:", unique_n)
