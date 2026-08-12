# E. Fable, E. Buglione-Ceresa (IPP-Garching) 2026
# try to run dkes interface (without rerunning vmec or boozer) and check if the coefficients agree with dkes output I have 
# when properly back converted.

# start to think that G11 in NEOTRANSP is just L11 of DKES with the geometric conversion, since the energy convolution seems to be done the same
# in PENTA code of stellopt that uses L11 of DKES.

# the sign of L13 has to be fixed

import os, sys, shutil, subprocess, glob, re
import numpy as np
from scipy.io import netcdf_file
from concurrent.futures import ProcessPoolExecutor, as_completed
from types import SimpleNamespace
from vmec import VMEC

stellopt_path = os.environ.get("STELLOPT_HOME")
os.environ["OMP_NUM_THREADS"] = "1"

MAX_WORKERS = min(96, os.cpu_count())

# =====================================
# USER SETTINGS
# =====================================
	
DKES_EXEC = os.path.join(stellopt_path, 'DKES', 'Release', 'xdkes')

#MAX_COUPLING_ORDER = 4   # user defined, change this to increase or decrease accuracy at low collisionality (higher = more computational time). 
MAX_COUPLING_ORDER = int(sys.argv[1]) if len(sys.argv) > 1 else 4

INPUT_FILE = "input.VMECoutput"
BOOZ_FILE = "boozmn_VMECoutput.nc"
BOOZ_INFILE = "inboozer.in"
VMEC_FILE = "wout_VMECoutput.nc"
VMEC_HEADER_FILE = "vmec_header_data_VMECoutput.txt"
with open(VMEC_HEADER_FILE, "r") as f:
    lines = f.readlines()
    third_line = lines[1].strip()
    psi_sign=np.sign(float(third_line))  # this is the one that the dkes interface is reading
print(third_line,psi_sign)

OUTPUT_NAME = "dkesout.VMECoutput"

#SURFACES = [13, 25, 37 ,49 ,61 ,73 ,85] #list(range(1, 11))
with open(BOOZ_INFILE, "r") as f:
    lines = f.readlines()

# pick the third line (index 2)
third_line = lines[2]

# split by whitespace and convert to integers
SURFACES = [int(x) for x in third_line.split()]
#SURFACES = [13, 85] #list(range(1, 11))

filename = "true_surfaces.txt"

# Check if the file exists
if not os.path.exists(filename):
    print(f"Error: file '{filename}' does not exist!")
    sys.exit(1)  # stop execution

# If it exists, read integers (space- or newline-separated)
with open(filename, "r") as f:
    next(f)
    TRUE_SURFACES = [int(x) for x in f.read().split()]

print("TRUE_SURFACES =", TRUE_SURFACES)

#TRUE_SURFACES = [ 2 , 5 , 13 , 25 , 41 , 60  ,85]   # this should come from user via file true_surfaces.txt where 
# a nr of surfaces and then a list of numbers should be given, boozer should have been run with ALL VMEC surfaces list, e.g. 
# 7
# 2 5 13 25 41 60 85

BOOZER_SURFACES = [s-1  for s in TRUE_SURFACES]

nr_surfaces = len(TRUE_SURFACES)
BOOZER_SURFACES2 = list(range(nr_surfaces))


SURFACES=TRUE_SURFACES

print(BOOZER_SURFACES2)
#exit()

CMUL_BLOCK_SIZE = 8   # number of (cmul,efield) rows per block

#read wout file
vmc = VMEC()
vmc.readWout(VMEC_FILE)
vmc.calcFtrap()

#open vmec and booz files
with netcdf_file(VMEC_FILE, 'r', mmap=False) as f:
    vmec = f.variables
with netcdf_file(BOOZ_FILE, 'r', mmap=False) as f:
    boox = f.variables

# VMEC quantities
Aminor = vmec["Aminor_p"].data  # this should be minor radius w7as ?
ns = vmec["ns"].data
print(ns)
#exit()
iotas = vmec["iotaf"].data
iotah = vmec["iotas"].data
phip = vmec["phi"].data
phip_p = vmec["phipf"].data
phip2 = booz["phi_b"].data
phip2_p = booz["phip_b"].data
ib = booz["iota_b"].data
iota_boozer=ib
iff = vmec["iotaf"].data
# need to find s_1 boozer and s_end boozer now!
s_vmec=phip/phip[-1]
s_boozer = s_vmec
s_boozer[0]=0
s_boozer[1:] = 0.5 * (s_vmec[1:] + s_vmec[:-1])

# note that in the boozer file, the half grid is defined for example so that (index starting from 1) iota_boozer(3)=mean(iota_vmec(2:3))
# iota_boozer(1) = 0, iota_boozer(2) = mean(iota_vmec(1:2)), etc. this means that in boozer file the surfaces are counted shifted by 1.



print(booz.keys(),s_vmec,s_boozer)
#exit()
# Note that boozer quantities are on the VMEC HALF GRID!, SO Boozer surface N is actually half surface between VMEC surface N and N+1, 
# for example boozer surface 2 is half grid between VMEC surface 2 and 3
print(phip/phip[-1],phip2/phip[-1],ib,iff, phip_p,phip2_p,iotas,iotah)
#exit()
#exit()

psia=phip[-1]/2./np.pi
psia_vmec=phip[-1]/2./np.pi
psiab=booz["phip_b"].data
# Boozer mode numbers
print(booz.keys())
#exit()
xm = vmec["xm"].data
xn = vmec["xn"].data
#print(xm,xn)
# Fourier amplitudes of |B|
bmnc = booz["bmnc_b"].data
rmnc = booz["rmnc_b"].data

# find index of m=0/1, n=0
xmb = booz["ixm_b"].data
xnb = booz["ixn_b"].data
nsb = booz["ns_b"].data
idx00 = np.where((xmb == 0) & (xnb == 0))[0][0]
idx01 = np.where((xmb == 1) & (xnb == 0))[0][0]

# choose flux surface index (example: last surface)

ib = booz["iota_b"].data
print(bmnc[:,idx00])
print(bmnc[:,idx01])
#exit()
B00_boozer = bmnc[BOOZER_SURFACES2, idx00]
R00_boozer = rmnc[BOOZER_SURFACES2, idx00]
B10_boozer = bmnc[BOOZER_SURFACES2, idx01]
print(B00_boozer,R00_boozer, B10_boozer,B10_boozer/B00_boozer,ib,nsb,len(ib),bmnc[:,idx00])

# note that bmnc and rmnc are length 98 instead of 99 like bvco since the last surface is skipped, 
# the first surface has rubbish reszults!!!!

#exit()
rmnc = vmec["rmnc"].data
zmns = vmec["zmns"].data

s = -1

R = rmnc[s]
Z = zmns[s]

a2 = 0.0

nmodes = len(xm)

for i in range(len(xm)):
    m = xm[i]

    if m == 0:
        continue

    a2 += m * R[i] * Z[i]

a_booz = np.sqrt(a2)

#a_booz = 0.49590 # benchmark vs old DKES case 20230214.055
Aminor = a_booz

#a_sum = 0.0
#parity = xn % 2
#for i in range(nmodes):
#    if xm[i] == 0:
#        continue
#    same_parity = (parity == parity[i])
#    a_sum += xm[i] * R[i] * np.sum(Z[same_parity])
#a_booz = np.sqrt(a_sum)

print(Aminor,psia,B00_boozer,a_booz,psia,psiab[-1])

#exit()

bmnc = booz["bmnc_b"][:]   # shape (ns, nmodes)
xn   = booz["ixn_b"][:]
xm   = booz["ixm_b"][:]
ns   = bmnc.shape[0]
print(ns)
#exit()
nmodes = bmnc.shape[1]

# Resolution for theta and phi grids
n_theta = 64
n_phi = 64
theta = np.linspace(0, 2*np.pi, n_theta, endpoint=False)
phi   = np.linspace(0, 2*np.pi, n_phi, endpoint=False)
theta_grid, phi_grid = np.meshgrid(theta, phi, indexing='ij')

def reconstruct_B(surface_index):
    """Reconstruct B(theta, phi) on the grid for a given surface"""
    B = np.zeros_like(theta_grid)
    for b, m, n in zip(bmnc[surface_index,:], xm, xn):
        B += b * np.cos(m*theta_grid - n*phi_grid)
    return B

def flux_surface_average(B):
    """Numerically compute <B^2> on the surface"""
    B2 = 1.0/B**2
    return np.mean(B2)  # simple average over theta, phi

# Compute <B^2> and beta for all surfaces
B2_avg_all = np.zeros(ns)
beta_num_all = np.zeros(ns)
for s in range(len(BOOZER_SURFACES2)):
    B_surf = reconstruct_B(s)
    B2_avg_all[s] = 1.0/flux_surface_average(B_surf)
    
print(B2_avg_all[BOOZER_SURFACES2]/B00_boozer**2)
print(B10_boozer/B00_boozer)

bmnc = booz["bmnc_b"][:]   # shape (ns, nmodes)
xn   = booz["ixn_b"][:]
xm   = booz["ixm_b"][:]
phi_booz   = booz["phi_b"][:]
ns   = bmnc.shape[0]

B2_avg = B2_avg_all[BOOZER_SURFACES2]

rrr=a_booz*np.sqrt(s_boozer[BOOZER_SURFACES])
epsilonz=rrr/R00_boozer
print(B10_boozer/B00_boozer,np.abs(B10_boozer/epsilonz/B00_boozer),rrr,R00_boozer,epsilonz,s_boozer[BOOZER_SURFACES])

# calcultae Er resonance, Er = iota * psip/R00 ; psip = B00*r
#Erresonance = 0.5*rrr/a_booz*iota_boozer[BOOZER_SURFACES]
Erresonance = np.abs(iota_boozer[BOOZER_SURFACES]/R00_boozer*rrr*B00_boozer**2.) # another B00 to make up for the normalization later on
print(Erresonance)

COMBINED_OUTPUT = "dkes_combined.dat"

CMUL_EFIELD_FILE = "cmul_efield_list.txt"

class Geomdat:
    pass

def safe_symlink(src, dst):
    """
    Create symbolic link, overwriting if it exists.
    """
    if os.path.exists(dst) or os.path.islink(dst):
        os.remove(dst)

    os.symlink(os.path.abspath(src), dst)


def split_cmul_efield_file(filename, block_size):

    with open(filename, "r") as f:
        lines = [l for l in f.readlines() if l.strip() != ""]

    blocks = []
    for i in range(0, len(lines), block_size):
        block = lines[i:i+block_size]
        ic = len(blocks)
        block_name = f"{filename}_{ic:02d}.txt"

        with open(block_name, "w") as out:
            out.writelines(block)

        blocks.append(block_name)

    return blocks

def parse_dkes_transport(text):
    results = []
    cmul = None
    efield = None
    Nmn = 0
    
    lines = text.splitlines()  # get all lines as a list
    for i, line in enumerate(lines):

        # detect CMUL / EFIELD header line
        if "MPOL" in line and "LALPHA" in line:
                numbers = lines[i+1].split()
                Lalpha = int(numbers[2])

        if "CHIP" in line and "PSIP" in line:
                numbers = lines[i+1].split()
                psip = float(numbers[1])

        if "M =" in line:
            m = re.search(r"M\s*=\s*(.*)", line)
            if m:
             m_values = m.group(1).split()
             Nmn += len(m_values)
		    
        if "CMUL" in line and "EFIELD" in line:
            # safely get the numeric line 3 lines below
            if i + 3 < len(lines):
                data_line = lines[i + 3].split()
                cmul = float(data_line[0])
                efield = float(data_line[1])


        # detect "+" transport row
        if line.strip().startswith("+"):
           parts = line.split()

           if len(parts) >= 5:
             L11p = float(parts[1])
             L33p = float(parts[2])
             L13p = float(parts[3])
             L31p = float(parts[4])
             next_parts = lines[i + 1].split()
             if len(next_parts) >= 5:
                 L11m = float(next_parts[1])
                 L33m = float(next_parts[2])
                 L13m = float(next_parts[3])
                 L31m = float(next_parts[4]) 

        elif line.strip().startswith("MKS FACT"):
           parts = line.split()

           m11 = float(parts[2])
           m33 = float(parts[3])
           m13 = float(parts[4])


        if "SPITZER" in line:
            # safely get the numeric line 3 lines below
            if i + 1 < len(lines):
                data_line = lines[i + 1].split()
                L33spitzer = float(data_line[2])

        # apply scaling
            L11 = 0.5*(L11p+L11m)
            L33 = 0.5*(L33p+L33m)
            L13 = 0.25*(L13p+L13m+L31p+L31m)
            DL11 = 0.5*np.abs(L11p-L11m)
            DL13 = 0.25*np.abs(L13p+L31p-L13m-L31m)
            DL33 = 0.5*np.abs(L33p-L33m)
            L33 = (L33spitzer-L33)
	    
            results.append((cmul, efield, L11, L33, L13, Lalpha, Nmn, DL11, DL13, DL33, psip))
            cmul = None
            efield = None
	
    return results
    
def get_dkdata(resultoo,conv11):
# needs cuml0, g110,g11e0, ug11_0, lg11_0, cmul_av, g11, efld
    """
    Collect DKES transport data for surface s into the dkdata structure.
    Assumes `results[s]` exists as a list of (ic, text) pairs per block.
    """

    transport_all = []

# parse all transport blocks for this surface
    for ic, text in resultoo:
     transport_block = parse_dkes_transport(text)
     transport_all.extend(transport_block)
    print(transport_all)

    if not transport_all:
     return None  # no data available

# --- CMUL_0 and L11 for zero E_r ---
    zero_efld = [t for t in transport_all if abs(t[1]) < 1e-12]
    if zero_efld:
     # collect all CMUL and L11 values for zero E_r
     cmul_0_list = [t[0] for t in zero_efld]
     g11_0_list  = [t[2]*conv11 for t in zero_efld]
     psip_list   = [t[10] for t in zero_efld]

    # approximate error from residual variations of L11
     g11e_0 = 1.e-11
     ug11_0 = max(g11_0_list)
     lg11_0 = min(g11_0_list)

    # optional averages
     cmul_0_avg = sum(cmul_0_list) / len(cmul_0_list)
     g11_0_avg  = sum(g11_0_list) / len(g11_0_list)
    else:
     cmul_0_list = g11_0_list = None
     g11e_0 = ug11_0 = lg11_0 = None
     cmul_0_avg = g11_0_avg = None

# --- average over all CMUL / EFIELD blocks ---
    cmul_vals = [t[0] for t in transport_all]
    g11_vals  = [t[2]*conv11 for t in transport_all]
    efld_vals = [t[1] for t in transport_all]
    psip = psip_list[0]
#    print(psip)
#    exit()
    cmul_av = cmul_vals
    g11_avg = g11_vals 
    efld_all = efld_vals  # list of all EFIELD values


    dkdata = {
     "cmul_0": cmul_0_list,       # list of CMUL values at zero field
     "g11_0": g11_0_list,         # list of L11 values at zero field
     "g11e_0": g11e_0,            # residual variation of L11
     "ug11_0": ug11_0,             # max L11
     "lg11_0": lg11_0,             # min L11
     "cmul_0_avg": cmul_0_avg,     # average CMUL at zero field
     "g11_0_avg": g11_0_avg,       # average L11 at zero field
     "cmul_av": cmul_av,           # average CMUL overall
     "g11": g11_avg,               # average L11 overall
     "efld": efld_all,             # all EFIELD values
     "psip": psip             # psip value
    }
    print(dkdata)
    return SimpleNamespace(**dkdata)
    
def dkes_fitd11_single(Geomdat, dkdata):

    R = Geomdat.R00
    epsilon = Geomdat.rnorm * Geomdat.minorradiusW7AS / R
    absiota = np.abs(Geomdat.iota)
    b0dk = 1.0

    cmul_0 = np.asarray(dkdata.cmul_0)
    g11_0  = np.asarray(dkdata.g11_0)

    if len(cmul_0) < 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    # DKES sign convention
    g11 = -g11_0 * cmul_0

    # -------------------------
    # 1/nu regime fit
    # -------------------------

    exp_fit = 0.5

    # choose last few points
    Nconsider = min(9, len(cmul_0) - 1)

    cmulexp = cmul_0[-Nconsider:] ** exp_fit
    g11_sub = g11[-Nconsider:]

    # detect inflection point
    dfdx = np.diff(g11_sub) / np.diff(cmulexp)
    x1 = (cmulexp[:-1] + cmulexp[1:]) / 2

    d2fdx2 = np.diff(dfdx) / np.diff(x1)
    x2 = (x1[:-1] + x1[1:]) / 2

    sinchange = np.where(np.diff(np.sign(d2fdx2)) != 0)[0]

    if sinchange.size == 0:
        Nfit = min(5, len(cmul_0))
    else:
        sci = sinchange[-1]
        inflection = x2[sci] - d2fdx2[sci] * (x2[sci+1] - x2[sci]) / (d2fdx2[sci+1] - d2fdx2[sci])
        indx = np.where(np.diff(np.sign(cmulexp - inflection)))[0]
        Nfit = max(4, Nconsider - indx[0]) if indx.size else 5

    # final fit
    x = cmul_0[-Nfit:] ** exp_fit
    y = g11[-Nfit:]

    weights = np.ones_like(x)

    y1, y0 = np.polyfit(x, y, 1, w=weights)


    while True:

        if Nfit < 2:
            # Emergency fallback
            y0 = g11[-1]
            y1 = 0.0
            break

        x = cmul_0[-Nfit:] ** exp_fit
        y = g11[-Nfit:]

        y1, y0 = np.polyfit(x, y, 1)

        if y0 >= 0:
            break

        # reduce number of points and retry
        Nfit -= 1



    g11fit = y1 * cmul_0**exp_fit + y0
    print(g11,Nfit,x,y,weights,y1,y0,g11fit,R,b0dk,g11_sub,cmulexp)
    eps_eff = (4.9982 * y0 * R**2 * b0dk**2)**(2.0/3.0)

    # -------------------------
    # sqrt(nu) regime fit
    # -------------------------

    cmul = np.asarray(dkdata.cmul_av)
    efld = np.asarray(dkdata.efld)
    efld=efld/Geomdat.B00   # DKES was run with Er not normalized to B00!
    g11_all = np.abs(np.asarray(dkdata.g11))

    er_res = absiota * epsilon
    er_u = 0.05 * er_res

    mask = (efld > 0) & (efld <= er_u)

    good_cmul = cmul[mask]
    good_efld = efld[mask]
    good_d11  = g11_all[mask]

    if len(good_d11) < 5:
        return eps_eff, y0, np.nan, np.nan, er_u

    # grid search
    ex_ers = np.linspace(0.5, 8.0, 120)
    g11_ers = np.logspace(-10, -2, 120)

    best_cost = np.inf
    best_ex = np.nan
    best_g11 = np.nan

    for ex_er in ex_ers:
        for g11_er in g11_ers:

            d11_sqrtnu = np.sqrt(good_cmul) * g11_er / good_efld**1.5
            d11_nuinv = y0 / good_cmul

            model = (d11_nuinv**(-ex_er) + d11_sqrtnu**(-ex_er))**(-1/ex_er)

            cost = np.sum((np.log(model) - np.log(good_d11))**2)

            if cost < best_cost:
                best_cost = cost
                best_ex = ex_er
                best_g11 = g11_er

    return eps_eff, y0, best_g11, best_ex, er_u
    
def extract_surface_quantities2(vmec_file, booz_file, s, resultoo):

    # VMEC quantities
    Aminor = a_booz  # this should be minor radius w7as ?
    iotas = iota_boozer

    idx = SURFACES[s] - 1

        # normalized toroidal flux
    s_norm = s_boozer[idx] #idx / (ns - 1)
    print(s_norm)
    
        # minor radius
    r = Aminor * np.sqrt(s_norm)

        # major radius
    R00 = R00_boozer[s]

        # B00
    B00 = B00_boozer[s]
    B00real = B00
    B00 = 1.0

        # rotational transform
    iota = iotas[idx]

    Bsq_avg = B2_avg[s]

        # flux surface average of B^2/B0^2
    Bsq = Bsq_avg/B00real**2. #bdotb[idx]/bmnc[idx, 0]**2.

        # ftrap
    ftrap = vmc.ftrap[idx]
	
        # Kn
    epsilonz=r/R00
#    print(idx,ind10)
    b10= B10_boozer[s]/B00_boozer[s]   # this should be VMEC Bnorm of m=1, n=0 at radius idx
    print(b10,epsilonz,b10/epsilonz,r,R00,epsilonz)
    kn = np.abs(b10/epsilonz)

        # geometric data
    geom=Geomdat()
    		
    geom.R00 = R00
    geom.rnorm = r/Aminor
    geom.minorradiusW7AS = Aminor
    geom.iota = iota
    geom.B00=B00real

    drtildedr=2.*psia/(B00real*Aminor**2)
    print(drtildedr)
    conv11=-1.0*B00real**2*np.abs(drtildedr)**0.
    conv13=-np.sqrt(Bsq)*B00real**1*np.sqrt(np.abs(drtildedr))**0.   # whats the correct sign here? not sure.
    conv33=-Bsq    # L33 is now PERFECT at the first CMUL
	
    dkdata = get_dkdata(resultoo,conv11) # get dkes data
    conv13=-conv13*np.sign(dkdata.psip)*psi_sign
    	
        # calculate the fits
    eps_eff, g11_ft, g11_er, ex_er, efield_u = dkes_fitd11_single(geom, dkdata)	




#collect results
    results = {
            "surface": s,
            "r": float(r),
            "R00": float(R00),
            "B00": float(B00),
            "iota": float(iota),
            "ftrap": float(ftrap),
            "kn": float(kn),
            "<B^2>": float(Bsq),
            "eps_eff": float(eps_eff),
            "g11_ft": float(g11_ft),
            "efield_u": float(efield_u),
            "g11_er": float(g11_er),
            "ex_er": float(ex_er),
            "conv11": float(conv11),
            "conv13": float(conv13),
            "conv33": float(conv33)
    }


    return results


def get_psip_for_dkes(surface):

    # VMEC quantities
    Aminor = a_booz  # this should be minor radius w7as ?

    psia=psia_vmec
    print(Aminor,psia)
#    exit()

    s = SURFACES.index(surface)
    idx = SURFACES[s]-1

        # normalized toroidal flux
    s_norm = s_boozer[idx]
    print(s_norm)
    
        # minor radius
    r = Aminor * np.sqrt(s_norm)

        # major radius
    R00 = R00_boozer[s]

        # B00
    B00 = B00_boozer[s]
    B00real = B00
    B00 = 1.0

        # rotational transform
    iota = iotas[idx]

    Bsq_avg = B2_avg[s]

        # flux surface average of B^2/B0^2
    Bsq = Bsq_avg/B00real**2. #bdotb[idx]/bmnc[idx, 0]**2.

#    psippo = r * 2 * psia / Aminor**2 * 1.0 /B00real
    psippo = r  * B00real
    

    return psippo, B00real


def run_dkes_surface_block(surface, ic, block_file):

    run_dir = f"surface_{surface:02d}_{ic:02d}"
    os.makedirs(run_dir, exist_ok=True)

    # Copy required files
    shutil.copy(INPUT_FILE, os.path.join(run_dir, INPUT_FILE))
    shutil.copy(block_file, os.path.join(run_dir, CMUL_EFIELD_FILE))
    shutil.copy(BOOZ_FILE, os.path.join(run_dir, BOOZ_FILE))
    shutil.copy(VMEC_FILE, os.path.join(run_dir, VMEC_FILE))
# Symlink required files
    safe_symlink(INPUT_FILE,
             os.path.join(run_dir, INPUT_FILE))

#    safe_symlink(block_file,
#             os.path.join(run_dir, CMUL_EFIELD_FILE))

    safe_symlink(BOOZ_FILE,
             os.path.join(run_dir, BOOZ_FILE))

    safe_symlink(VMEC_FILE,
             os.path.join(run_dir, VMEC_FILE))
    legend_file = os.path.join(run_dir, "legendrnm.txt")

    # Example: define Lalpha, M, N here (or compute them from your data)
    Lalpha = 300   # pitch-angle grid
    M_values = [44]  # example M modes
    N_values = [44]  # example N modes

    with open(legend_file, "w") as f:
       f.write(f"{Lalpha:6d}\n")
       
    psip_file = os.path.join(run_dir, "psi_def_for_dkes.txt")

    # get psip in new definition
    surface_index = TRUE_SURFACES.index(surface)  # 0-based index    
    psippo, B00real = get_psip_for_dkes(surface)   
    actual_js = TRUE_SURFACES[surface_index]

    with open(psip_file, "w") as f:
       f.write(f"{psippo:7.4f} {B00real:7.4f} {MAX_COUPLING_ORDER}  {actual_js}\n")


    log_file = os.path.join(run_dir, "dkes.log")


#rescale electric field
    file_path = os.path.join(run_dir, CMUL_EFIELD_FILE)
    data = np.loadtxt(file_path)
    data[:, 1] *= Erresonance[surface_index]
    np.savetxt(file_path, data)   # cmul_efield_list.txt full from dat has to have Er going from 0. to 1.
    
    with open(log_file, "w") as log:
        subprocess.run(
            [DKES_EXEC, INPUT_FILE, str(surface_index+2)],
            cwd=run_dir,
            stdout=log,
            stderr=log
        )

    output_file = os.path.join(run_dir, OUTPUT_NAME)

    if not os.path.exists(output_file):
        print(f"Surface {surface} block {ic} failed")
        return surface, ic, None

    with open(output_file, "r") as f:
        data = f.read()

    return surface, ic, data

def collect_existing_dkes_outputs():

    results = {}

    for run_dir in sorted(glob.glob("surface_*_*")):

        m = re.match(r"surface_(\d+)_(\d+)", os.path.basename(run_dir))
        if not m:
            continue

        surface = int(m.group(1))
        ic = int(m.group(2))

        output_file = os.path.join(run_dir, OUTPUT_NAME)

        if not os.path.exists(output_file):
            print(f"Missing {output_file}")
            continue

        with open(output_file, "r") as f:
            data = f.read()

        if surface not in results:
            results[surface] = []

        results[surface].append((ic, data))

    # sort blocks
    for s in results:
        results[s].sort(key=lambda x: x[0])

    return results


# =====================================
# MAIN PARALLEL BLOCK
# =====================================

def main(run_dkes):

    if run_dkes:
# remove directories surface_*
     for d in glob.glob("surface_*"):
      if os.path.isdir(d):
        shutil.rmtree(d)

# remove files cmul_efield_list.txt_*
     for f in glob.glob("cmul_efield_list.txt_*"):
      if os.path.isfile(f):
        os.remove(f)    
    # ---------------------------------
    # Split cmul/efield file
    # ---------------------------------

     block_files = split_cmul_efield_file(
         CMUL_EFIELD_FILE,
         CMUL_BLOCK_SIZE
     )

     tasks = []
     for s in SURFACES:
         for ic, bf in enumerate(block_files):
             tasks.append((s, ic, bf))

     results = {}

    # ---------------------------------
    # Parallel execution
    # ---------------------------------
     with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
#    with ProcessPoolExecutor() as executor:

         futures = [
             executor.submit(run_dkes_surface_block, s, ic, bf)
             for (s, ic, bf) in tasks
         ]

         for future in as_completed(futures):

            surface, ic, data = future.result()

            if data is not None:

                if surface not in results:
                    results[surface] = []

                results[surface].append((ic, data))

    # ---------------------------------
    # Sort blocks per surface
    # ---------------------------------
     for s in results:
        results[s].sort(key=lambda x: x[0])

#    eqdata = extract_surface_quantities(VMEC_FILE, BOOZ_FILE, SURFACES)


    else:
     results = collect_existing_dkes_outputs()


    # ---------------------------------
    # Write combined file
    # ---------------------------------
    with open(COMBINED_OUTPUT, "w") as out:

        out.write("cc   ref274: configuration nickname\n")
        out.write("cc   DKES database generated automatically\n")
        out.write(f"cc   VMEC file: {VMEC_FILE}\n")
        out.write("cc   DKES database generated automatically\n")
        out.write("cc   DKES database generated automatically\n")
        out.write("cc   DKES database generated automatically\n")

        for i, surface in enumerate(sorted(results.keys())):

        #    eq = eqdata[i]
            eq = extract_surface_quantities2(VMEC_FILE, BOOZ_FILE, i, results[surface])

            # Collect transport from all blocks
            transport_all = []

            for ic, text in results[surface]:
                transport_block = parse_dkes_transport(text)
                transport_all.extend(transport_block)

            # FIRST LINE
            out.write(
                f"{eq['r']:7.4f}  "
                f"{eq['R00']:7.4f}  "
                f"{eq['B00']:4.1f}  "
                f"{eq['iota']:7.4f} "
                f"{eq['kn']:7.5f} "
                f"{eq['ftrap']:7.5f} "
                f"{eq['<B^2>']:7.5f}  "
                " r,R,B,io,xkn,ft,<b^2>\n"
            )

            if i == 0:
                out.write(
                    " c         eps_eff     g11_ft      efield_u    g11_er    ex_er\n"
                )
    
            out.write(
              f" cfit  "
              f"{eq['eps_eff']:10.3E}  "
              f"{eq['g11_ft']:10.3E}  "
              f"{eq['efield_u']:10.3E}  "
              f"{eq['g11_er']:10.3E}  "
              f"{eq['ex_er']:6.3f}\n"
            )
            # Write combined transport
            for row in transport_all:

                cmul, efield, L11, L33, L13, Lalpha, Nmn, dl11, dl13, dl33, psip = row

# Need to convert L11,L33,L31 into D11,D13,D33 of Beidler et al. as explained in the notes of H. Smith.
                L11 = L11*eq["conv11"]
                L13 = L13*eq["conv13"]
                L33 = L33*eq["conv33"]
                dl11= dl11*np.abs(eq["conv11"])
                dl13= dl13*np.abs(eq["conv13"])
                dl33= dl33*np.abs(eq["conv33"])
                efield = efield/B00_boozer[i]

                out.write(
                    f"{cmul:12.3E}"
                    f"{efield:12.3E}"
                    f"{L11:14.5E}"
                    f"{L13:14.5E}"
                    f"{L33:14.5E}\n"
                )
                out.write(
                    f">3"
                    f"{Nmn:12d}"
                    f"{Lalpha:6d}"
                    f"{dl11:14.5E}"
                    f"{dl13:14.5E}"
                    f"{dl33:14.5E}\n"
                )

            out.write(" e\n")

    print("Finished")


# =====================================

if __name__ == "__main__":
#    main(False)    # just do calculation of fits
    main(True)    # do entire DKES calculation
