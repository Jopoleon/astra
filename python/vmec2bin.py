import sys, os, argparse, logging, time
import numpy as np
from scipy.interpolate import interp1d
from pellet import PELLET
from vmec import VMEC, f2h
from parse_fortran_nml import parse_fortran_namelist

GP2 = 2.*np.pi

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('vmec2a')
if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)
#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

# Interpolation
linterp = lambda x, y, xi: interp1d(x, y, kind='linear'   , fill_value='extrapolate')(xi)
qinterp = lambda x, y, xi: interp1d(x, y, kind='quadratic', fill_value='extrapolate')(xi)


def vmec2astra(vmc, NA1, metric_file, b00_output_file):
    
    mn_index_b00 = np.where((vmc.xm_nyq == 0) & (vmc.xn_nyq == 0))[0]

    if mn_index_b00.size == 0:
        logger.error("ERROR: Could not find the (m=0, n=0) mode for 'bmnc'.")
        sys.exit(1)

    B00_physical_profile = vmc.bmnc[:, mn_index_b00].squeeze()

#    Ip_contrib = -vmc.phi_abs[-1]*s31[-1]/(GP2*vmc.BTOR*vmc.RHOVMECmesh)
#    F0_contrib = -vmc.phi_abs[-1]*s32[-1]/(GP2*vmc.BTOR*vmc.RHOVMECmesh)
    Ip_contrib = np.abs(0.*vmc.rbtor)
    F0_contrib = np.abs(0.*vmc.rbtor)

# Constants from ASTRA

    HROX = 1./(NA1 - 0.5)

# Define ASTRA arrays
    SXHO = HROX*np.arange(1, NA1+1)
    XRHO = SXHO - 0.5*HROX

#some quantities from VMEC output
    lasym  = vmc.lasym
    RTOR   = vmc.rmajor
    ABC    = vmc.aminor
    volume = vmc.volume

# Astra grids
    ROC  = vmc.RHOVMECmesh[-1]
    RHO  = XRHO*ROC
    SRHO = SXHO*ROC
    HRO  = RHO[1] - RHO[0]

    GRADROVMEC = vmc.avg_grad_rho * ROC

# interpolate to ASTRA grid
    SG11   = linterp(vmc.RHOVMECmesh, vmc.S11a, SRHO)
    SG12   = linterp(vmc.RHOVMECmesh, vmc.S12a, SRHO)
    SG21   = linterp(vmc.RHOVMECmesh, vmc.S21a, SRHO)
    SG22   = linterp(vmc.RHOVMECmesh, vmc.S22a, SRHO)
    VR     = linterp(vmc.RHOVMECmesh, vmc.Vpf, RHO)
    VRS    = linterp(vmc.RHOVMECmesh, vmc.Vpf, SRHO)
    GRADRO = linterp(vmc.RHOVMECmesh, GRADROVMEC, SRHO)
    G11    = linterp(vmc.RHOVMECmesh, vmc.g11, SRHO)
    Rmaj   = linterp(vmc.RHOVMECmesh, vmc.Rmaj, RHO)
    VOLUM  = linterp(vmc.RHOVMECmesh, vmc.Vol, RHO)
    AMETR  = linterp(vmc.RHOVMECmesh, vmc.Amineff, RHO)
    FTPT   = linterp(vmc.RHOVMECmesh, vmc.ftrap, RHO)
    CU     = linterp(vmc.RHOVMECmesh, vmc.jpar, RHO)
    MU     = linterp(vmc.RHOVMECmesh, vmc.iotaf, SRHO)
    SLAT   = VR*GRADRO

# Approximate values at NA1 for shifted grid
    def extrapolate_end(arr):
        return 1.5 * arr[-2] - 0.5 * arr[-3]

    SG11[-1] = extrapolate_end(SG11)
    SG12[-1] = extrapolate_end(SG12)
    SG21[-1] = extrapolate_end(SG21)
    SG22[-1] = extrapolate_end(SG22)
    MV = -SG12/SG11
    VRS[-1]    = extrapolate_end(VRS)
    GRADRO[-1] = extrapolate_end(GRADRO)
    G11[-1]    = extrapolate_end(G11)
    MU[-1]     = extrapolate_end(MU)
    xm = vmc.xm.astype(np.int32)   # poloidal mode numbers
    xn = vmc.xn.astype(np.int32)   # toroidal mode numbers

    nn = len(vmc.raxis_cc)

    logger.info(f"Major radius (RTOR): {vmc.rmajor:.4f} m")
    logger.info(f"Minor radius (ABC and AB):  {vmc.aminor:.4f} m")
    logger.info(f"Toroidal field on axis (BTOR): {vmc.b0:.4f} T")

# Writing bin file for ASTRA
    with open(metric_file, 'wb') as f:

        np.array(
            [
                HROX, HRO, ROC, RTOR, ABC, vmc.BTOR,
                volume, vmc.GVAC, vmc.Fboundary,
                vmc.nfp, vmc.dphidsb,
                Ip_contrib, F0_contrib
            ],
            dtype=np.float64
        ).tofile(f)

        for arr in [
            RHO, SRHO, SG11, SG12,
            SG21, SG22, MV, VR,
            VRS, GRADRO, G11, Rmaj,
            VOLUM, AMETR, SLAT, FTPT
        ]:
            arr.tofile(f) # They are np.float64, checked

        np.array([vmc.mnmax], dtype=np.int32).tofile(f)

        xm.tofile(f)
        xn.tofile(f)

        vmc.rmnc[-1, :].tofile(f)
        vmc.zmns[-1, :].tofile(f)
        np.array(lasym, dtype=np.int32).tofile(f)
        if lasym:
            vmc.rmns[-1, :].tofile(f)
            vmc.zmnc[-1, :].tofile(f)

        np.array([nn], dtype=np.int32).tofile(f)
        vmc.raxis_cc.tofile(f)
        vmc.zaxis_cs.tofile(f)

        ns_all = vmc.rmnc.shape[0]   # number of flux surfaces
        np.array([ns_all], dtype=np.int32).tofile(f)

        vmc.rmnc.tofile(f)
        vmc.zmns.tofile(f)
        if lasym:
            vmc.rmns.tofile(f)
            vmc.zmnc.tofile(f)
        vmc.phi.tofile(f)

    logger.info(f'Written {metric_file}')

#----------------------
    B00_ASTRA_grid = linterp(vmc.RHOVMECmesh, B00_physical_profile, RHO)
    np.savetxt(b00_output_file, B00_ASTRA_grid)
    logger.info(f'Written B00 profile to {b00_output_file}')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "astra_nrad", type=int, nargs="?", default=91, help="Number of ASTRA radial grid points")

    args = parser.parse_args()

    NA1 = args.astra_nrad

    namelist_path = f'{awd}/vmec_io/stell_files.nml'
    nl_vmec = parse_fortran_namelist(namelist_path, 'VMEC_TO_ASTRA_INPUTS')
    nl_dkes = parse_fortran_namelist(namelist_path, 'ASTRA_DKES_INTERFACE')
    vmec_wd     = f'{awd}/{nl_vmec["vmec_wd"]}'
    wout_file   = f'{vmec_wd}/{nl_vmec["vmec_wout_file"]}'
    metric_file = f'{vmec_wd}/{nl_vmec["vmec2a_metric"]}'
    b00_output_file = f'{vmec_wd}/{nl_dkes["b00_profile_file"]}'
    header_file     = f'{vmec_wd}/{nl_dkes["vmec_header_file"]}'
    radius_out      = f'{vmec_wd}/{nl_dkes["minor_radius_w7as_file"]}'

    logger.info('Reading VMEC NetCDF output')
    vmc = VMEC(wout_file)
    logger.info('Calculating VMEC moments')
    vmc.calcMoms(nu=64, nv=128)
    logger.info('Calculating VMEC algebra')
    vmc.calcAll()
#---------------------- VMEC output files
    vmc.write_header(header_file)
    vmc.write_amin(radius_out)

# NGS pellet: if a launch is configured, trace the straight-line chord for
# this (updated) equilibrium and refresh the table the Fortran model reads.
    pel = PELLET()
    pel.parse_pellet_nml()
    pel.chords(vmc)

    vmec2astra(vmc, NA1, metric_file, b00_output_file)
