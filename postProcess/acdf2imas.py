import os, datetime, logging, argparse
from scipy.io import netcdf_file
import imas
import numpy as np

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('acdf2imas')
if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)
#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

# Units conversion factors

keV_m3_to_Pa = 1.602*1e-16
e19m3_to_m3  = 1e19
keV_to_eV    = 1e3


def ACDF2IMAS(args, write_ids=True):

    logger.info('Creating IDS structure')

    if args.ids_backend == 'HDF5':
        db = imas.DBEntry(imas.imasdef.HDF5_BACKEND   , 'aug', args.shot, args.ids_run, os.getenv('IMASDB'), '3')
    elif args.ids_backend == 'MDS+':
        db = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND, 'aug', args.shot, args.ids_run, os.getenv('IMASDB'), '3')
    elif args.ids_backend == 'ASCII':
        db = imas.DBEntry(imas.imasdef.ASCII_BACKEND  , 'aug', args.shot, args.ids_run, os.getenv('IMASDB'), '3')
    status, _ = db.create()

    cv = netcdf_file(args.fcdf, 'r', mmap=False).variables
    cp = fill_core_profiles(cv)
    eq = fill_equilibrium(cv)

    if write_ids:
        logger.info('Dumping IDS file %s', args.ids_backend)
        logger.info('in dir %s' %os.getenv('IMASDB'))
        logger.info('Putting core_profiles')
        db.put(cp)
        logger.info('Putting equilibrium')
        db.put(eq)
        logger.info('Closed IMAS file')


def fill_core_profiles(cv):

    cp = imas.core_profiles()
    cp.ids_properties.homogeneous_time = 1   # same timebase for core_profiles IDS
    cp.ids_properties.creation_date = datetime.datetime.today().strftime("%d/%m/%y")

    cp.time = np.atleast_1d(cv['TIME'].data)
    nt_cp = len(cp.time)

#    cp.vacuum_toroidal_field.r0 = cv['RMAJ'].data
    cp.vacuum_toroidal_field.b0 = cv['BTOR'].data

# profiles_1d

    cp1d = cp.profiles_1d
    cp1d.resize(nt_cp)

    nimp_flag = [0, 0, 0]
    for jimp in range(3):
        if np.max(cv['NIZ%d' %(jimp+1)].data) > 1e-8:
            nimp_flag[jimp] = 1
    n_ions_used = np.sum(nimp_flag) + 1
    n_ions_max = 4

    print(cv['ZMAIN'].shape, cv['ZIM1'].shape)
    for jt in range(nt_cp):
        qprof = 1./cv['MU'][jt, :]
        ne = e19m3_to_m3*cv['NE'][jt, :]
        ni = e19m3_to_m3*cv['NI'][jt, :]
        nd = e19m3_to_m3*cv['NDEUT'][jt, :]
        te = keV_to_eV*cv['TE'][jt, :]
        ti = keV_to_eV*cv['TI'][jt, :]
        cp1d[jt].grid.rho_tor_norm = cv['XRHO'].data
        cp1d[jt].grid.volume       = cv['VOLUM'][jt, :]
        cp1d[jt].grid.psi          = cv['FP']   [jt, :]
        cp1d[jt].electrons.temperature     = te
        cp1d[jt].electrons.density         = ne
        cp1d[jt].electrons.density_thermal = ne
        cp1d[jt].electrons.pressure         = keV_m3_to_Pa*ne*te
        cp1d[jt].electrons.pressure_thermal = keV_m3_to_Pa*ne*te
        cp1d[jt].q = 1./cv['MU'][jt, :]
        cp1d[jt].zeff = cv['ZEF'][jt, :]
        cp1d[jt].t_i_average = ti
        
        cp1d[jt].ion.resize(n_ions_used)
        cp1d[jt].ion[0].z_ion = np.float64(np.around(cv['ZMAIN'][jt, 0]))
        cp1d[jt].ion[0].density         = nd + e19m3_to_m3*cv['NIBM'][jt, :]
        cp1d[jt].ion[0].density_thermal = nd
        cp1d[jt].ion[0].density_fast    = e19m3_to_m3*cv['NIBM'][jt, :]
        cp1d[jt].ion[0].temperature     = ti
        cp1d[jt].ion[0].pressure_fast_perpendicular = cv['PBPER'][jt, :]
        cp1d[jt].ion[0].pressure_fast_parallel      = cv['PBLON'][jt, :]
        cp1d[jt].ion[0].element.resize(1)
        cp1d[jt].ion[0].element[0].a   = np.mean(cv['AMAIN'][jt, :])
        cp1d[jt].ion[0].element[0].z_n = np.mean(np.around(cv['ZMAIN'][jt, :]))
        jimp = 1
        for jion in range(n_ions_max):
            if jion > 0:
                if nimp_flag[jimp-1]:
                    cp1d[jt].ion[jimp].temperature = ti
                    cp1d[jt].ion[jimp].z_ion = np.float64(np.around(cv['ZIM%d' %jimp][jt, 0]))
                    cp1d[jt].ion[jimp].density         = e19m3_to_m3*cv['NIZ%d' %jimp][jt, :]
                    cp1d[jt].ion[jimp].density_thermal = e19m3_to_m3*cv['NIZ%d' %jimp][jt, :]
                    jimp += 1

    return cp


def fill_equilibrium(cv):

    eq = imas.equilibrium()
    eq.code.name = "astra"
    eq.code.version = "2024.06.01"
    eq.ids_properties.homogeneous_time = 1   # same timebase for core_profiles IDS
    eq.ids_properties.creation_date = datetime.datetime.today().strftime("%d/%m/%y")

    eq.time = np.atleast_1d(cv['TIME'].data)
    nt_eq = len(eq.time)
    logger.debug('%d', nt_eq)
    print(cv['TIME'].data)

    eq.vacuum_toroidal_field.r0 = cv['RTOR'].data[0]
    eq.vacuum_toroidal_field.b0 = cv['BTOR'].data
    
    prof_map = {'psi': 'psi', 'phi': 'phi', \
        'pressure': 'pressure', 'dpressure_dpsi': 'pprime', \
        'f_df_dpsi': 'ffprime', 'q': 'q', \
        'volume': 'volume', 'area': 'areat'}

    eqt = eq.time_slice
    eqt.resize(nt_eq)

    for itim in range(nt_eq):
 
        eqt[itim].boundary.outline.r = cv['r'][itim, -1, :]
        eqt[itim].boundary.outline.z = cv['z'][itim, -1, :]
        eqt[itim].boundary_separatrix.outline.r = cv['r'][itim, -1, :]
        eqt[itim].boundary_separatrix.outline.z = cv['z'][itim, -1, :]
        eqt[itim].profiles_1d.rho_tor_norm = cv['RHO_SURF'].data
        for imas_lbl, sf_lbl in prof_map.items():
            eqt[itim].profiles_1d.__dict__[imas_lbl] = \
                np.array(cv[sf_lbl][itim])

# global quantities
# this all assumes standard AUG shotfile writing starting from separatrix and moving towards axis
# This should be checked in future, as, I believe, IDE does not do this...
        eqt[itim].global_quantities.volume       = cv['volume'][itim, -1]
        eqt[itim].global_quantities.area         = cv['areat' ][itim, -1]
        eqt[itim].global_quantities.psi_axis     = cv['psiaxis'][itim]
        eqt[itim].global_quantities.psi_boundary = cv['psibound'][itim]
        eqt[itim].global_quantities.ip           = cv['i_plasma'][itim]
        eqt[itim].global_quantities.magnetic_axis.r = cv['r'][itim, 1, 1]
        eqt[itim].global_quantities.magnetic_axis.z = cv['z'][itim, 1, 1]
        eqt[itim].global_quantities.magnetic_axis.b_field_tor = cv['b0'][itim]

    return eq

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='CDF to IMAS conversion')
    parser.add_argument('-f', '--fcdf', help='ASTRA NetCDF output (full path)', required=False, default='aug34954fluxes.CDF')
    parser.add_argument('-i', '--ids_backend', help='IDS backend: HDF5, MDSPLUS or ASCII', required=False, default='HDF5')
    parser.add_argument('-s', '--shot'   , type=int, help='Shot number', required=False, default=34954)
    parser.add_argument('-r', '--ids_run', type=int, help='IDS run'    , required=False, default=6)
    args = parser.parse_args()

    ACDF2IMAS(args)
