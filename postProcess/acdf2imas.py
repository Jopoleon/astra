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
MW_to_W = 1e6
MA_to_A = 1e6

spec_lbl = {(1, 1): 'H', (1, 2): 'D', (2, 4): 'He', (6, 12): 'C'}


def ACDF2IMAS(args, db='aug', v_major='3'):

    logger.info('Creating IDS structure')

    if args.ids_backend == 'HDF5':
        db = imas.DBEntry(imas.imasdef.HDF5_BACKEND   , db, args.shot, args.ids_run, os.getenv('IMASDB'), v_major)
    elif args.ids_backend == 'MDS+':
        db = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND, db, args.shot, args.ids_run, os.getenv('IMASDB'), v_major)
    elif args.ids_backend == 'ASCII':
        db = imas.DBEntry(imas.imasdef.ASCII_BACKEND  , db, args.shot, args.ids_run, os.getenv('IMASDB'), v_major)
    status, _ = db.create()

    cv = netcdf_file(args.fcdf, 'r', mmap=False).variables
    cp = fill_core_profiles(cv)
    eq = fill_equilibrium(cv)
    cs = fill_core_sources(cv)

    logger.info('Dumping IDS file %s', args.ids_backend)
    logger.info('in dir %s' %os.getenv('IMASDB'))
    logger.info('Putting core_profiles')
    db.put(cp)
    logger.info('Putting equilibrium')
    db.put(eq)
    logger.info('Putting core_sources')
    db.put(cs)
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
        cp1d[jt].j_total = MA_to_A*cv['CU'][jt, :]
        cp1d[jt].j_bootstrap = MA_to_A*cv['CUBS'][jt, :]
        
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
    eq.code.version = "2026.03.01"
    eq.ids_properties.homogeneous_time = 1   # same timebase for core_profiles IDS
    eq.ids_properties.creation_date = datetime.datetime.today().strftime("%d/%m/%y")

    eq.time = np.atleast_1d(cv['TIME'].data)
    nt_eq = len(eq.time)
    logger.debug('%d', nt_eq)

    eq.vacuum_toroidal_field.r0 = cv['RTOR'].data[0]
    eq.vacuum_toroidal_field.b0 = cv['BTOR'].data

    imas_profs = ['psi', 'phi', 'pressure', 'dpressure_dpsi', 'f_df_dpsi', 'q', \
                  'volume', 'area', 'r_inboard', 'r_outboard']
    prof_map = {'dpressure_dpsi': 'pprime', 'f_df_dpsi': 'ffprime', 'area': 'areat'}

    eqt = eq.time_slice
    eqt.resize(nt_eq)

    for itim in range(nt_eq):
 
        eqt[itim].boundary.outline.r = cv['r'][itim, -1, :]
        eqt[itim].boundary.outline.z = cv['z'][itim, -1, :]
        eqt[itim].boundary_separatrix.outline.r = cv['r'][itim, -1, :]
        eqt[itim].boundary_separatrix.outline.z = cv['z'][itim, -1, :]
        eqt[itim].profiles_1d.rho_tor_norm = cv['RHO_SURF'].data
        for imas_lbl in imas_profs:
            if imas_lbl in prof_map:
                sf_lbl = prof_map[imas_lbl]
            else:
                sf_lbl = imas_lbl
            eqt[itim].profiles_1d.__dict__[imas_lbl] = \
                np.array(cv[sf_lbl][itim])
        eqt[itim].profiles_1d.phi = cv['RHO_SURF'].data**2 * np.pi * cv['BTOR'][itim]*cv['ROC'][itim]

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


def fill_core_sources(cv):

    from identifiers.core_source_identifier import core_source_identifier

    cs = imas.core_sources()
    cs.ids_properties.creation_date = datetime.datetime.today().strftime("%d/%m/%y")
    cs.ids_properties.homogeneous_time = 1
    nsrc = 2
    cs.time = np.atleast_1d(cv['TIME'].data)
    nt = len(cs.time)

    cs.source.resize(nsrc)

#-----#
# NBI #
#-----#

    jsrc = 0
    cs.source[jsrc].identifier.name = 'nbi'
    cs.source[jsrc].identifier.index = core_source_identifier['nbi']['index']
    cs.source[jsrc].identifier.description = core_source_identifier['nbi']['description']
    cs.source[jsrc].code.name = 'ASTRA/RABBIT'
    cs_nb = cs.source[jsrc].profiles_1d
    cs_nb.resize(nt)
    spec_label = spec_lbl[(1, 2)] # Hardcoded to Deuterium

    for jt in range(nt):
# Filling
        cs_nb[jt].grid.rho_tor_norm = cv['XRHO'].data
        cs_nb[jt].grid.volume       = cv['VOLUM'][jt, :]
        cs_nb[jt].grid.area         = cv['AREAT'][jt, :]
        #Rabbit provides toroidal current, not parallel current
        cs_nb[jt].j_parallel          = MA_to_A*cv['CUBM'][jt, :]
        cs_nb[jt].momentum_tor        = cv['RTOR'][jt]*cv['SCUBM'][jt, :]
        cs_nb[jt].electrons.particles = e19m3_to_m3*cv['SNEBM'][jt, :]
        cs_nb[jt].total_ion_energy    = MW_to_W*cv['PIBM'][jt, :]
        cs_nb[jt].electrons.energy    = MW_to_W*cv['PEBM'][jt, :]
        cs_nb[jt].ion.resize(1)
        cs_nb[jt].ion[0].element.resize(1)
        cs_nb[jt].ion[0].element[0].a   = 1.0
        cs_nb[jt].ion[0].element[0].z_n = 2.0
        cs_nb[jt].ion[0].element[0].atoms_n = 1 # Number of atoms of this element in the molecule
        cs_nb[jt].ion[0].z_ion     = 2.
        cs_nb[jt].ion[0].label     = spec_label
        cs_nb[jt].ion[0].particles = e19m3_to_m3*cv['SNEBM'][jt, :]
        cs_nb[jt].ion[0].energy    = MW_to_W*cv['PIBM'][jt, :]

#------#
# ECRH #
#------#

    jsrc = 1
    cs.source[jsrc].identifier.name = 'ec'
    cs.source[jsrc].identifier.index = core_source_identifier['ec']['index']
    cs.source[jsrc].identifier.description = core_source_identifier['ec']['description']
    cs.source[jsrc].code.name = 'ASTRA/TORBEAM'
    cs_ec = cs.source[jsrc].profiles_1d
    cs_ec.resize(nt)

    for jt in range(nt):
# Filling
        cs_ec[jt].grid.rho_tor_norm = cv['XRHO'].data
        cs_ec[jt].grid.volume       = cv['VOLUM'][jt, :]
        cs_ec[jt].grid.area         = cv['AREAT'][jt, :]
        cs_ec[jt].electrons.energy  = MW_to_W*cv['PEECR'][jt, :]
        cs_ec[jt].j_parallel        = MA_to_A*cv['CUECR'][jt, :]

    return cs

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='CDF to IMAS conversion')
    parser.add_argument('-f' , '--fcdf', help='ASTRA NetCDF output (full path)', required=False, default='../ncdf_out/aug34954fluxes.CDF')
    parser.add_argument('-i' , '--ids_backend', help='IDS backend: HDF5, MDSPLUS or ASCII', required=False, default='HDF5')
    parser.add_argument('-s' , '--shot', type=int, help='Shot number', required=False, default=34954)
    parser.add_argument('-r' , '--ids_run', type=int, help='IDS run', required=False, default=6)
    parser.add_argument('-db', '--database', help='IDS DataBase', required=False, default='aug')
    args = parser.parse_args()

    ACDF2IMAS(args)
