import os, argparse
from trview import plasma_state

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
astra_ext = os.getenv('ASTRA_EXT')

aug_nml = \
'''&species
   Aimp = 200.00
   Zimp = 100.00
/

&output_settings
   it_orbout = -1
/

&physics
   jumpcor = 2
   table_path   = '%s/rabbit/tables_highRes'
   limiter_file = '%s/rabbit/limiters/limiter_aug.dat'
   Rlim = 2.20 !AUG defaults, overridden if limiter_file exists
   zlim = 0.20
   Rmax = 2.26 !AUG defaults
   Rmin = 1.08
   torqjxb_model = 3
/

&numerics
   norbits = 20
   orbit_dt_fac = 1.56
   distfun_nv = 200
   distfun_vmax = 4.5e6
/
''' %(astra_ext, astra_ext)


class PROF:
    pass


def sf2astra(shot, tbeg, tend):

    gpar = type('', (), {})()
    gpar.shot    = shot
    gpar.tbeg    = tbeg
    gpar.tend    = tend
    gpar.dt_in   = 0.05
    gpar.eq_exp  = 'AUGD'
    gpar.eq_diag = 'EQH'
    gpar.eq_ed   = 0
    gpar.nrho_in = 51
    gpar.nthe_in = 91
    gpar.nmom    = 6
    gpar.fit_tol = 0.3
    gpar.sigmas  = 2.5
    gpar.pel_trig= 0.25
    gpar.nml     = aug_nml
    gpar.angf    = True
    gpar.no_elm  = False

    sflists = {'Ne': 'AUGD:IDA:ne:0', 'Te': 'AUGD:IDA:Te:0', \
               'Ti': 'AUGD:CEZ:Ti_c:0', 'Angf': 'AUGD:CEZ:vrot:0 AUGD:CMZ:vrot:0'}
    gprof = {}
    for key in sflists.keys():
        gprof[key] = PROF()
        gprof[key].method = 'Rec. spline'
        gprof[key].sigma  = 1.
        gprof[key].sflist = sflists[key]
        if key in ('Ne', 'Te'):
            gprof[key].fit_tol = 0.
        else:
            gprof[key].fit_tol = 0.3
    plasma = plasma_state.PLASMA_STATE(shot)
    plasma.fromShotfile(gpar, gprof)
    plasma.toASTRA(awd=awd)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Reading IMAS output and writing ASTRA input')
    parser.add_argument('-s', '--shot', type=int, help='Shot number', required=False, default=38384)
    parser.add_argument('-tb', '--tbeg', type=float, help='t beg [s]', required=False, default=0.)
    parser.add_argument('-te', '--tend', type=float, help='t end [s]', required=False, default=10.)
    args = parser.parse_args()

    sf2astra(args.shot, args.tbeg, args.tend)
