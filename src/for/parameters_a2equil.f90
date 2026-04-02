module parameters_a2equil       ! declaration of code parameters

use imas_ids, only: type_equilibrium

integer :: fix_adapgrid, s_fazt=0
integer, parameter, private :: DP=kind(1.0D0)
double precision, parameter :: GP=3.14159265359, GP2=2.*GP, GP4=GP2**2.0, muvac=0.4*GP*1.E-06
double precision :: epsf_tol, epss_tol, ydiff, ydiff2, epsv_tol, sorparam, epsg_tol, &
    epstol, itertol, lambdatol, urelax, urelax2, time_fix_eqpff, murelax2, cheb_coefs(5, 5)
character(len=80) :: name_gsefdir = 'exp/equ/' ! working directory path

integer :: max_iter, miter_ext, diagnostic_gsef, interp_routine, &
   do_adcmp, interp_method_rect, cheb_degree, four_degree, &
   advanced_methods, iter_one_only_fbe, spidat_yes, i3method, &
   fix_eqpf_eqff
integer :: key_no_startz  ! if 0, if key_start=1, skip iterations
integer :: key_no_refits  ! if 0, does not overwrite coil.dat with new refit currents
type(type_equilibrium) :: equil_now

!feqis
double precision :: err_circ_in                 = 1.e-9  ! err circ
double precision :: err_find_oxpoints_in        = 1.e-13 ! err find oxpoints
double precision :: err_find_oxpoints_derivs_in = 1.e-13 ! err find oxpoints deriv
double precision :: err_find_psistab_in         = 1.e-8  ! err find psistab
double precision :: err_find_delr_in            = 1.e-10 ! err find delr
double precision :: err_find_biquad_in          = 1.e-12 ! err find biquad
double precision :: err_epsilon_in              = 1.e-12 ! epsilon
double precision :: err_gaptolez_in             = 1.e-5  ! err gap tolez
double precision :: err_fix_boundary_in         = 1.e-9  ! fix boundary tolerance

end module parameters_a2equil
