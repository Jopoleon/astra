module parameters_a2equil       ! declaration of code parameters

use imas_ids, only: type_equilibrium

integer :: fix_adapgrid, s_fazt=0
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

end module parameters_a2equil
