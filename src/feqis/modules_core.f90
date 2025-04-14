module errors_params

implicit none

double precision :: err_circ_plasma_iter, err_find_oxpoints, &
    err_find_oxpoints_derivs, err_find_psistab, err_find_delr, &
    err_find_biquad, err_epsilon, err_gaptolez, err_fix_boundary

end module errors_params

!---------------------------------------------------------------------
module global_params

implicit none

double precision :: iplasma, btor0, rgeom0

end module global_params

!---------------------------------------------------------------------
module green_function

implicit none

double precision, dimension(:, :, :), allocatable :: greeni
end module green_function

!---------------------------------------------------------------------
module fft_mod_eff

implicit none

integer, parameter :: dp=selected_real_kind(15, 300)
double precision, dimension(:, :), allocatable :: costable
double precision, dimension(:, :), allocatable :: sintable

end module fft_mod_eff

!---------------------------------------------------------------------
module pi_vars

implicit none

double precision, parameter :: GPI=3.141592653589793, GPI2=2.*GPI, &
    GPI4=GPI2**2, muvac=4.e-7*GPI, mu0=0.4*GPI

end module pi_vars
