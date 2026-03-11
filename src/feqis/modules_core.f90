module errors_params

implicit none

double precision, parameter :: err_circ_plasma_iter     = 1.e-9  ! err circ
double precision, parameter :: err_find_oxpoints        = 1.e-13 ! err find oxpoints
double precision, parameter :: err_find_oxpoints_derivs = 1.e-13 ! err find oxpoints deriv
double precision, parameter :: err_find_psistab         = 1.e-8  ! err find psistab
double precision, parameter :: err_find_delr            = 1.e-10 ! err find delr
double precision, parameter :: err_find_biquad          = 1.e-12 ! err find biquad
double precision, parameter :: err_epsilon              = 1.e-12 ! epsilon
double precision, parameter :: err_gaptolez             = 1.e-5  ! err gap tolez
double precision, parameter :: err_fix_boundary         = 1.e-9  ! fix boundary tolerance

end module errors_params

!---------------------------------------------------------------------
module green_function

implicit none

double precision, dimension(:, :, :), allocatable :: greeni
end module green_function

!---------------------------------------------------------------------
module pi_vars

implicit none

double precision, parameter :: GPI=3.141592653589793, GPI2=2.*GPI, &
    GPI4=GPI2**2, muvac=4.e-7*GPI, mu0=0.4*GPI

end module pi_vars
