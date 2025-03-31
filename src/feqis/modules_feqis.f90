module errors_params

double precision :: err_circ_plasma_iter, err_find_oxpoints, &
    err_find_oxpoints_derivs, err_find_psistab, err_find_delr, &
    err_find_biquad, err_epsilon, err_gaptolez, err_fix_boundary

end module errors_params

!---------------------------------------------------------------------
module ferromagstructure

integer, parameter, private :: DP=kind(1.0D0)

type type_position   ! Structure for list of R,Z positions (1D)
    integer :: npoints
    integer :: sigma_surface
    real(DP), pointer :: R(:)       => null()  ! /r - Major radius [m]. Vector(max_npoints). Time-dependent
    real(DP), pointer :: Z(:)       => null()  ! /z - Altitude [m]. Vector(max_npoints). Time-dependent
    real(DP), pointer :: tanangl(:) => null()  ! /z - Altitude [m]. Vector(max_npoints). Time-dependent
    real(DP), pointer :: length(:)  => null()  ! /z - Altitude [m]. Vector(max_npoints). Time-dependent
    real(DP), pointer :: MagnetizationChi(:)  => null()  ! /calculated chi for this element
    real(DP), pointer :: Btangfield(:)  => null()  ! /calculated chi for this element
    real(DP), pointer :: Current(:)  => null()  ! /I_s of this element in MA
endtype type_position
 
type type_mutmatrix  ! Structure for list of R,Z positions (1D)
    real(DP), pointer :: Mij(:, :) => null()    ! /mutual induction matrix
endtype type_mutmatrix
 
type type_magnetiz   ! Structure for list of R,Z positions (1D)
    integer :: nvalues
    real(DP), pointer :: Chi(:) => null()   ! /magnetic suceptibility
    real(DP), pointer :: H(:)   => null()   ! /Bvacuum_tangent/mu0
endtype

type type_ferromag
    type (type_position)  :: position      ! /ferromag/MHrelation - RZ description 
    type (type_magnetiz)  :: MHrelation    ! /ferromag/MHrelation - 
    type (type_mutmatrix) :: mutual_matrix ! /ferromag/mutual_matrix - 
endtype

end module ferromagstructure

!---------------------------------------------------------------------
module global_params

double precision :: iplasma, btor0, rgeom0, psplex, li3, li_aug, betapol, betapol_iter, &
 wkin, bpkin

end module global_params

!---------------------------------------------------------------------
module transfer_functions

double precision, dimension(:), allocatable :: dpsidvbez, psibez, &
    g2bez, g2ibez, gm1bez, routbez, rinbez, vbez, g1bez, gm41bez, &
    ggrhobez, bmaxbez, bminbez, gm4bez, bdb0bez, gm5bez, fofbbez, &
    areatbez, perimbez, shifbez, triaubez, trialbez, kbez, surfbez, phibez, &
    qbez, t2dbez, rbp2_b2bez, &
    ffprimebez, pprimebez, pressbez, ipolbez, shivbez, squarebez

double precision, dimension(:, :), allocatable :: rpbez, zpbez, &
    rminbez, bpcellbez, bcellbez, rmin2dbez, jrhobez, &
    bpcell2dbez, bcell2dbez

end module transfer_functions

!---------------------------------------------------------------------
module metric_coefficients_pbe

! metric coefficients in polar coordinates 
double precision :: R_curr_0D, Z_curr_0D
double precision, dimension(:, :), allocatable :: lambda2d, dator, fsa_kernel

end module metric_coefficients_pbe

!---------------------------------------------------------------------
module green_matrix

double precision, dimension(:, :, :), allocatable :: greeni, &
    dgreenirpl, dgreenizpl
double precision, dimension(:, :), allocatable :: dgreenirj, dgreenizj

end module green_matrix

!---------------------------------------------------------------------
module fft_mod_eff

use pi_vars, only: GPI2

implicit none

integer, parameter :: dp=selected_real_kind(15, 300)
double precision, dimension(:, :), allocatable :: sintable, costable

end module fft_mod_eff
