module ferromagstructure

implicit none

type type_position   ! Structure for list of R,Z positions (1D)
    integer :: npoints
    integer :: sigma_surface
    double precision, pointer :: R(:)       => null()  ! /r - Major radius [m]. Vector(max_npoints). Time-dependent
    double precision, pointer :: Z(:)       => null()  ! /z - Altitude [m]. Vector(max_npoints). Time-dependent
    double precision, pointer :: tanangl(:) => null()  ! /z - Altitude [m]. Vector(max_npoints). Time-dependent
    double precision, pointer :: length(:)  => null()  ! /z - Altitude [m]. Vector(max_npoints). Time-dependent
    double precision, pointer :: MagnetizationChi(:)  => null()  ! /calculated chi for this element
    double precision, pointer :: Btangfield(:)  => null()  ! /calculated chi for this element
    double precision, pointer :: Current(:)  => null()  ! /I_s of this element in MA
endtype type_position
 
type type_mutmatrix  ! Structure for list of R,Z positions (1D)
    double precision, pointer :: Mij(:, :) => null()    ! /mutual induction matrix
endtype type_mutmatrix
 
type type_magnetiz   ! Structure for list of R,Z positions (1D)
    integer :: nvalues
    double precision, pointer :: Chi(:) => null()   ! /magnetic suceptibility
    double precision, pointer :: H(:)   => null()   ! /Bvacuum_tangent/mu0
endtype

type type_ferromag
    type (type_position)  :: position      ! /ferromag/MHrelation - RZ description 
    type (type_magnetiz)  :: MHrelation    ! /ferromag/MHrelation - 
    type (type_mutmatrix) :: mutual_matrix ! /ferromag/mutual_matrix - 
endtype

end module ferromagstructure

!---------------------------------------------------------------------
module feqis_scalars

implicit none

double precision :: psplex, li3, li_aug, betapol, betapol_iter, &
    wkin, bpkin, iplasma, btor0, rgeom0

end module feqis_scalars

!---------------------------------------------------------------------
module transfer_functions

implicit none

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

implicit none

! metric coefficients in polar coordinates 
double precision :: R_curr_0D, Z_curr_0D
double precision, dimension(:, :), allocatable :: lambda2d, dator, fsa_kernel

end module metric_coefficients_pbe

!---------------------------------------------------------------------
module green_matrix

implicit none

double precision, dimension(:, :, :), allocatable :: dgreenirpl, dgreenizpl
double precision, dimension(:, :), allocatable :: dgreenirj, dgreenizj

end module green_matrix
