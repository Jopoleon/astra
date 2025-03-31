module pbe_core

implicit none

!grids
double precision :: dteta, psibndp, psiaxisp
double precision, dimension(:), allocatable :: teta, psigrid
double precision, dimension(:, :), allocatable :: rho, &
    rpol, zpol, rpul, zpul, psirhoteta

! r(z)pol: R, Z in polar coordinates half radial grid
! r(z)pul: R, Z in polar coordinates full radial grid
double precision, dimension(:), allocatable :: psia_1d, ffp_1d, ppp_1d

! boundary and axis FBE, PBE
double precision :: raxp, zaxp
double precision, dimension(:), allocatable :: green_bnd_f
double precision, dimension(:), allocatable :: rbndp, zbndp, rexp, zexp, tetaexp

! plasma parameters
double precision, dimension(:), allocatable :: pprime, ffprime, pressure, psigrida, ipol
double precision, dimension(:, :), allocatable :: jrhoteta

integer :: nrho, nteta

end module pbe_core
