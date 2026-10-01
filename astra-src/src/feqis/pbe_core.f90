module pbe_core

implicit none

!grids
double precision :: dtheta, psibndp, psiaxisp
double precision, dimension(:), allocatable :: theta, psigrid
double precision, dimension(:, :), allocatable :: rho, &
    rpol, zpol, rpul, zpul, psirhotheta

! r(z)pol: R, Z in polar coordinates half radial grid
! r(z)pul: R, Z in polar coordinates full radial grid
double precision, dimension(:), allocatable :: psia_1d, ffp_1d, ppp_1d

! boundary and axis FBE, PBE
double precision :: raxp, zaxp
double precision, dimension(:), allocatable :: green_bnd_f
double precision, dimension(:), allocatable :: rbndp, zbndp, rexp, zexp, thetaexp

! plasma parameters
double precision, dimension(:), allocatable :: pprime, ffprime, pressure, psigrida, ipol
double precision, dimension(:, :), allocatable :: jrhotheta

integer :: nrho, ntheta

end module pbe_core
