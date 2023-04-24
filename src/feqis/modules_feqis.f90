module dimensions_ef_parameters

implicit none

integer, parameter :: i_dim1=300, i_dim2=300 ! #coils, #rad.grid

end module dimensions_ef_parameters

!-------------------------------------------
module pi_vars

implicit none

double precision, parameter :: GPI=3.141592653589793, &
    GPI2=2.*GPI, mu0=0.4*GPI

end module pi_vars

!-------------------------------------------
module ef_circuit       ! declaration of minimal CPOs

use dimensions_ef_parameters, only: i_dim2

implicit none 

!generic
character(len=120) :: data_dir

! coordinates:
! r, z --> rectangular grid in meters
! rho --> distance of point from magnetic axis in meters
! teta --> angle from LFS midplane which is 0, coincides with Zmag plane
! psigrid --> in poloidal flux equispaced

!grids
integer :: nr, nz, nrho, nteta, nr2, nz2, ncoils

double precision :: rmin, rmax, zmin, zmax, psibndp, psiaxisp, raxp, zaxp, &
    iplasma, btor0, rgeom0, psplex, li3, betapol
double precision, dimension(i_dim1) :: voltage, psiplasmatoconduc, psi_cur_old
double precision, dimension(i_dim2) :: teta, rbndp, zbndp, rexp, zexp, &
    psia_2d, ffp_2d, ppp_2d, ipol_2d, pres_2d, pprime, ffprime, pressure, psigrida, ipol
double precision, dimension(i_dim2, i_dim2) :: rho, jrhoteta, psiextrz, psirz

end module ef_circuit
