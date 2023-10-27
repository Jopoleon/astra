!---------------------------------------------------------------------
module feqis_circuit

use feqis_dimensions, only: i_dim1, i_dim2, i_dim5, nrho2d

implicit none 

integer :: use_limiter_yesno

!time stepping
double precision :: tau_old, tau_new
double precision, dimension(i_dim1) :: psi_cur_old, dpc

!circuits
integer :: ncoils, nreseqcoil, nlimiter
integer, dimension(i_dim1) :: mequivalence
double precision, dimension(500) :: limiterR, limiterZ
double precision :: lim_maxR, lim_minR, lim_maxZ, lim_minZ
integer :: ilim_maxR, ilim_minR, ilim_maxZ, ilim_minZ
double precision, dimension(i_dim1) :: Rcoil, Zcoil, drcoil, dzcoil, &
    anglecoil

integer :: nconduc, nblocks, npassive, nactive
double precision, dimension(i_dim1) :: curconduc, voltage, voltage_old, &
    cur_con_old, r_cond, z_cond
double precision, dimension(i_dim1, i_dim1) :: resconduc, indconduc
double precision :: psiplasmatoconduc(i_dim1) !plasma --> conduc at t

! coordinates:
! r, z --> rectangular grid in meters
! rho --> distance of point from magnetic axis in meters
! teta --> angle from LFS midplane which is 0, coincides with Zmag plane
! psigrid --> in poloidal flux equispaced

!grids
integer :: nr, nz, nrho, nteta, nr2, nz2, nr1, nz1, nbnd, ngbnd, &
    redo_bnd

double precision :: rmin, rmax, zmin, zmax, dr, dz, dteta, &
    zbot, ztop, raus, rinner
double precision, dimension(i_dim2) :: r, z, teta, rcomp, zcomp, psigrid  
double precision, dimension(i_dim2, i_dim2) :: rho, area_eff, &
    rpol, zpol, rpul, zpul, u_n, omega_pl, zlimpotential, &
    psirz, psirhoteta, psiextrz, psiplasrz
! r(z)pol: R, Z in polar coordinates half radial grid
! r(z)pul: R, Z in polar coordinates full radial grid
double precision, dimension(nrho2d) :: psia_2d, ffp_2d, ppp_2d
double precision, dimension(i_dim2) :: psia_1d, ffp_1d, ppp_1d
double precision, dimension(8) :: derivpsi

! boundary and axis FBE, PBE
integer, parameter :: max_xpoints=500
integer :: i_plasmatype !(0-limited, 1-single null, 2-double null)
integer :: iaxis, jaxis, n_of_xpoints
double precision :: psibnd, psiaxis, psibndp, psiaxisp, &
    rax, zax, trax, tzax, raxp, zaxp, &
    alpsep, psistabR, psistabZ, dr_factor_init, dz_factor_init
double precision, dimension(i_dim5) :: rbnd, zbnd
double precision, dimension(max_xpoints) :: r_xpoint, z_xpoint
double precision :: deriv_x(5, max_xpoints)
double precision :: green_bnd_f(16*i_dim2**2)
double precision, dimension(i_dim2) :: rbndp, zbndp, rexp, zexp, tetaexp

! plasma parameters
double precision :: iplasma, btor0, rgeom0, psplex, li3, betapol
double precision, dimension(i_dim2) :: pprime, ffprime, pressure, psigrida, ipol
double precision, dimension(i_dim2, i_dim2) :: jrz, jrhoteta

end module feqis_circuit
