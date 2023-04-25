module circuit

use dimensions_ef_parameters, only: i_dim1, i_dim2, i_dim5

implicit none 

character(len=120) :: data_dir
integer :: data_dir_k, iteration_step, use_limiter_yesno, &
    nr_of_fit_parameters, ncoils, nsubcoils, nreseqcoil, &
    nfirstwall, nlimiter, &
    nblanket, nblanketpc, nelemblanketpc, nelemblanket, &
    ilim_maxR, ilim_minR, ilim_maxZ, ilim_minZ, &
    nconduc, nblocks, npassive, nactive
integer, dimension(i_dim1) :: nelemcoil, mturns, mequivalence
double precision :: rmag_fit, zmag_fit, k_fit, rxp_fit, zxp_fit, &
    tau_old, tau_new, &
    lim_maxR, lim_minR, lim_maxZ, lim_minZ, &
    reswall, resblan, widthblan, psiconductoplasma
double precision, dimension(i_dim1) :: psi_cur_old, dpc, &
    Rcoil, Zcoil, drcoil, dzcoil, anglecoil, indcoil, curcoil, &
    Rwall, Zwall, areawall, curwall, Rblan, Zblan, areablan, curblan, &
    Rblanpc, Zblanpc, resblanpc, areablanpc, curblanpc, indblanpc, &
    curconduc, voltage, voltage_old, psiplasmatoconduc, &
    cur_con_old, r_cond, z_cond
double precision, dimension(300) :: indwall, indblan
double precision, dimension(500) :: limiterR, limiterZ
double precision, dimension(i_dim1, i_dim1) :: rescoil, resconduc, indconduc

! coordinates:
! r,z --> rectangular grid in meters
! rho --> distance of point from magnetic axis in meters
! teta --> angle from LFS midplane which is 0, coincides with Zmag plane
! psigrid --> in poloidal flux equispaced

integer :: nr, nz, nrho, nteta, nr2, nz2, nr1, nz1  ! nteta+1 is the periodic point. nrho is the plasma boundary

double precision :: rmin, rmax, dr, zmin, zmax, dz, dteta, &
    zbot, ztop, raus, rinner, &
    iplasma, btor0, rgeom0, psplex, li3, betapol
double precision, dimension(8) :: derivpsi
double precision, dimension(i_dim2) :: r, z, teta, rcomp, zcomp, &
    psigrid, psia_2d, ffp_2d, ppp_2d, ipol_2d, pres_2d, &
    psia_1d, ffp_1d, ppp_1d, &
    pprime, ffprime, pressure, psigrida, ipol, &
    rbndp, zbndp, rexp, zexp, tetaexp
double precision, dimension(i_dim2, i_dim2) :: rho, drho, rpol, zpol, &
    rpul, zpul, area_eff, psirz, psirhoteta, u_n, omega_pl, &
    psiextrz, psiplasrz, phirhoteta, zlimpotential, &
    jrz, jrhoteta

! boundary and axis FBE
integer :: nbnd, ngbnd, redo_bnd, i_plasmatype, iaxis, jaxis, &
    n_of_xpoints, max_xpoints=50
integer, dimension(i_dim5) :: ibnd, jbnd
double precision :: alpsep, psibnd, psiaxis, &
    psibndp, psiaxisp, psistabR, psistabZ, &
    rax, raxp, trax, dr_factor_init, &
    zax, zaxp, tzax, dz_factor_init
double precision, dimension(50) :: r_xpoint, z_xpoint
double precision, dimension(i_dim5) :: rbnd, zbnd
double precision, dimension(5, 50) :: deriv_x
double precision :: green_bnd_f(16*i_dim2**2) ! Memory!!

contains

!---------------------------------------------------------------------
    subroutine interp_j_fromrhotorz

    integer :: i, j

    jrz = 0.
    do j=1, nz2
        do i=1, nr2
            call curinterp_ef(r(i), z(j), jrhoteta(1:nrho-1, 1:nteta), & 
                rho(1:nrho-1, 1:nteta), teta(1:nteta), raxp, zaxp, nrho-1, nteta, jrz(i, j))
        enddo
    enddo

    return
    end subroutine interp_j_fromrhotorz

!---------------------------------------------------------------------
    subroutine plasma_psi_to_coils_ef

    use green_matrix, only: greeni

    implicit none

    integer :: i

    do i=1, nconduc
        psiplasmatoconduc(i) = sum(jrz(1:nr2, 1:nz2) * &
            area_eff(1:nr2, 1:nz2)*greeni(1:nr2, 1:nz2, i))
    enddo

    return
    end subroutine plasma_psi_to_coils_ef

!---------------------------------------------------------------------
    subroutine psi_external_calc_ef

    use green_matrix, only: greeni

    implicit none

    integer :: i, j

    do j=1, nz2
        do i=1, nr2
            psiextrz(i, j) = sum(curconduc(1:nconduc)*greeni(i, j, 1:nconduc))
        enddo
    enddo

    return
    end subroutine psi_external_calc_ef

!---------------------------------------------------------------------
    subroutine green_bnd_integral(bgintsol, g, jcounty)

! calculates  integral_over_boundary of -Green * dg/dn * dl for point r0, z0
    use dimensions_ef_parameters, only: i_dim2

    implicit none

    double precision, intent(in), dimension(i_dim2, i_dim2) :: g
    integer, intent(out) :: jcounty
    double precision, intent(out) :: bgintsol

    integer :: j
    double precision :: greenf
    double precision, dimension(i_dim2) :: dgdn

    bgintsol = 0.
! lower side
    do j=1, nr1
        jcounty = jcounty + 1
        greenf = green_bnd_f(jcounty)
        dgdn(j) = -(g(j, 2) - g(j, 1) + g(j+1, 2) - g(j+1, 1))/2./dz*greenf*dr/(r(j) + dr/2.)
    enddo
    bgintsol = bgintsol - sum(dgdn(1:nr1))
! right side
    do j=1, nz1
        jcounty = jcounty + 1
        greenf = green_bnd_f(jcounty)
        dgdn(j) = (g(nr2, j) - g(nr1, j) + g(nr2, j+1) - g(nr1, j+1))/2./dr*greenf*dz/(r(nr2) + r(nr1))*2.
    enddo
    bgintsol = bgintsol - sum(dgdn(1:nz1))
! upper side
    do j=1, nr1
        jcounty = jcounty + 1
        greenf = green_bnd_f(jcounty)
        dgdn(j) = (g(j, nz2) - g(j, nz1) + g(j+1, nz2) - g(j+1, nz1))/2./dz*greenf*dr/(r(j) + dr/2.)
    enddo
    bgintsol = bgintsol - sum(dgdn(1:nr1))
! left side
    do j=1, nz1
        jcounty = jcounty + 1
        greenf = green_bnd_f(jcounty)
        dgdn(j) = -(g(2, j) - g(1, j) + g(2, j+1) - g(1, j+1))/2./dr*greenf*dz/(r(1) + r(2))*2.
    enddo
    bgintsol = bgintsol - sum(dgdn(1:nz1))

    return
    end subroutine green_bnd_integral

end module circuit
