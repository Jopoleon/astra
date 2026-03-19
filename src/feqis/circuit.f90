module circuit

use fbe_core, only: nr, nr1, nr2, nz, nz1, nz2, nbnd, nlimiter, &
    nconduc, n_of_xpoints, nrho2d, &
    Rrect, Zrect, dr, dz, Rbnd, Zbnd, rax, zax, &
    iaxis, jaxis, r_xpoint, z_xpoint, limiterR, limiterZ, &
    jrz, curconduc, psia_2d, ffp_2d, ppp_2d, &
    psiextrz, psiplasrz, psirz, psiaxis, psibnd, psistabR, psistabZ, psiferro, &
    get_psiplasrz, find_new_axis, find_psi_boundary, new_jrz, &
    compound_psi, psi_external_calc
use pbe_core, only: raxp, zaxp, rbndp, zbndp, rpol, zpol, ntheta, nrho, &
    rho, jrhotheta, theta, pprime, ffprime, psigrida, psirhotheta
use feqis_scalars, only: iplasma

implicit none

!time stepping
double precision :: tau_old, tau_new
double precision, dimension(:), allocatable :: psi_cur_old, dpc

!circuits
integer :: ncoils, nreseqcoil
integer, dimension(:), allocatable :: mequivalence
double precision, dimension(:), allocatable :: Rcoil, Zcoil, drcoil, dzcoil, &
    anglecoil, anglehcoil

integer :: nblocks, npassive, nactive, nferromag
double precision, dimension(:), allocatable :: voltage, voltage_old, &
    cur_con_old, r_cond, z_cond
double precision, dimension(:, :), allocatable :: resconduc, indconduc
double precision, dimension(:), allocatable :: psiplasmatoconduc !plasma --> conduc at t

! coordinates:
! r, z --> rectangular grid in meters
! rho --> distance of point from magnetic axis in meters
! theta --> angle from LFS midplane which is 0, coincides with Zmag plane
! psigrid --> in poloidal flux equispaced

!grids

! r(z)pol: R, Z in polar coordinates half radial grid
! r(z)pul: R, Z in polar coordinates full radial grid

! boundary and axis FBE, PBE

! plasma parameters

!ferromag
double precision, dimension(100,100) :: matrix_ferro_to_invert, matrix_ferro_inverse

contains

!---------------------------------------------------------------------
    double precision function curinterp(r_in, z_in, jrho, rho, theta, rax, zax, nrho, ntheta)

    use pi_const, only: GP2
    use numerical_tools, only: linterp
    use feqis_tools, only: pol_angle

    integer,  intent(in) :: nrho, ntheta
    double precision, intent(in) :: r_in, z_in, rax, zax
    double precision, intent(in), dimension(ntheta) :: theta
    double precision, intent(in), dimension(nrho, ntheta) :: jrho, rho

    integer :: i, k1, k2, k3, k4, j1, j2
    double precision :: anglr, rho0, r1, r2, r3, r4, z1, z2, z3, z4, d1, d2, d3, d4
    double precision, dimension(4) :: jj1

    anglr = pol_angle(rax, zax, r_in, z_in)
    rho0 = sqrt((r_in - rax)**2 + (z_in - zax)**2)

    if (anglr < theta(1)) anglr = anglr + GP2

    j1 = 1
    do i=1, ntheta
        if (anglr >= theta(i)) j1 = i
    enddo
    j2 = j1 + 1

    z1 = rho(nrho, j1)
    z2 = rho(nrho, j2)
    if (rho0 > z1 .or. rho0 > z2) then
        curinterp = 0.
        return
    endif

    k1 = 1
    do i=1, nrho-1
        if (rho0 >= rho(i, j1)) k1 = i
    enddo
    k2 = k1 + 1

    k3 = 1
    do i=1, nrho-1
        if (rho0 >= rho(i, j2)) k3 = i
    enddo
    k4 = k3 + 1

    r1 = rax + rho(k1, j1)*cos(theta(j1))
    r2 = rax + rho(k2, j1)*cos(theta(j1))
    r3 = rax + rho(k4, j2)*cos(theta(j2))
    r4 = rax + rho(k3, j2)*cos(theta(j2))
    z1 = zax + rho(k1, j1)*sin(theta(j1))
    z2 = zax + rho(k2, j1)*sin(theta(j1))
    z3 = zax + rho(k4, j2)*sin(theta(j2))
    z4 = zax + rho(k3, j2)*sin(theta(j2))
    jj1(1) = jrho(k1, j1)
    jj1(2) = jrho(k2, j1)
    jj1(3) = jrho(k4, j2)
    jj1(4) = jrho(k3, j2)
    d1 = sqrt((r1 - r_in)**2 + (z1 - z_in)**2)
    d2 = sqrt((r2 - r_in)**2 + (z2 - z_in)**2)
    d3 = sqrt((r3 - r_in)**2 + (z3 - z_in)**2)
    d4 = sqrt((r4 - r_in)**2 + (z4 - z_in)**2)

    curinterp = (jj1(1)*d2*d3*d4 + jj1(2)*d1*d3*d4 + jj1(3)*d1*d2*d4 + jj1(4)*d1*d2*d3) / &
        (d1*d2*d3 + d1*d3*d4 + d2*d3*d4 + d1*d2*d4)

    end function curinterp

!---------------------------------------------------------------------
    subroutine psi_mutual_effect_conductors(psi_to_conductors)

    use green_function, only: greeni
    use fbe_core, only: dr, dz  

    double precision, intent(out) :: psi_to_conductors(nconduc)
    integer :: i

    do i=1, nconduc
        psi_to_conductors(i) = sum(jrz(1: nr2, 1: nz2) * dr * dz * greeni(1: nr2, 1: nz2, i))
    enddo

    end subroutine psi_mutual_effect_conductors

!---------------------------------------------------------------------
    subroutine psi_mutual_effect_conductors_simple(psi_to_conductors)

    use green_function, only: greeni

    double precision, intent(out) :: psi_to_conductors(nconduc)
    integer :: i, j, k
    
    i = nint((raxp - Rrect(1))/dr) + 1
    j = nint((zaxp - Zrect(1))/dz) + 1
    do k=1, nconduc
        psi_to_conductors(k) = iplasma * greeni(i, j, k)
    enddo

    end subroutine psi_mutual_effect_conductors_simple

!---------------------------------------------------------------------
    subroutine get_zccurb(rc_cur, zc_cur, z2c_cur, rgeoc, zgeoc, ahorc)

    use metric_coefficients_pbe, only: R_curr_0D, Z_curr_0D, dator

    double precision, intent(out) :: rc_cur, zc_cur, z2c_cur, rgeoc, zgeoc, ahorc

    integer :: i, j
    double precision :: perimz, ahorc2, avgelem

    perimz = 0.
    rgeoc  = 0.
    zgeoc  = 0.
    ahorc  = 0.
    ahorc2 = 0.
    zc_cur  = Z_curr_0D
    z2c_cur = zc_cur
    rc_cur  = R_curr_0D

    do i=1, nrho-1
        do j=1, ntheta
            avgelem = dator(i, j) !/bpcell(i, j-1)
            perimz = perimz + avgelem
            rgeoc  = rgeoc + rpol(i, j)*avgelem
            zgeoc  = zgeoc + zpol(i, j)*avgelem
        enddo
    enddo
    rgeoc = rgeoc/perimz
    zgeoc = zgeoc/perimz

    do i=1, nrho-1
        do j=1, ntheta
            avgelem = dator(i, j) !/bpcell(i, j-1)
            ahorc2  = ahorc2 + (rpol(i, j) - rgeoc)**2 * avgelem
        enddo
    enddo
    ahorc2 = ahorc2/perimz
    ahorc = 2.*sqrt(ahorc2)

    end subroutine get_zccurb

!---------------------------------------------------------------------
subroutine estimate_boundary_to_pbe(rbnd, zbnd, nthetaz)

    use pi_const, only: GP2
    use fbe_core, only: nr, nr2, nz, iaxis, jaxis, &
        raus, rinner, zbot, ztop, &
        Rrect, Zrect, dr, dz, rax, zax,  &
        psiaxis, psibnd, psirz
    use pbe_core, only: theta
    use feqis_tools, only: pol_angle, interp2d_psi

    integer, intent(IN) :: nthetaz
    double precision, intent(OUT), dimension(nthetaz) :: rbnd, zbnd

    integer :: i, j, k, j4, ntheta
    double precision :: x1, x2, t1, t2, t3, z1, z2, z3, x11, dx, dtheta
    double precision, dimension(500) :: theta_fbe

    ntheta = nthetaz
    rbnd = 0.
    zbnd = 0.

    do i=1, ntheta + 1
        theta(i) = GP2*(i - 1.)/(ntheta + 0.)
    enddo
    dtheta = theta(2) - theta(1)

    dx = sqrt(dr**2 + dz**2)

! Find boundary
    j = jaxis
    do i=iaxis, nr2
        if (psirz(i, j) >= psibnd) k = i
        if (psirz(i, j) <= psibnd) EXIT
    enddo

    if (psirz(k, j) == psibnd) then
        rbnd(1) = Rrect(k)
    else
        rbnd(1) = Rrect(k) - (psirz(k, j) - psibnd)/(psirz(k, j) - psirz(k-1, j))*dr
    endif
    zbnd(1) = Zrect(j)
    theta_fbe(1) = pol_angle(rax, zax, rbnd(1), zbnd(1))

    theta_loop: do i=2, ntheta
        dx = sqrt(dr**2 + dz**2)
        x1 = sqrt((rbnd(i-1) - rax)**2 + (zbnd(i-1) - zax)**2)
        theta_fbe(i) = theta_fbe(i-1) + dtheta

        t1 = rax + x1*cos(theta_fbe(i))
        t2 = zax + x1*sin(theta_fbe(i))

        if (t1 >= raus) then
            x1 = (raus - rax)/cos(theta_fbe(i))
        endif
        if (t1 <= rinner) then
            x1 = (rinner - rax)/cos(theta_fbe(i))
        endif
        if (t2 >= ztop) then
            x1 = (ztop - zax)/sin(theta_fbe(i))
        endif
        if (t2 <= zbot) then
            x1 = (zbot - zax)/sin(theta_fbe(i))
        endif
        t1 = rax + x1*cos(theta_fbe(i))
        t2 = zax + x1*sin(theta_fbe(i))
        t3 = interp2d_psi(t1, t2, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))

        if (t3 == psibnd) then
            rbnd(i) = t1
            zbnd(i) = t2
        elseif (t3 < psibnd) then
            do
                z1 = rax + (x1 - dx)*cos(theta_fbe(i))
                z2 = zax + (x1 - dx)*sin(theta_fbe(i))
                z3 = interp2d_psi(z1, z2, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                if (z3 < psibnd) then
                    dx = 1.1*dx
                else
                    EXIT
                endif
            enddo
            x11 = x1 - (t3 - psibnd)/(t3 - z3)*dx
            rbnd(i) = rax + x11*cos(theta_fbe(i))
            zbnd(i) = zax + x11*sin(theta_fbe(i))
        else
            do
                j4 = 0
                x2 = x1 + dx
                z1 = rax + x2*cos(theta_fbe(i))
                z2 = zax + x2*sin(theta_fbe(i))
                if (z1 >= raus) then
                    x2 = (raus - rax)/cos(theta_fbe(i))
                    j4 = 1
                endif
                if (z1 <= rinner) then
                    x2 = (rinner - rax)/cos(theta_fbe(i))
                    j4 = 1
                endif
                if (z2 >= ztop) then
                    x2 = (ztop - zax)/sin(theta_fbe(i))
                    j4 = 1
                endif
                if (z2 <= zbot) then
                    x2 = (zbot - zax)/sin(theta_fbe(i))
                    j4 = 1
                endif
                z1 = rax + x2*cos(theta_fbe(i))
                z2 = zax + x2*sin(theta_fbe(i))
                z3 = interp2d_psi(z1, z2, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                if (z3 > psibnd) then
                    if (j4 == 1) then
                        rbnd(i) = z1
                        zbnd(i) = z2
                        CYCLE theta_loop
                    endif
                    dx = dx*1.1
                else
                    EXIT
                endif
            enddo
            x11 = x1 + (t3 - psibnd)/(t3 - z3)*dx
            rbnd(i) = rax + x11*cos(theta_fbe(i))
            zbnd(i) = zax + x11*sin(theta_fbe(i))
       endif

    enddo theta_loop

    end subroutine estimate_boundary_to_pbe

!---------------------------------------------------------------------
    subroutine coil_forces(ncoilz, force_R, force_Z, plasma_state)
! this one is only between coils and coils

    use fbe_core, only: jrz, nr2, nz2, dr, dz, curconduc
    use green_matrix, only: dgreenirpl, dgreenizpl, dgreenirj, dgreenizj

    integer, intent(in) :: ncoilz, plasma_state
    double precision, intent(out), dimension(ncoilz) :: force_R, force_Z

    integer :: i, j, nblock_a
    double precision :: x1

    force_R = 0.
    force_Z = 0.
    nblock_a = nblocks - npassive

    if (plasma_state == 1) then !not sure about the plasma response...
        do i=1, nblock_a
            x1 =  sum(jrz(1:nr2, 1:nz2) * dr * dz * dgreeniRpl(1:nr2, 1:nz2, i))
            force_R(i) = force_R(i) + curconduc(mequivalence(i)) * x1
            x1 =  sum(jrz(1:nr2, 1:nz2) * dr * dz * dgreeniZpl(1:nr2, 1:nz2, i))
            force_Z(i) = force_Z(i) + curconduc(mequivalence(i)) * x1
        enddo
    endif

!block-to-block
    do i=1, nblock_a
        do j=1, nblock_a
            if (i /= j) then
                force_R(i) = force_R(i) + curconduc(mequivalence(j)) * curconduc(mequivalence(i)) * dgreeniRj(i, j)
                force_Z(i) = force_Z(i) + curconduc(mequivalence(j)) * curconduc(mequivalence(i)) * dgreeniZj(i, j)
            endif
        enddo
    enddo

! force_R and force_Z are F_R and F_Z components in [N] for each block , does not include forces from the passive elements or on the passive elements.

    end subroutine coil_forces

!---------------------------------------------------------------------
    subroutine get_zccurfbe(rc_cur, zc_cur)

    double precision, intent(out) :: rc_cur, zc_cur

    integer :: i, j
    double precision :: curtotal

    curtotal = sum(jrz)
    zc_cur = 0.
    rc_cur = 0.
    do j=1, nz2
        do i=1, nr2
            zc_cur = zc_cur + jrz(i,j)*Zrect(j)
            rc_cur = rc_cur + jrz(i,j)*Rrect(i)
        enddo
    enddo
    zc_cur = zc_cur/curtotal
    rc_cur = rc_cur/curtotal

    end subroutine get_zccurfbe

!---------------------------------------------------------------------
    double precision function psib_ext
! Returns external flux on plasma boundary

    use feqis_tools, only: interp2d_psi

    integer :: i
    double precision :: psiext_out, dlt, dllt, psi_ext_1, psiext1, psiext2

!cycle over boundary
    psiext_out = 0.
    dllt = 0.
    psi_ext_1 = interp2d_psi(rbnd(1), zbnd(1), Rrect(1:nr), Zrect(1:nz), psiextrz(1:nr, 1:nz))
    psiext1 = psi_ext_1
    do i=1, nbnd-1
        psiext2 = interp2d_psi(rbnd(i+1), zbnd(i+1), Rrect(1:nr), Zrect(1:nz), psiextrz(1:nr, 1:nz))
        dlt = sqrt((rbnd(i+1) - rbnd(i))**2 + (zbnd(i+1) - zbnd(i))**2)
        psiext_out = psiext_out + 0.5*(psiext1 + psiext2)*dlt
        dllt = dllt + dlt
        psiext1 = psiext2
    enddo
    psiext1 = psi_ext_1
    dlt = sqrt((rbnd(1) - rbnd(nbnd))**2 + (zbnd(1) - zbnd(nbnd))**2)
    psiext_out = psiext_out + 0.5*(psiext1 + psiext2)*dlt
    dllt = dllt + dlt

    psib_ext = psiext_out/dllt

    end function psib_ext

!---------------------------------------------------------------------
    double precision function find_l_gap(psibnd, l_ref_in, gapmin, gapmax, geom)

    use errors_params, only: err_gaptolez
    use feqis_tools, only: interp2d_psi

    double precision, intent(in) :: psibnd, l_ref_in, gapmin, gapmax, geom(3)

    double precision :: dur1, dur2, u001, u002, x1, x2, y1, y2, l_ref

    l_ref = l_ref_in
    dur1 = 0.
    dur2 = l_ref

    do
        x1 = geom(1) + dur1*cos(geom(3))
        y1 = geom(2) + dur1*sin(geom(3))
        x2 = geom(1) + dur2*cos(geom(3))
        y2 = geom(2) + dur2*sin(geom(3))
        u001 = interp2d_psi(x1, y1, Rrect(1:nr), Zrect(1:nz), psirz(1: nr, 1:nz))
        u002 = interp2d_psi(x2, y2, Rrect(1:nr), Zrect(1:nz), psirz(1: nr, 1:nz))
        if (abs(l_ref) < err_gaptolez) then
            find_l_gap = 0.5*(dur1 + dur2)
            EXIT
        endif
        if (u001 == psibnd) then
            find_l_gap = dur1
            EXIT
        endif
        if (u002 == psibnd) then
            find_l_gap = dur2
            EXIT
        endif
        if ((u002 > psibnd .and. u001 < psibnd) .or. (u002 < psibnd .and. u001 > psibnd)) then
            l_ref = -0.5*l_ref
            dur1 = dur2
            dur2 = dur1 + l_ref
            CYCLE
        endif
! Case with no intersection: larger 1. m
        find_l_gap = 0.5*(dur1 + dur2)
        if (find_l_gap >= gapmax .or. find_l_gap < gapmin) then
            find_l_gap = -5000.
            EXIT
        endif
        dur1 = dur1 + l_ref
        dur2 = dur2 + l_ref
    enddo

    end function find_l_gap

!---------------------------------------------------------------------
    function find_demo_gaps(ngaps, demo_gaps) result(geom1d)

    integer, intent(in) :: ngaps
    double precision, intent(in) :: demo_gaps(ngaps, 4)
    double precision :: geom1d(ngaps)

    integer :: i, onlypos
    double precision :: d_step, gapmin, gapmax, l_gap, l_gap_pos, l_gap_neg

    d_step = 0.01 ! advance in 1 cm steps
    gapmin = -1.5
    gapmax = 2.5

    do i=1, ngaps
        onlypos = nint(demo_gaps(i, 4))
        l_gap_pos = find_l_gap(psibnd,  d_step, gapmin, gapmax, demo_gaps(i, 1:3))
        l_gap_neg = find_l_gap(psibnd, -d_step, gapmin, gapmax, demo_gaps(i, 1:))
! choose minimum of absolute values
        if (onlypos == 0) then
            if (abs(l_gap_neg) < abs(l_gap_pos)) then
                l_gap = l_gap_neg
            else
                l_gap = l_gap_pos
            endif
            if (isnan(l_gap_pos)) l_gap = l_gap_neg
            if (isnan(l_gap_neg)) l_gap = l_gap_pos
            geom1d(i) = min(gapmax, max(gapmin, l_gap))
        else
            geom1d(i) = min(gapmax, max(gapmin, l_gap_pos))
        endif
    enddo

    end function find_demo_gaps

!---------------------------------------------------------------------
    subroutine psi_external_calc_position(r_target, z_target, psi_target)

    use green_function, only: greeni
    
    double precision, intent(IN) :: r_target, z_target
    double precision, intent(OUT) :: psi_target
    integer :: i, j

    i = nint((r_target - Rrect(1))/dr) + 1   
    j = nint((z_target - Zrect(1))/dz) + 1   
    psi_target = sum(curconduc(1: nconduc)*greeni(i, j, 1: nconduc))

    end subroutine psi_external_calc_position

!--------------------------------------------------------------------
    subroutine restab_F_function_full_fonfit
! Refits all currents

    use errors_params, only: err_find_psistab
    use transport2fbe, only: sigma_coils, sigma_b, sigma_axis
    use green_function, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix

    integer :: i, j, k, j_iter, iax, jax
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, raxref, zaxref
    double precision, dimension(9) :: bub
    double precision, dimension(500) :: rbref, zbref
    double precision, dimension(nconduc) :: G_00c, G_00r, G_00z, &
        Fderiv, curref, curnow, curdiff
    double precision, dimension(ntheta) :: psicorr
    double precision, dimension(nconduc, ntheta) :: G_00
    double precision, dimension(nconduc, nconduc) :: matrix, invmatrix

    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhotheta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, Rrect(1), dr)
    jaxis = closest_index(zax, Zrect(1), dz)
    iax = iaxis
    jax = jaxis

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    sigma_coils(nactive + 1:nconduc) = sigma_coils(nactive + 1)

    curref(1:nconduc) = curconduc(1:nconduc)
    curnow(1:nconduc) = curconduc(1:nconduc)
    curdiff = 0.

    raxref = raxp
    zaxref = zaxp
    rbref(1:ntheta) = rbndp(1:ntheta)
    zbref(1:ntheta) = zbndp(1:ntheta)

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.
    do j=1, nconduc
        do k=1, ntheta
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:ntheta))/(0. + ntheta)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
    enddo

    do j=1, nconduc
        do i=1, nconduc
            if (i == j) matrix(i, j) = matrix(i, j) + 2.*sigma_coils(i)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:ntheta) - G_00c(i))*(G_00(j, 1:ntheta) - G_00c(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j))
       enddo
    enddo

! Calculate inverse
    invmatrix = inv_matrix(matrix, nconduc)

    do j_iter=1, 300000 !iterations to find currents
        if (j_iter > 150) stop
        psiplasrz = get_psiplasrz()

! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nconduc
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        do j=1, ntheta
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(ntheta + 0.) !average psi on the boundary

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz

        Ffunc = sigma_B*sum((psicorr - x1)**2) + sum(sigma_coils(1:nconduc)*curdiff**2) + sigma_axis*(x2**2 + x3**2)

! Calculate F derivative
        do i=1, nconduc
            Fderiv(i) = 2.*(sigma_coils(i)*curdiff(i) + &
                sigma_B*sum((psicorr - x1)*(G_00(i, 1:ntheta) - G_00c(i))) + &
                sigma_axis*(x2*G_00r(i) + x3*G_00z(i)) )
        enddo

! Calculate new currents

        do i=1, nconduc
            curnow(i) = curnow(i) - sum(invmatrix(i, 1:nconduc)*Fderiv)
        enddo
        curdiff = curnow - curref

! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nconduc
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        call find_new_axis
        call find_psi_boundary
        call new_jrz  ! calculate new right hand side

        temp_err = abs(Ffunc - Ffunc_old)
        Ffunc_old = Ffunc
        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nconduc
        curconduc(i) = curnow(i)
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis
    call find_psi_boundary
    call new_jrz  ! calculate new right hand side

    end subroutine restab_F_function_full_fonfit

!--------------------------------------------------------------------
    subroutine restab_F_function_full_fonfit_xpoints ! valid only if n_xpoint_fit > 0
! Refits all currents and x points

    use errors_params, only: err_find_psistab
    use transport2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_xpoint, &
        r_xpoint_fit, z_xpoint_fit, n_xpoint_fit
    use green_function, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix

    integer :: i, j, k, j_iter, iax, jax
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, raxref, zaxref
    double precision, dimension(9) :: bub
    double precision, dimension(500) :: rbref, zbref
    double precision, dimension(nconduc) :: G_00c, G_00r, G_00z, &
        Fderiv, curref, curnow, curdiff

    double precision, dimension(nconduc, n_xpoint_fit) :: G_00xr, G_00xz
    double precision, dimension(n_xpoint_fit) :: dummyx

    double precision, dimension(ntheta) :: psicorr
    double precision, dimension(nconduc, ntheta) :: G_00
    double precision, dimension(nconduc, nconduc) :: matrix, invmatrix

    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhotheta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, Rrect(1), dr)
    jaxis = closest_index(zax, Zrect(1), dz)
    iax = iaxis
    jax = jaxis

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    sigma_coils(nactive + 1:nconduc) = sigma_coils(nactive + 1)

    curref(1:nconduc) = curconduc(1:nconduc)
    curnow(1:nconduc) = curconduc(1:nconduc)
    curdiff = 0.

    raxref = raxp
    zaxref = zaxp
    rbref(1:ntheta) = rbndp(1:ntheta)
    zbref(1:ntheta) = zbndp(1:ntheta)

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.
    do j=1, nconduc
        do k=1, ntheta
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:ntheta))/(0. + ntheta)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
        do k=1, n_xpoint_fit
            bub(1) = interp2d_psi(r_xpoint_fit(k) - dr/2., z_xpoint_fit(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            bub(2) = interp2d_psi(r_xpoint_fit(k) + dr/2., z_xpoint_fit(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(k) - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(k) + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            G_00xr(j,k) = (bub(2) - bub(1))/dr
            G_00xz(j,k) = (bub(4) - bub(3))/dz
        enddo
    enddo

    do j=1, nconduc
        do i=1, nconduc
            if (i == j) matrix(i, j) = matrix(i, j) + 2.*sigma_coils(i)
            dummyx(1:n_xpoint_fit) = G_00xr(i,1:n_xpoint_fit)*G_00xr(j,1:n_xpoint_fit) + G_00xz(i,1:n_xpoint_fit)*G_00xz(j,1:n_xpoint_fit)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:ntheta) - G_00c(i))*(G_00(j, 1:ntheta) - G_00c(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j)) + &
                2.*sigma_xpoint*sum(dummyx)
        enddo
    enddo

! Calculate inverse
    invmatrix = inv_matrix(matrix, nconduc)

    do j_iter=1, 300000 !iterations to find currents
        if (j_iter > 150) stop
        psiplasrz = get_psiplasrz()

! Construct correction
        call compound_psi
    do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nconduc
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        do j=1, ntheta
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(ntheta + 0.) !average psi on the boundary

        Ffunc = 0.
        Ffunc = sigma_B*sum((psicorr - x1)**2) + sum(sigma_coils(1:nconduc)*curdiff**2)

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz
        Ffunc = Ffunc +  sigma_axis*(x2**2 + x3**2)
        do k=1, n_xpoint_fit
            bub(1) = interp2d_psi(r_xpoint_fit(k) - 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(2) = interp2d_psi(r_xpoint_fit(k) + 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            x2 = (bub(2) - bub(1))/dr ! dPsi/dr
            x3 = (bub(4) - bub(3))/dz ! dPsi/dz
            Ffunc = Ffunc +  sigma_xpoint*(x2**2 + x3**2)
        enddo

! Calculate F derivative
        do i=1, nconduc
            Fderiv(i) = 2.*(sigma_coils(i)*curdiff(i) + sigma_B*sum((psicorr - x1)*(G_00(i, 1:ntheta) - G_00c(i))))
            bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            x2 = (bub(2) - bub(1))/dr ! dPsi/dr
            x3 = (bub(4) - bub(3))/dz ! dPsi/dz
            Fderiv(i) = Fderiv(i) + 2.*sigma_axis*(x2*G_00r(i) + x3*G_00z(i))
            do k=1, n_xpoint_fit
                bub(1) = interp2d_psi(r_xpoint_fit(k) - 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                bub(2) = interp2d_psi(r_xpoint_fit(k) + 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                x2 = (bub(2) - bub(1))/dr ! dPsi/dr
                x3 = (bub(4) - bub(3))/dz ! dPsi/dz
                Fderiv(i) = Fderiv(i) +  2.*sigma_xpoint*(x2*G_00xr(i,k) + x3*G_00xz(i,k))
            enddo
        enddo

! Calculate new currents

        do i=1, nconduc
            curnow(i) = curnow(i) - sum(invmatrix(i, 1:nconduc)*Fderiv)
        enddo
        curdiff = curnow - curref

! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nconduc
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        call find_new_axis
        call find_psi_boundary
        call new_jrz  ! calculate new right hand side

        temp_err = abs(Ffunc - Ffunc_old)
        Ffunc_old = Ffunc
        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nconduc
        curconduc(i) = curnow(i)
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis
    call find_psi_boundary
    call new_jrz  ! calculate new right hand side

    end subroutine restab_F_function_full_fonfit_xpoints

!--------------------------------------------------------------------
    subroutine restab_2_timepoints_evolution

! Finds active currents from scratch including evolution from time point t1 to time point t2

    use errors_params, only: err_find_psistab
    use transport2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy
    use green_function, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix
    use pi_const, only: GP2, GP2_sq, mu0
    use numerical_tools, only: linterp

    integer, parameter :: n_evol = 2
    integer :: i, j, k, j_iter, jt
    integer, dimension(n_evol) :: iax_ev, jax_ev, nxp_ev
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, lambda, deltapsiext !lagrange multplier lambda
    double precision, dimension(9) :: bub
    double precision, dimension(n_evol) :: rax_ev,zax_ev,ip_ev, V_loop, L_ext, delta_t, psi_ext_ev
    double precision, dimension(ntheta) :: psicorr
    double precision, dimension(:), allocatable :: raxref_ev, zaxref_ev, &
        ffunc_ev, result_vector, Fderiv
    double precision, dimension(nrho, n_evol) :: psia_ev,pprim_ev,ffprim_ev
    double precision, dimension(nactive, n_evol) :: i_totev
    double precision, dimension(1000, n_evol) :: rxp_ev, zxp_ev
    double precision, dimension(:, :), allocatable :: rbref_ev, zbref_ev, G_00c, G_00r, G_00z, &
        curref, curnow, curdiff, matrix, invmatrix
    double precision, dimension(:, :, :), allocatable :: j_ev, psi_ev, G_00
    character(len=80) :: file_time

    allocate(rbref_ev(ntheta, n_evol))
    allocate(zbref_ev(ntheta, n_evol))
    allocate(raxref_ev(n_evol))
    allocate(zaxref_ev(n_evol))
    allocate(ffunc_ev(n_evol))
    allocate(j_ev(nr2, nz2, n_evol))
    allocate(psi_ev(nr2, nz2, n_evol))
    allocate(G_00(nactive, ntheta, n_evol))
    allocate(G_00c(nactive, n_evol))
    allocate(G_00r(nactive, n_evol))
    allocate(G_00z(nactive, n_evol))
    allocate(curref(nactive, n_evol))
    allocate(curnow(nactive, n_evol))
    allocate(curdiff(nactive, n_evol))
    allocate(result_vector(nactive*n_evol+1))
    allocate(Fderiv(nactive*n_evol+1))
    allocate(matrix(nactive*n_evol+1, nactive*n_evol+1))
    allocate(invmatrix(nactive*n_evol+1, nactive*n_evol+1))

    do jt=1, 2
        if (jt == 1) file_time = 'dat/feqis_time1_file.dat' !contains the prescribed boundary data for time t1: jrhotheta, rho, theta, rb, zb, rho augmented to ntheta+1, theta augmented to ntheta+1, raxp, zaxp, iplasma, psia,pprim,ffprim
        if (jt == 2) file_time = 'dat/feqis_time2_file.dat' !contains the prescribed boundary data for time t2: jrhotheta, rho, theta, rho augmented to ntheta+1, theta augmented to ntheta+1, raxp, zaxp, iplasma, psia2,pprim2,ffprim2, Lext, deltaT, Vloop
        open(32, file=TRIM(file_time))
            read(32, *) jrhotheta(1:nrho, 1:ntheta)
            read(32, *) rho(1:nrho, 1:ntheta)
            read(32, *) theta(1:ntheta)
            read(32, *) rbndp(1:ntheta)
            read(32, *) zbndp(1:ntheta)
            read(32, *) raxp, zaxp, ip_ev(jt)
            read(32, *) psia_ev(1:nrho, jt), pprim_ev(1:nrho, jt), ffprim_ev(1:nrho, jt) !pprime is Pascal / grad(FP), ffprime is F dF/dFP
            read(32, *) L_ext(jt), delta_t(jt), V_loop(jt)
        close(32)
        rho(1:nrho, ntheta+1) = rho(1:nrho, 1)
        theta(ntheta+1) = theta(1) + GP2
        call interp_j_fromrhotorz
        curr = SUM(jrz)*dr*dz
        j_ev(:, :, jt) = jrz/curr*ip_ev(jt)
        psia_ev(:, jt) = (psia_ev(:, jt) - psia_ev(1, jt))/(psia_ev(nrho, jt) - psia_ev(1, jt))
        rax_ev(jt) = raxp
        iax_ev(jt) = closest_index(rax_ev(jt), Rrect(1), dr)
        zax_ev(jt) = zaxp
        jax_ev(jt) = closest_index(zax_ev(jt), Zrect(1), dz)
        rbref_ev(1: ntheta, jt) = rbndp(1: ntheta)
        zbref_ev(1: ntheta, jt) = zbndp(1: ntheta)
    enddo

    deltapsiext = -(0.5*sum(V_loop)*(delta_t(2) - delta_t(1)) + 0.5*sum(L_ext)*(ip_ev(2) - ip_ev(1)))

    psicorr   = 0.
    Ffunc_old = 1.e6
    Fderiv = 0.

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    curref = 0.
    curnow = 0.
    curdiff = 0.
    result_vector = 0.

    raxref_ev = rax_ev
    zaxref_ev = zax_ev

    nxp_ev = 0
    rxp_ev = 0.
    zxp_ev = 0.

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.
    psiextrz = 0.  !assume total vacuum, no eddy currents
    lambda = 0.

! build up bordered Hessian matrix (1...I1, 1....I2, lambda ; same)
    do jt=1,n_evol
        do j=1, nactive
            do k=1, ntheta
                G_00(j, k, jt) = interp2d_psi(rbref_ev(k,jt), zbref_ev(k,jt), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            enddo
            G_00c(j, jt) = sum(G_00(j, 1:ntheta, jt))/(0. + ntheta)
            bub(1) = interp2d_psi(raxref_ev(jt) - dr/2., zaxref_ev(jt), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            bub(2) = interp2d_psi(raxref_ev(jt) + dr/2., zaxref_ev(jt), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            bub(3) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            bub(4) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            G_00r(j,jt) = (bub(2) - bub(1))/dr
            G_00z(j,jt) = (bub(4) - bub(3))/dz
        enddo

        do j=1, nactive
            do i=1, nactive
                if (i == j) matrix((jt-1)*nactive+i, (jt-1)*nactive+j) = matrix((jt-1)*nactive+i, (jt-1)*nactive+j) + sigma_coils(i)*sigma_energy*indconduc(i, i)
                matrix((jt-1)*nactive+i, (jt-1)*nactive+j) = matrix((jt-1)*nactive+i, (jt-1)*nactive+j) +  &
                    2.*sigma_B*sum((G_00(i, 1:ntheta, jt) - G_00c(i, jt))*(G_00(j, 1:ntheta, jt) - G_00c(j, jt))) +  &
                    2.*sigma_axis*(G_00r(i, jt)*G_00r(j, jt) + G_00z(i, jt)*G_00z(j, jt))
            enddo
        enddo
    enddo

! Borders with lambda
    do i=1, nactive
        matrix(2*nactive+1, i) = G_00c(i, 1)
        matrix(i, 2*nactive+1) = G_00c(i, 1)
    enddo
    do i=nactive+1, 2*nactive
        matrix(2*nactive+1, i) = -G_00c(i-nactive, 2)
        matrix(i, 2*nactive+1) = -G_00c(i-nactive, 2)
    enddo
    matrix(2*nactive+1, 2*nactive+1) = 0.

! Calculate inverse
    invmatrix = inv_matrix(matrix, 2*nactive+1)

    do j_iter=1, 300000 !iterations to find currents
        if (j_iter > 150) stop
        do jt=1, n_evol
            jrz = j_ev(:, :, jt)   ! assign previous current density
            psiplasrz = get_psiplasrz()
            psi_ev(:, :, jt) = psiplasrz
! Construct correction
!  call compound_psi
            do j=1, nz2
                do i=1, nr2
                    f_correction = 0.
                    do k=1, nactive
                        f_correction = f_correction + curdiff(k, jt)*greeni(i, j, k)
                    enddo
                   psirz(i, j) = psi_ev(i, j, jt) + f_correction !total flux
               enddo
           enddo
           do j=1, ntheta
               psicorr(j) = interp2d_psi(rbref_ev(j,jt), zbref_ev(j,jt), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
           enddo
           x1 = sum(psicorr)/(ntheta + 0.) !average psi on the boundary
! Derivative at ref axis
           bub(1) = interp2d_psi(raxref_ev(jt) - 0.5*dr, zaxref_ev(jt), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
           bub(2) = interp2d_psi(raxref_ev(jt) + 0.5*dr, zaxref_ev(jt), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
           bub(3) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
           bub(4) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
           x2 = (bub(2) - bub(1))/dr ! dPsi/dr
           x3 = (bub(4) - bub(3))/dz ! dPsi/dz

           psi_ext_ev(jt) = sum(G_00c(:,jt)*curdiff(:,jt))
           Ffunc_ev(jt) = sigma_B*sum((psicorr - x1)**2) + sigma_axis*(x2**2 + x3**2)
           do i=1, nactive
               Ffunc_ev(jt) = Ffunc_ev(jt) + 0.5*indconduc(i, i)*sigma_coils(i)*sigma_energy*curdiff(i, jt)**2
           enddo

! Calculate F derivative
           do i=1, nactive
               Fderiv((jt-1)*nactive+i) = 2.*(0.5*indconduc(i, i)*sigma_coils(i)*sigma_energy*curdiff(i, jt) + &
                   sigma_B*sum((psicorr - x1)*(G_00(i, 1:ntheta, jt) - G_00c(i, jt))) + &
                   sigma_axis*(x2*G_00r(i, jt) + x3*G_00z(i, jt)) )
           enddo
       enddo !end time loop

! add lambda contributions
       Ffunc = sum(Ffunc_ev) + lambda*(deltapsiext - (psi_ext_ev(2) - psi_ext_ev(1)))
       do i=1, nactive
           Fderiv(i) = Fderiv(i) + lambda*G_00c(i, 1)
       enddo
       do i=1, nactive
           Fderiv(nactive+i) = Fderiv(nactive+i) - lambda*G_00c(i, 2)
       enddo
       Fderiv(2*nactive+1)=deltapsiext-(psi_ext_ev(2)-psi_ext_ev(1)) ! dF/dlambda

! Calculate new currents
       do i=1, n_evol*nactive+1
           result_vector(i) = result_vector(i) - sum(invmatrix(i, 1:n_evol*nactive+1)*Fderiv)
       enddo
       curdiff(:, 1) = result_vector(1:nactive)
       curdiff(:, 2) = result_vector(nactive+1:n_evol*nactive)
       lambda = result_vector(n_evol*nactive+1)

! Update plasma current density field
       do jt=1,n_evol
! Construct correction
           do j=1, nz2
               do i=1, nr2
                   f_correction = 0.
                   do k=1, nactive
                       f_correction = f_correction + curdiff(k, jt)*greeni(i, j, k)
                   enddo
                   psirz(i, j) = psi_ev(i, j, jt) + f_correction ! total flux
               enddo
           enddo

! axis block
           iaxis = iax_ev(jt)
           jaxis = jax_ev(jt)
           call find_new_axis
           iax_ev(jt) = iaxis
           jax_ev(jt) = jaxis

! boundary block
           n_of_xpoints = nxp_ev(jt)
           if (n_of_xpoints > 0) then
               r_xpoint(1: n_of_xpoints) = rxp_ev(1: n_of_xpoints, jt)
               z_xpoint(1: n_of_xpoints) = zxp_ev(1: n_of_xpoints, jt)
           endif
           call find_psi_boundary
           nxp_ev(jt) = n_of_xpoints
           if (n_of_xpoints > 0) then
               rxp_ev(1: n_of_xpoints, jt) = r_xpoint(1: n_of_xpoints)
               zxp_ev(1: n_of_xpoints, jt) = z_xpoint(1: n_of_xpoints)
           endif

! current block
           call linterp(psia_ev(:, jt), ffprim_ev(:, jt), nrho, psia_2d, ffp_2d, nrho2d)
           call linterp(psia_ev(:, jt), pprim_ev (:, jt), nrho, psia_2d, ppp_2d, nrho2d)
           ffp_2d = -GP2/mu0*ffp_2d
           ppp_2d = -GP2*1.e-6*ppp_2d
           call new_jrz  ! calculate new right hand side
           j_ev(:, :, jt) = jrz
           i_totev(1: nactive, jt) = curdiff(1: nactive, jt)
       enddo

       temp_err = sum(abs(Fderiv)) !error
       Ffunc_old = Ffunc
       if (temp_err <= err_find_psistab) EXIT
    enddo

    open(32, file='dat/output_timefit.dat')
        write(32,*) 't1 currents: ', i_totev(1: nactive, 1)
        write(32,*) 't2 currents: ', i_totev(1: nactive, 2)
    close(32)
    pause

    end subroutine restab_2_timepoints_evolution

!--------------------------------------------------------------------
    subroutine restab_2_timepoints_evolution_limits

! Finds active currents from scratch including evolution from time point t1 to time point t2

    use errors_params, only: err_find_psistab
    use transport2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy, &
        current_limit
    use green_function, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix
    use pi_const, only: GP2, GP2_sq, mu0
    use numerical_tools, only: linterp

    integer, parameter :: n_evol = 2
    integer :: i, j, k, j_iter, jt, iactive, jactive
    integer, dimension(n_evol) :: iax_ev, jax_ev, nxp_ev
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, lambda, deltapsiext !lagrange multplier lambda
    double precision, dimension(9) :: bub
    double precision, dimension(ntheta) :: psicorr
    double precision, dimension(n_evol) :: rax_ev, zax_ev, ip_ev, &
         V_loop, L_ext, delta_t, psi_ext_ev
    double precision, dimension(:), allocatable :: raxref_ev, zaxref_ev, &
         ffunc_ev, result_vector, Fderiv
    double precision, dimension(nrho, n_evol) :: psia_ev,pprim_ev,ffprim_ev
    double precision, dimension(1000, n_evol) :: rxp_ev, zxp_ev
    double precision, dimension(nactive, n_evol) :: i_totev
    double precision, dimension(:, :), allocatable :: rbref_ev, zbref_ev, &
        G_00c, G_00r, G_00z, curref, curnow, curdiff, matrix, invmatrix
    double precision, dimension(:, :, :), allocatable :: j_ev, psi_ev, G_00
    character(len=80) :: file_time

    allocate(rbref_ev(ntheta, n_evol))
    allocate(zbref_ev(ntheta, n_evol))
    allocate(raxref_ev(n_evol))
    allocate(zaxref_ev(n_evol))
    allocate(ffunc_ev(n_evol))
    allocate(j_ev(nr2, nz2, n_evol))
    allocate(psi_ev(nr2, nz2, n_evol))
    allocate(G_00(nactive, ntheta, n_evol))
    allocate(G_00c(nactive, n_evol))
    allocate(G_00r(nactive, n_evol))
    allocate(G_00z(nactive, n_evol))
    allocate(curref(nactive, n_evol))
    allocate(curnow(nactive, n_evol))
    allocate(curdiff(nactive, n_evol))
    allocate(result_vector(nactive*n_evol+1))
    allocate(Fderiv(nactive*n_evol+1))
    allocate(matrix(nactive*n_evol+1, nactive*n_evol+1))
    allocate(invmatrix(nactive*n_evol+1, nactive*n_evol+1))

    do jt=1, 2
        if (jt==1) file_time='dat/feqis_time1_file.dat'   !contains the prescribed boundary data for time t1: jrhotheta, rho, theta, rb, zb, rho augmented to ntheta+1, theta augmented to ntheta+1, raxp, zaxp, iplasma, psia,pprim,ffprim
        if (jt==2) file_time='dat/feqis_time2_file.dat'   !contains the prescribed boundary data for time t2: jrhotheta, rho, theta, rho augmented to ntheta+1, theta augmented to ntheta+1, raxp, zaxp, iplasma, psia2,pprim2,ffprim2, Lext, deltaT, Vloop
        open(32, file=TRIM(file_time))
            read(32, *) jrhotheta(1: nrho, 1: ntheta)
            read(32, *) rho(1: nrho, 1: ntheta)
            read(32, *) theta(1: ntheta)
            read(32, *) rbndp(1: ntheta)
            read(32, *) zbndp(1: ntheta)
            read(32, *) raxp, zaxp, ip_ev(jt)
            read(32, *) psia_ev(1: nrho,jt), pprim_ev(1: nrho, jt), ffprim_ev(1: nrho, jt) !pprime is Pascal / grad(FP), ffprime is F dF/dFP
            read(32,*) L_ext(jt), delta_t(jt), V_loop(jt)
        close(32)
        rho(1: nrho, ntheta+1) = rho(1: nrho, 1)
        theta(ntheta+1) = theta(1) + GP2
        call interp_j_fromrhotorz
        curr = SUM(jrz) *dr*dz
        j_ev(:, :, jt) = jrz/curr*ip_ev(jt)
        psia_ev(:, jt) = (psia_ev(:, jt) - psia_ev(1, jt))/(psia_ev(nrho, jt) - psia_ev(1, jt))
        rax_ev(jt) = raxp
        iax_ev(jt) = closest_index(rax_ev(jt), Rrect(1), dr)
        zax_ev(jt) = zaxp
        jax_ev(jt) = closest_index(zax_ev(jt), Zrect(1), dz)
        rbref_ev(1: ntheta, jt) = rbndp(1: ntheta)
        zbref_ev(1: ntheta, jt) = zbndp(1: ntheta)
    enddo

    deltapsiext = -(0.5*sum(V_loop)*(delta_t(2) - delta_t(1)) + 0.5*sum(L_ext) * (ip_ev(2) - ip_ev(1)))

    psicorr   = 0.
    Ffunc_old = 1.e6
    Fderiv = 0.

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    curref = 0.
    curnow = 0.
    curdiff = 0.
    result_vector = 0.
    raxref_ev = rax_ev
    zaxref_ev = zax_ev
    nxp_ev = 0
    rxp_ev = 0.
    zxp_ev = 0.

! Calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.
    psiextrz = 0.  !assume total vacuum, no eddy currents
    lambda = 0.

! Build up bordered Hessian matrix (1...I1, 1....I2, lambda ; same)
    do jt=1, n_evol
        do j=1, nactive
            do k=1, ntheta
                G_00(j, k, jt) = interp2d_psi(rbref_ev(k, jt), zbref_ev(k, jt), Rrect(1: nr), Zrect(1: nz), greeni(1: nr, 1: nz, j))
            enddo
            G_00c(j, jt) = sum(G_00(j, 1:ntheta, jt))/(0. + ntheta)
            bub(1) = interp2d_psi(raxref_ev(jt) - dr/2., zaxref_ev(jt), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            bub(2) = interp2d_psi(raxref_ev(jt) + dr/2., zaxref_ev(jt), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            bub(3) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            bub(4) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            G_00r(j, jt) = (bub(2) - bub(1))/dr
            G_00z(j, jt) = (bub(4) - bub(3))/dz
        enddo

        do j=1, nactive
            jactive = (jt - 1)*nactive + j
            do i=1, nactive
               iactive = (jt - 1)*nactive + i
               if (i == j) then
                   matrix(iactive, jactive) = matrix(iactive, jactive) + sigma_coils(i)*sigma_energy*indconduc(i, i)
               endif
               matrix(iactive, jactive) = matrix(iactive, jactive) +  &
                   2.*sigma_B*sum((G_00(i, 1:ntheta, jt) - G_00c(i, jt))*(G_00(j, 1:ntheta, jt) - G_00c(j, jt))) +  &
                   2.*sigma_axis*(G_00r(i, jt)*G_00r(j, jt) + G_00z(i, jt)*G_00z(j, jt))
            enddo
        enddo
    enddo

! Borders with lambda
    do i=1,nactive
        matrix(2*nactive+1, i) = G_00c(i, 1)
        matrix(i, 2*nactive+1) = G_00c(i, 1)
    enddo
    do i=nactive+1, 2*nactive
        matrix(2*nactive+1, i) = -G_00c(i-nactive, 2)
        matrix(i, 2*nactive+1) = -G_00c(i-nactive, 2)
    enddo
    matrix(2*nactive+1, 2*nactive+1) = 0.

! Calculate inverse
    invmatrix = inv_matrix(matrix, 2*nactive+1)

    do j_iter=1, 300000 ! iterations to find currents
        if (j_iter > 150) stop
        do jt=1, n_evol
            jrz = j_ev(:, :, jt)    ! assign previous current density
            psiplasrz = get_psiplasrz()
            psi_ev(:, :, jt) = psiplasrz
! Construct correction
            do j=1, nz2
                do i=1, nr2
                    f_correction = 0.
                    do k=1, nactive
                        f_correction = f_correction + curdiff(k, jt)*greeni(i, j, k)
                    enddo
                    psirz(i, j) = psi_ev(i, j, jt) + f_correction ! total flux
                enddo
            enddo
            do j=1, ntheta
                psicorr(j) = interp2d_psi(rbref_ev(j, jt), zbref_ev(j, jt), Rrect(1: nr), Zrect(1: nz), psirz(1: nr, 1: nz))
            enddo
            x1 = sum(psicorr)/(ntheta + 0.) !average psi on the boundary
! Derivative at ref axis
            bub(1) = interp2d_psi(raxref_ev(jt) - 0.5*dr, zaxref_ev(jt), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(2) = interp2d_psi(raxref_ev(jt) + 0.5*dr, zaxref_ev(jt), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(3) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(4) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            x2 = (bub(2) - bub(1))/dr ! dPsi/dr
            x3 = (bub(4) - bub(3))/dz ! dPsi/dz
            psi_ext_ev(jt) = sum(G_00c(:,jt)*curdiff(:,jt))
            Ffunc_ev(jt) = sigma_B*sum((psicorr - x1)**2) + sigma_axis*(x2**2 + x3**2)
            do i=1, nactive
                Ffunc_ev(jt) = Ffunc_ev(jt) + 0.5*indconduc(i, i)*sigma_coils(i)*sigma_energy*curdiff(i,jt)**2
            enddo

! Calculate F derivative
            do i=1, nactive
                Fderiv((jt-1)*nactive+i) = 2.*(0.5*indconduc(i, i)*sigma_coils(i)*sigma_energy*curdiff(i, jt) + &
                    sigma_B*sum((psicorr - x1)*(G_00(i, 1:ntheta, jt) - G_00c(i, jt))) + &
                    sigma_axis*(x2*G_00r(i, jt) + x3*G_00z(i, jt)) )
            enddo
        enddo !end time loop

! Add lambda contributions
        Ffunc = sum(Ffunc_ev) + lambda*(deltapsiext - (psi_ext_ev(2) - psi_ext_ev(1)))
        do i=1, nactive
            Fderiv(i) = Fderiv(i) + lambda*G_00c(i, 1)
        enddo
        do i=1, nactive
            Fderiv(nactive+i) = Fderiv(nactive+i) - lambda*G_00c(i, 2)
        enddo
        Fderiv(2*nactive+1) = deltapsiext - (psi_ext_ev(2) - psi_ext_ev(1)) ! dF/dlambda

! Calculate new currents
        do i=1, n_evol*nactive+1
            result_vector(i) = result_vector(i) - sum(invmatrix(i, 1: n_evol*nactive+1)*Fderiv)
        enddo
        curdiff(:, 1) = result_vector(1: nactive)
        curdiff(:, 2) = result_vector(nactive+1: n_evol*nactive)
        lambda = result_vector(n_evol*nactive+1)
! Apply limits
        do jt=1, n_evol
            do i=1, nactive
                curdiff(i, jt) = min(curdiff(i, jt), current_limit(i, 1))
                curdiff(i, jt) = max(curdiff(i, jt), current_limit(i, 2))
            enddo
        enddo

! Update plasma current density field
        do jt=1, n_evol
! Construct correction
            do j=1, nz2
                do i=1, nr2
                    f_correction = 0.
                    do k=1, nactive
                        f_correction = f_correction + curdiff(k, jt)*greeni(i, j, k)
                    enddo
                    psirz(i, j) = psi_ev(i, j, jt) + f_correction ! total flux
                enddo
            enddo

! Axis block
            iaxis = iax_ev(jt)
            jaxis = jax_ev(jt)
            call find_new_axis
            iax_ev(jt) = iaxis
            jax_ev(jt) = jaxis

! Boundary block
            n_of_xpoints = nxp_ev(jt)
            if (n_of_xpoints > 0) then
                r_xpoint(1: n_of_xpoints) = rxp_ev(1: n_of_xpoints, jt)
                z_xpoint(1: n_of_xpoints) = zxp_ev(1: n_of_xpoints, jt)
            endif
            call find_psi_boundary
            nxp_ev(jt) = n_of_xpoints
            if (n_of_xpoints > 0) then
                rxp_ev(1: n_of_xpoints, jt) = r_xpoint(1: n_of_xpoints)
                zxp_ev(1: n_of_xpoints, jt) = z_xpoint(1: n_of_xpoints)
            endif

! Current block
            call linterp(psia_ev(:, jt), ffprim_ev(:, jt), nrho, psia_2d, ffp_2d, nrho2d)
            call linterp(psia_ev(:, jt), pprim_ev (:, jt), nrho, psia_2d, ppp_2d, nrho2d)
            ffp_2d = -GP2/mu0*ffp_2d
            ppp_2d = -GP2*1.e-6*ppp_2d
            call new_jrz  ! calculate new right hand side
            j_ev(:, :, jt) = jrz

            i_totev(1: nactive, jt) = curdiff(1: nactive, jt)
        enddo

        temp_err = sum(abs(Fderiv)) !error
        Ffunc_old = Ffunc
        if (temp_err <= err_find_psistab) EXIT
    enddo

    open(32, file='dat/output_timefit.dat')
        write(32, *) 't1 currents: ', i_totev(1: nactive, 1)
        write(32, *) 't2 currents: ', i_totev(1: nactive, 2)
    close(32)

    pause

    end subroutine restab_2_timepoints_evolution_limits

!--------------------------------------------------------------------
    subroutine restab_j_timepoints_evolution_limits_xpoints_boundariz  ! this one does everthing, boundary and isoflux 4 and more isoflux

! Finds active currents from scratch including evolution from time point j-1 to point j

! to calculate forces using this mode, call: call coil_forces(ncoil_blocks, force_R, force_Z, 0)   !the 0 at the end means no plasma contribution. also eddy currents are ignored.

! The cost function should be:

! F = 1/2 * G

! G = sigma_b*sum_b(psib-psibavg)^2 + 
!     sigma_ax*grad(psi)_ax^2 + 
!     sigma_X*sum_E grad(psi)_E^2 + 
!     lambda*(dpsiext - Vloop dt) + 
!     sum_j ( sigma_coils_j*I_j^2 + 0.5*sigma_energy*L_j*I_j^2 + sigma_coils_ref_j*(I_j -I_j_ref)^2 )

    use errors_params, only: err_find_psistab
    use transport2fbe, only: sigma_coils, sigma_coils_ref, sigma_b, sigma_axis, sigma_energy, &
        current_limit, sigma_xpoint, r_xpoint_fit, z_xpoint_fit, &
        n_xpoint_fit, vloop_avg, L_ext, &
        dIp_dt, tau_gseq, time_astra, &
        use_isoflux, n_isoflux, r_isoflux, z_isoflux, which_x_point, &
        voltage_limits_active_coils, cur_init, sigma_isoflux
    use green_function, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix
    use pi_const, only: GP2, GP2_sq
    use numerical_tools, only: linterp

    integer :: i, j, k, j_iter, iax, jax, j_time, ntheta_temp
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, lambda, deltapsiext, psibexto, psibext !lagrange multplier lambda
    double precision, dimension(9) :: bub
    double precision, dimension(ntheta) :: psicorr, rbref, zbref
    double precision, dimension(nactive+1) :: Fderiv
    double precision, dimension(nactive, ntheta) :: G_00, G_002
    double precision, dimension(nactive) :: G_00c, G_00c2, G_00r, G_00z, curref, curnow, curdiff
    double precision, dimension(nactive+1,nactive+1) :: matrix, invmatrix
    double precision, dimension(:, :), allocatable :: G_00xr, G_00xz
    double precision, dimension(:), allocatable :: dummyx, rbndtemp, zbndtemp
    double precision, dimension(nblocks-npassive) :: force_r, force_z
    double precision :: raxref, zaxref, r_norm_ref, sigma_boundary_points(ntheta), boundary_weight
    double precision :: currents_limits_adds(nactive, 2)
    integer, dimension(:), allocatable :: gridpoint_tipe
    character(len=80) :: file_time
    logical :: file_exists
    data j_time/0/
    save j_time

    deltapsiext = tau_gseq*(vloop_avg+L_ext*dIp_dt)   ! jump in psiext

    r_norm_ref = 0.5*(Rrect(1) + Rrect(nr2))
    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhotheta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, Rrect(1), dr)
    jaxis = closest_index(zax, Zrect(1), dz)
    iax = iaxis
    jax = jaxis
    raxref = raxp
    zaxref = zaxp
    rbref(1:ntheta) = rbndp(1:ntheta)
    zbref(1:ntheta) = zbndp(1:ntheta)
    ntheta_temp = ntheta
    nbnd = ntheta
    rbnd(1:nbnd) = rbref(1:ntheta)
    zbnd(1:nbnd) = zbref(1:ntheta)
    allocate(rbndtemp(ntheta),zbndtemp(ntheta))
    
    sigma_boundary_points = 1. ! later everything gets multiplied by sigma_b which is the global weight of the boundary
    
    if (use_isoflux > 0) then
    !the code assumes that isoflux 1 is raus, isoflux 2 is ztop, isoflux 3 is rin, isoflux 4 is zbot. if one or more than one are x points is given in which_x_point
        ntheta_temp = n_isoflux
        n_xpoint_fit = n_isoflux
        if (n_xpoint_fit > 0) then
            allocate(G_00xr(nactive, n_xpoint_fit))
            allocate(G_00xz(nactive, n_xpoint_fit))
            allocate(dummyx(n_xpoint_fit))
        endif
        allocate(gridpoint_tipe(n_isoflux))
        rbref(1:ntheta_temp) = r_isoflux(1:ntheta_temp)          
        zbref(1:ntheta_temp) = z_isoflux(1:ntheta_temp)          
        r_xpoint_fit(1:ntheta_temp)=rbref(1:ntheta_temp)
        z_xpoint_fit(1:ntheta_temp)=zbref(1:ntheta_temp)
        gridpoint_tipe = -1
        do i = 1, n_isoflux
            if (which_x_point(i) == -1) gridpoint_tipe(i) = -1
            if (which_x_point(i) == 0)  gridpoint_tipe(i) =  0
            if (which_x_point(i) == 1)  gridpoint_tipe(i) =  1
            if (which_x_point(i) == 2)  gridpoint_tipe(i) =  2
        enddo
	sigma_boundary_points(1:n_isoflux) = sigma_isoflux(1:n_isoflux)
    else
        if (n_xpoint_fit > 0) then
            allocate(G_00xr(nactive, n_xpoint_fit))
            allocate(G_00xz(nactive, n_xpoint_fit))
            allocate(dummyx(n_xpoint_fit))
            allocate(gridpoint_tipe(n_xpoint_fit))
        endif
    endif
    
    boundary_weight = sum(sigma_boundary_points(1:ntheta_temp))

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)
    if (j_time == 0) curnow(1:nactive) = cur_init(1:nactive)
    if (j_time  > 0) curnow(1:nactive) = curconduc(1:nactive)
    curref(1:nactive) = curnow(1:nactive)
    curdiff = 0.

! Calculate additional current limits
    do i=1, nactive
         currents_limits_adds(i, 1) = curref(i) + tau_gseq*(voltage_limits_active_coils(i, 1) - resconduc(i, i)*curref(i))/indconduc(i, i)
         currents_limits_adds(i, 2) = curref(i) + tau_gseq*(voltage_limits_active_coils(i, 2) - resconduc(i,i)*curref(i))/indconduc(i, i)
    enddo

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.
    lambda = 0.

    do j=1, nactive
        do k=1, ntheta_temp
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
            G_002(j, k) = sigma_boundary_points(k)*interp2d_psi(rbref(k), zbref(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:ntheta_temp))/(0. + ntheta_temp)
        G_00c2(j) = sum(G_002(j, 1:ntheta_temp))/boundary_weight
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
        if (n_xpoint_fit > 0) then
            do k=1, n_xpoint_fit
                bub(1) = interp2d_psi(r_xpoint_fit(k) - dr/2., z_xpoint_fit(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
                bub(2) = interp2d_psi(r_xpoint_fit(k) + dr/2., z_xpoint_fit(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
                bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(k) - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
                bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(k) + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
                SELECT CASE(gridpoint_tipe(k))
                CASE(-1)
                    G_00xr(j,k) = 0.
                    G_00xz(j,k) = 0.
                CASE(0)
                    G_00xr(j,k) = 0.
                    G_00xz(j,k) = (bub(4) - bub(3))/dz
                CASE(1)
                    G_00xr(j,k) = (bub(2) - bub(1))/dr
                    G_00xz(j,k) = 0.
                CASE(2)
                    G_00xr(j,k) = (bub(2) - bub(1))/dr
                    G_00xz(j,k) = (bub(4) - bub(3))/dz
                END SELECT
            enddo
        endif
    enddo

    psibexto = sum(G_00c2*curnow)

    do j=1, nactive
        do i=1, nactive
            if (i == j) matrix(i, j) = matrix(i, j) + sigma_energy*indconduc(i, i) + 2.*(sigma_coils(i) + sigma_coils_ref(i))
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_002(i, 1:ntheta_temp) - G_00c2(i))*(G_002(j, 1:ntheta_temp) - G_00c2(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j))
            if (n_xpoint_fit > 0) then
                dummyx(1:n_xpoint_fit) = G_00xr(i, 1:n_xpoint_fit)*G_00xr(j, 1:n_xpoint_fit) + &
                                         G_00xz(i, 1:n_xpoint_fit)*G_00xz(j, 1:n_xpoint_fit)
                matrix(i, j) = matrix(i, j) + 2.*sigma_xpoint*sum(dummyx)
            endif
        enddo
    enddo

! Borders with lambda
    do i=1, nactive
        matrix(nactive+1, i) = G_00c2(i)
        matrix(i, nactive+1) = G_00c2(i)
    enddo
    matrix(nactive+1, nactive+1) = 0.

! Calculate inverse
    if (j_time == 0) then
        invmatrix(1:nactive, 1:nactive) = inv_matrix(matrix(1:nactive, 1:nactive), nactive)
    else
        invmatrix = inv_matrix(matrix, nactive+1)
    endif

! Iteration to find currents
    do j_iter=1, 300000
        psiplasrz = get_psiplasrz()
! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nactive
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        do j=1, ntheta_temp
            psicorr(j) = sigma_boundary_points(j)*interp2d_psi(rbref(j), zbref(j), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr(1:ntheta_temp))/boundary_weight !average psi on the boundary

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz

        Ffunc = sigma_B*sum((psicorr(1:ntheta_temp) - x1)**2) + sigma_axis*(x2**2 + x3**2)
        do i=1, nactive
            Ffunc = Ffunc+0.5*indconduc(i, i)*sigma_energy*curnow(i)**2 + sigma_coils(i)*curnow(i)**2 + &
                    sigma_coils_ref(i)*curdiff(i)**2
        enddo

        if (n_xpoint_fit > 0) then
            do k=1, n_xpoint_fit
                bub(1) = interp2d_psi(r_xpoint_fit(k) - 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                bub(2) = interp2d_psi(r_xpoint_fit(k) + 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                SELECT CASE(gridpoint_tipe(k))
                CASE(-1)
                    x2 = 0.
                    x3 = 0.
                CASE(0)
                    x2 = 0.
                    x3 = (bub(4) - bub(3))/dz
                CASE(1)
                    x2 = (bub(2) - bub(1))/dr
                    x3 = 0.
                CASE(2)
                    x2 = (bub(2) - bub(1))/dr
                    x3 = (bub(4) - bub(3))/dz
                END SELECT
                Ffunc = Ffunc +  sigma_xpoint*(x2**2 + x3**2)
             enddo
        endif

        psibext = sum(G_00c2*curnow)

        if (j_time > 0) Ffunc = Ffunc + lambda*(psibext - psibexto + deltapsiext/GP2)

! Calculate F derivative
        Fderiv = 0.
        do i=1, nactive
            Fderiv(i) = 2.*(0.5*indconduc(i,i)*sigma_energy*curnow(i) + &
                sigma_coils_ref(i)*curdiff(i) + &
                sigma_coils(i)*curnow(i) + sigma_B*sum((psicorr(1:ntheta_temp) - x1)*(G_002(i, 1:ntheta_temp) - G_00c2(i))))
            bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            x2 = (bub(2) - bub(1))/dr ! dPsi/dr
            x3 = (bub(4) - bub(3))/dz ! dPsi/dz
            Fderiv(i) = Fderiv(i) + 2.*sigma_axis*(x2*G_00r(i) + x3*G_00z(i))
            if (n_xpoint_fit > 0) then
                do k=1, n_xpoint_fit
                    bub(1) = interp2d_psi(r_xpoint_fit(k) - 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                    bub(2) = interp2d_psi(r_xpoint_fit(k) + 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                    bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                    bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                    SELECT CASE(gridpoint_tipe(k))
                    CASE(-1)
                        x2 = 0.
                        x3 = 0.
                    CASE(0)
                        x2 = 0.
                        x3 = (bub(4) - bub(3))/dz
                    CASE(1)
                        x2 = (bub(2) - bub(1))/dr
                        x3 = 0.
                    CASE(2)
                        x2 = (bub(2) - bub(1))/dr
                        x3 = (bub(4) - bub(3))/dz
                    END SELECT
                    Fderiv(i) = Fderiv(i) +  2.*sigma_xpoint*(x2*G_00xr(i, k) + x3*G_00xz(i, k))
                enddo
            endif
        enddo

        if (j_time > 0) Fderiv(1:nactive) = Fderiv(1:nactive) + lambda*G_00c2(1:nactive)
        if (j_time > 0) Fderiv(nactive+1) = (psibext-psibexto+deltapsiext/GP2)

! Calculate new currents

        do i=1, nactive
           if (j_time == 0) then
               curnow(i) = curnow(i) - sum(invmatrix(i, 1:nactive)*Fderiv(1:nactive))
           else
               curnow(i) = curnow(i) - sum(invmatrix(i, 1:nactive+1)*Fderiv(1:nactive+1))
           endif
        enddo
        lambda = lambda - sum(invmatrix(nactive+1, 1:nactive+1)*Fderiv(1:nactive+1))

! cut new currents to limits due to current
        do i=1, nactive
            curnow(i) = min(curnow(i), current_limit(i, 1))
            curnow(i) = max(curnow(i), current_limit(i, 2))
        enddo
! cut new currents to limits due to voltage if jtime> 0
        if (j_time > 0) then 
            do i=1, nactive
                curnow(i) = min(curnow(i), currents_limits_adds(i, 1))
                curnow(i) = max(curnow(i), currents_limits_adds(i, 2))
            enddo
        endif

        curdiff = curnow - curref

! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nactive
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        call find_new_axis
        call find_psi_boundary
        call new_jrz  ! calculate new right hand side

        temp_err = sum(abs(Fderiv))
        Ffunc_old = Ffunc

        if (j_iter > 15000) EXIT
        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nactive
        curconduc(i) = curnow(i)
    enddo

    write(file_time,'(1A, I0)') 'dat/file_fit_currents_output_', j_time

    inquire(FILE=file_time, EXIST=file_exists)
    if (file_exists) call execute_command_line('rm '//trim(file_time), wait=.true.)

! Evaluate forces
    call coil_forces(nblocks - npassive, force_R, force_Z, 1)

    call psi_mutual_effect_conductors(psiplasmatoconduc)

    call estimate_boundary_to_pbe(rbndtemp, zbndtemp, ntheta)

    open(32, file = file_time)
    write(32, *) nactive, curconduc(1:nactive)*1e3, -GP2*psibext, rax, zax, &
        tau_gseq, vloop_avg, L_ext*dIp_dt, &
        nr2, nz2, psirz(1:nr2, 1:nz2), nblocks - npassive, force_R, force_Z, &
        psibnd, psiaxis, indconduc(1:nactive, 1:nactive), resconduc(1:nactive, 1:nactive), &
        nlimiter, limiterR(1:nlimiter), limiterZ(1:nlimiter), Rrect(1), Rrect(nr2), Zrect(1), Zrect(nz2), ntheta_temp, &
        rbref(1:ntheta_temp), zbref(1:ntheta_temp), time_astra, tau_gseq, GP2*psiplasmatoconduc(1:nactive), &  !saved in kA
        ntheta, rbndtemp(1:ntheta), zbndtemp(1:ntheta), rbndp(1:ntheta), zbndp(1:ntheta)
    close(32)

!update time
    j_time = j_time +1

    if (n_xpoint_fit > 0) then
        deallocate(G_00xr)
        deallocate(G_00xz)
        deallocate(dummyx)
    endif
    if (allocated(gridpoint_tipe)) deallocate(gridpoint_tipe)
    deallocate(rbndtemp, zbndtemp)

    end subroutine restab_j_timepoints_evolution_limits_xpoints_boundariz
  
!--------------------------------------------------------------------
    subroutine restab_1_timepoint_limits_xpoints_boundariz  ! this one does everthing, boundary and isoflux 4 and more isoflux

! Finds active currents from scratch for a fixed time point

    use errors_params, only: err_find_psistab
    use transport2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy, &
        sigma_xpoint, r_xpoint_fit, z_xpoint_fit, &
        n_xpoint_fit, vloop_avg, L_ext, &
        dIp_dt, tau_gseq, time_astra, &
        use_isoflux, n_isoflux, r_isoflux, z_isoflux, which_x_point, &
        cur_init
    use green_function, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix
    use pi_const, only: GP2, GP2_sq
    use numerical_tools, only: linterp

    integer :: i, j, k, j_iter, iax, jax, ntheta_temp
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, psibexto, psibext
    double precision, dimension(9) :: bub
    double precision, dimension(ntheta) :: psicorr, rbref, zbref
    double precision, dimension(nactive) :: Fderiv
    double precision, dimension(nactive, ntheta) :: G_00
    double precision, dimension(nactive) :: G_00c, G_00r, G_00z, curref, curnow, curdiff
    double precision, dimension(nactive,nactive) :: matrix, invmatrix
    double precision, dimension(:, :), allocatable :: G_00xr, G_00xz
    double precision, dimension(:), allocatable :: dummyx, rbndtemp, zbndtemp
    double precision, dimension(nblocks-npassive) :: force_r, force_z
    double precision :: raxref, zaxref, r_norm_ref
    integer, dimension(:), allocatable :: gridpoint_tipe
    character(len=80) :: file_time
    logical :: file_exists

    r_norm_ref = 0.5*(Rrect(1) + Rrect(nr2))
    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhotheta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, Rrect(1), dr)
    jaxis = closest_index(zax, Zrect(1), dz)
    iax = iaxis
    jax = jaxis
    raxref = raxp
    zaxref = zaxp
    rbref(1:ntheta) = rbndp(1:ntheta)
    zbref(1:ntheta) = zbndp(1:ntheta)
    ntheta_temp = ntheta
    nbnd = ntheta
    rbnd(1:nbnd) = rbref(1:ntheta)
    zbnd(1:nbnd) = zbref(1:ntheta)
    allocate(rbndtemp(ntheta),zbndtemp(ntheta))

    if (use_isoflux > 0) then
    !the code assumes that isoflux 1 is raus, isoflux 2 is ztop, isoflux 3 is rin, isoflux 4 is zbot. if one or more than one are x points is given in which_x_point
        ntheta_temp = n_isoflux
        n_xpoint_fit = n_isoflux
        if (n_xpoint_fit > 0) then
            allocate(G_00xr(nactive, n_xpoint_fit))
            allocate(G_00xz(nactive, n_xpoint_fit))
            allocate(dummyx(n_xpoint_fit))
        endif
        allocate(gridpoint_tipe(n_isoflux))
        rbref(1:ntheta_temp) = r_isoflux(1:ntheta_temp)          
        zbref(1:ntheta_temp) = z_isoflux(1:ntheta_temp)          
        r_xpoint_fit(1:ntheta_temp)=rbref(1:ntheta_temp)
        z_xpoint_fit(1:ntheta_temp)=zbref(1:ntheta_temp)
        gridpoint_tipe = -1
        do i = 1, n_isoflux
            if (which_x_point(i) == -1) gridpoint_tipe(i) = -1
            if (which_x_point(i) == 0)  gridpoint_tipe(i) =  0
            if (which_x_point(i) == 1)  gridpoint_tipe(i) =  1
            if (which_x_point(i) == 2)  gridpoint_tipe(i) =  2
        enddo
    else
        if (n_xpoint_fit > 0) then
            allocate(G_00xr(nactive, n_xpoint_fit))
            allocate(G_00xz(nactive, n_xpoint_fit))
            allocate(dummyx(n_xpoint_fit))
            allocate(gridpoint_tipe(n_xpoint_fit))
        endif
    endif

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)
		curnow(1:nactive) = cur_init(1:nactive)
    curref(1:nactive) = curnow(1:nactive)
    curdiff = 0.

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.

    do j=1, nactive
        do k=1, ntheta_temp
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:ntheta_temp))/(0. + ntheta_temp)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
        if (n_xpoint_fit > 0) then
            do k=1, n_xpoint_fit
                bub(1) = interp2d_psi(r_xpoint_fit(k) - dr/2., z_xpoint_fit(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
                bub(2) = interp2d_psi(r_xpoint_fit(k) + dr/2., z_xpoint_fit(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
                bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(k) - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
                bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(k) + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
                SELECT CASE(gridpoint_tipe(k))
                CASE(-1)
                    G_00xr(j,k) = 0.
                    G_00xz(j,k) = 0.
                CASE(0)
                    G_00xr(j,k) = 0.
                    G_00xz(j,k) = (bub(4) - bub(3))/dz
                CASE(1)
                    G_00xr(j,k) = (bub(2) - bub(1))/dr
                    G_00xz(j,k) = 0.
                CASE(2)
                    G_00xr(j,k) = (bub(2) - bub(1))/dr
                    G_00xz(j,k) = (bub(4) - bub(3))/dz
                END SELECT
            enddo
        endif
    enddo

    psibexto = sum(G_00c*curnow)

    do j=1, nactive
        do i=1, nactive
            if (i == j) matrix(i, j) = matrix(i, j) + sigma_energy/sigma_coils(i)*indconduc(i, i) + 2.*sigma_coils(i)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:ntheta_temp) - G_00c(i))*(G_00(j, 1:ntheta_temp) - G_00c(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j))
            if (n_xpoint_fit > 0) then
                dummyx(1:n_xpoint_fit) = G_00xr(i, 1:n_xpoint_fit)*G_00xr(j, 1:n_xpoint_fit) + &
                                         G_00xz(i, 1:n_xpoint_fit)*G_00xz(j, 1:n_xpoint_fit)
                matrix(i, j) = matrix(i, j) + 2.*sigma_xpoint*sum(dummyx)
            endif
        enddo
    enddo

! Calculate inverse
        invmatrix(1:nactive, 1:nactive) = inv_matrix(matrix(1:nactive, 1:nactive), nactive)

! Iteration to find currents
    do j_iter=1, 300000
        psiplasrz = get_psiplasrz()
! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nactive
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        do j=1, ntheta_temp
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr(1:ntheta_temp))/(ntheta_temp + 0.) !average psi on the boundary

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz

        Ffunc = sigma_B*sum((psicorr(1:ntheta_temp) - x1)**2) + sigma_axis*(x2**2 + x3**2)
        do i=1, nactive
            Ffunc = Ffunc+0.5*indconduc(i, i)*sigma_energy/sigma_coils(i)*curnow(i)**2 + sigma_coils(i)*curdiff(i)**2
        enddo

        if (n_xpoint_fit > 0) then
            do k=1, n_xpoint_fit
                bub(1) = interp2d_psi(r_xpoint_fit(k) - 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                bub(2) = interp2d_psi(r_xpoint_fit(k) + 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                SELECT CASE(gridpoint_tipe(k))
                CASE(-1)
                    x2 = 0.
                    x3 = 0.
                CASE(0)
                    x2 = 0.
                    x3 = (bub(4) - bub(3))/dz
                CASE(1)
                    x2 = (bub(2) - bub(1))/dr
                    x3 = 0.
                CASE(2)
                    x2 = (bub(2) - bub(1))/dr
                    x3 = (bub(4) - bub(3))/dz
                END SELECT
                Ffunc = Ffunc +  sigma_xpoint*(x2**2 + x3**2)
             enddo
        endif

        psibext = sum(G_00c*curnow)

! Calculate F derivative
        Fderiv = 0.
        do i=1, nactive
            Fderiv(i) = 2.*(0.5*indconduc(i,i)*sigma_energy/sigma_coils(i)*curnow(i) + &
                sigma_coils(i)*curdiff(i) + sigma_B*sum((psicorr(1:ntheta_temp) - x1)*(G_00(i, 1:ntheta_temp) - G_00c(i))))
            bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
            x2 = (bub(2) - bub(1))/dr ! dPsi/dr
            x3 = (bub(4) - bub(3))/dz ! dPsi/dz
            Fderiv(i) = Fderiv(i) + 2.*sigma_axis*(x2*G_00r(i) + x3*G_00z(i))
            if (n_xpoint_fit > 0) then
                do k=1, n_xpoint_fit
                    bub(1) = interp2d_psi(r_xpoint_fit(k) - 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                    bub(2) = interp2d_psi(r_xpoint_fit(k) + 0.5*dr, z_xpoint_fit(K), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                    bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                    bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
                    SELECT CASE(gridpoint_tipe(k))
                    CASE(-1)
                        x2 = 0.
                        x3 = 0.
                    CASE(0)
                        x2 = 0.
                        x3 = (bub(4) - bub(3))/dz
                    CASE(1)
                        x2 = (bub(2) - bub(1))/dr
                        x3 = 0.
                    CASE(2)
                        x2 = (bub(2) - bub(1))/dr
                        x3 = (bub(4) - bub(3))/dz
                    END SELECT
                    Fderiv(i) = Fderiv(i) +  2.*sigma_xpoint*(x2*G_00xr(i, k) + x3*G_00xz(i, k))
                enddo
            endif
        enddo

! Calculate new currents

        do i=1, nactive
           curnow(i) = curnow(i) - sum(invmatrix(i, 1:nactive)*Fderiv(1:nactive))
        enddo

! cut new currents to limits due to current
! cut new currents to limits due to voltage if jtime> 0

        curdiff = curnow - curref

! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nactive
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        call find_new_axis
        call find_psi_boundary
        call new_jrz  ! calculate new right hand side

        temp_err = sum(abs(Fderiv))
        Ffunc_old = Ffunc

        if (j_iter > 15000) EXIT
        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nactive
        curconduc(i) = curnow(i)
    enddo

    write(file_time,'(1A)') 'dat/file_fit_currents_output_fixed_time'

    inquire(FILE=file_time, EXIST=file_exists)
    if (file_exists) call execute_command_line('rm '//trim(file_time), wait=.true.)

! Evaluate forces
    call coil_forces(nblocks - npassive, force_R, force_Z, 1)

    call psi_mutual_effect_conductors(psiplasmatoconduc)

    call estimate_boundary_to_pbe(rbndtemp, zbndtemp, ntheta)

    open(32, file = file_time)
    write(32, *) nactive, curconduc(1:nactive)*1e3, -GP2*psibext, rax, zax, &
        tau_gseq, vloop_avg, L_ext*dIp_dt, &
        nr2, nz2, psirz(1:nr2, 1:nz2), nblocks - npassive, force_R, force_Z, &
        psibnd, psiaxis, indconduc(1:nactive, 1:nactive), resconduc(1:nactive, 1:nactive), &
        nlimiter, limiterR(1:nlimiter), limiterZ(1:nlimiter), Rrect(1), Rrect(nr2), Zrect(1), Zrect(nz2), ntheta_temp, &
        rbref(1:ntheta_temp), zbref(1:ntheta_temp), time_astra, tau_gseq, GP2*psiplasmatoconduc(1:nactive), &  !saved in kA
        ntheta, rbndtemp(1:ntheta), zbndtemp(1:ntheta)
    close(32)


    if (n_xpoint_fit > 0) then
        deallocate(G_00xr)
        deallocate(G_00xz)
        deallocate(dummyx)
    endif
    if (allocated(gridpoint_tipe)) deallocate(gridpoint_tipe)
    deallocate(rbndtemp, zbndtemp)
    write(*,*) 'stopping'
    stop

    end subroutine restab_1_timepoint_limits_xpoints_boundariz
  
!--------------------------------------------------------------------
    subroutine diagnose(filename)  !call it in sbr/assign_geom.f90 in astra for example, it has access to this module.

    use metric_coefficients_pbe, only: dator

    integer, parameter :: unit=32
    character(*), intent(in) :: filename

    open(unit, file=trim(filename))
!first PBE stuff
        write(unit, *) nrho, ntheta
        write(unit, *) pprime, ffprime, psigrida  ! pprim, ffprim, psigrid
        write(unit, *) rpol, zpol, jrhotheta, psirhotheta  !R, Z, jrhotheta, PSI
        write(unit, *) dator !area elements
!now free boundary
        if (allocated(Rrect)) then
            write(unit, *) nr2, nz2, Rrect, Zrect !grid
            write(unit, *) psiextrz, psiplasrz, psirz !psivacuum, psiplasma, psitotal maps [radiants]
            write(unit, *) jrz !current densiy, area elements
            write(unit, *) nconduc, r_cond, z_cond !conductors positions
            write(unit, *) curconduc, voltage, psiplasmatoconduc !currents, voltages, psiplasma mutual induct
            write(unit, *) nlimiter, limiterr, limiterz !limiter
        else
            write(unit, *) '-1'
        endif
    close(unit)

    end subroutine diagnose
  
!--------------------------------------------------------------------
    subroutine restab_F_function_full_currents
! Finds active currents from scratch. Passive currents are given.

    use errors_params, only: err_find_psistab
    use transport2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy
    use green_function, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix

    integer :: i, j, k, j_iter, iax, jax
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, raxref, zaxref
    double precision, dimension(9) :: bub
    double precision, dimension(500) :: rbref, zbref
    double precision, dimension(nactive) :: G_00c, G_00r, G_00z, &
        Fderiv, curref, curnow, curdiff
    double precision, dimension(ntheta) :: psicorr
    double precision, dimension(nactive, ntheta) :: G_00
    double precision, dimension(nactive, nactive) :: matrix, invmatrix

    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhotheta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, Rrect(1), dr)
    jaxis = closest_index(zax, Zrect(1), dz)
    iax = iaxis
    jax = jaxis

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    curref(1:nactive) = 0.
    curnow(1:nactive) = 0.
    curconduc(1:nactive) = 0.
    curdiff = 0.

    raxref = raxp
    zaxref = zaxp
    rbref(1:ntheta) = rbndp(1:ntheta)
    zbref(1:ntheta) = zbndp(1:ntheta)

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.

    do j=1, nactive
        do k=1, ntheta
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:ntheta))/(0. + ntheta)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
    enddo

    do j=1, nactive
        do i=1, nactive
            if (i == j) matrix(i, j) = matrix(i, j) + sigma_energy*sigma_coils(i)*indconduc(i,i)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:ntheta) - G_00c(i))*(G_00(j, 1:ntheta) - G_00c(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j))
        enddo
    enddo

! Calculate inverse
    invmatrix = inv_matrix(matrix, nactive)

! Iteration to find currents
    do j_iter=1, 300000 
        if (j_iter > 150) stop
        psiplasrz = get_psiplasrz()
! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nactive
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        do j=1, ntheta
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(ntheta + 0.) !average psi on the boundary

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz

        Ffunc = sigma_B*sum((psicorr - x1)**2) + sigma_axis*(x2**2 + x3**2)
        do i=1,nactive
            Ffunc = Ffunc + 0.5*indconduc(i, i)*sigma_energy*sigma_coils(i)*curdiff(i)**2
        enddo

! Calculate F derivative
        do i=1, nactive
            Fderiv(i) = 2.*(0.5*indconduc(i,i)*sigma_energy*sigma_coils(i)*curdiff(i) + &
                sigma_B*sum((psicorr - x1)*(G_00(i, 1:ntheta) - G_00c(i))) + &
                sigma_axis*(x2*G_00r(i) + x3*G_00z(i)) )
        enddo

! Calculate new currents

        do i=1, nactive
            curnow(i) = curnow(i) - sum(invmatrix(i, 1:nactive)*Fderiv)
        enddo
        curdiff = curnow - curref

! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nactive
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        call find_new_axis
        call find_psi_boundary
        call new_jrz  ! calculate new right hand side

        temp_err = sum(abs(Fderiv))
        Ffunc_old = Ffunc

        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nactive
        curconduc(i) = curnow(i)
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis
    call find_psi_boundary
    call new_jrz  ! calculate new right hand side

    end subroutine restab_F_function_full_currents

!--------------------------------------------------------------------
    subroutine restab_F_function_full_currents_limits
! Finds active currents from scratch. Passive currents are given.

    use errors_params, only: err_find_psistab
    use transport2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy, &
      current_limit
    use green_function, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix

    integer :: i, j, k, j_iter, iax, jax
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, raxref, zaxref
    double precision, dimension(9) :: bub
    double precision, dimension(500) :: rbref, zbref
    double precision, dimension(nactive) :: G_00c, G_00r, G_00z, &
        Fderiv, curref, curnow, curdiff
    double precision, dimension(ntheta) :: psicorr
    double precision, dimension(nactive, ntheta) :: G_00
    double precision, dimension(nactive, nactive) :: matrix, invmatrix

    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhotheta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, Rrect(1), dr)
    jaxis = closest_index(zax, Zrect(1), dz)
    iax = iaxis
    jax = jaxis

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    curref(1:nactive) = 0.
    curnow(1:nactive) = 0.
    curconduc(1:nactive) = 0.
    curdiff = 0.

    raxref = raxp
    zaxref = zaxp
    rbref(1:ntheta) = rbndp(1:ntheta)
    zbref(1:ntheta) = zbndp(1:ntheta)

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.

    do j=1, nactive
        do k=1, ntheta
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:ntheta))/(0. + ntheta)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
    enddo

    do j=1, nactive
        do i=1, nactive
            if (i == j) matrix(i, j) = matrix(i, j) + sigma_energy*sigma_coils(i)*indconduc(i,i)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:ntheta) - G_00c(i))*(G_00(j, 1:ntheta) - G_00c(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j))
       enddo
    enddo

! Calculate inverse
    invmatrix = inv_matrix(matrix, nactive)
! Iteration to find currents
    do j_iter=1, 300000
        if (j_iter > 150) stop
        psiplasrz = get_psiplasrz()
! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nactive
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        do j=1, ntheta
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(ntheta + 0.) !average psi on the boundary

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz

        Ffunc=sigma_B*sum((psicorr - x1)**2) + sigma_axis*(x2**2 + x3**2)
        do i=1, nactive
            Ffunc = Ffunc + 0.5*indconduc(i,i)*sigma_energy*sigma_coils(i)*curdiff(i)**2
        enddo

! Calculate F derivative
        do i=1, nactive
            Fderiv(i) = 2.*(0.5*indconduc(i,i)*sigma_energy*sigma_coils(i)*curdiff(i) + &
                sigma_B*sum((psicorr - x1)*(G_00(i, 1:ntheta) - G_00c(i))) + &
                sigma_axis*(x2*G_00r(i) + x3*G_00z(i)) )
        enddo

! Calculate new currents

        do i=1, nactive
            curnow(i) = curnow(i) - sum(invmatrix(i, 1:nactive)*Fderiv)
        enddo
        curdiff = curnow - curref
! cut new currents to limits
        do i=1, nactive
            curdiff(i) = min(curdiff(i), current_limit(i, 1))
            curdiff(i) = max(curdiff(i), current_limit(i, 2))
        enddo

! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nactive
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        call find_new_axis
        call find_psi_boundary
        call new_jrz  ! calculate new right hand side

        temp_err = sum(abs(Fderiv))
        Ffunc_old = Ffunc

        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nactive
        curconduc(i) = curnow(i)
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis
    call find_psi_boundary
    call new_jrz  ! calculate new right hand side

    end subroutine restab_F_function_full_currents_limits

!--------------------------------------------------------------------
    subroutine restab_F_function_full_currents_forces
! Finds active currents from scratch. Passive currents are given. Forces get minimized too.

! how to include forces???

    use errors_params, only: err_find_psistab
    use transport2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy
    use green_function, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix

    integer :: i, j, k, j_iter, iax, jax
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, raxref, zaxref
    double precision, dimension(9) :: bub
    double precision, dimension(500) :: rbref, zbref
    double precision, dimension(nactive) :: G_00c, G_00r, G_00z, &
        Fderiv, curref, curnow, curdiff
    double precision, dimension(ntheta) :: psicorr
    double precision, dimension(nactive, ntheta) :: G_00
    double precision, dimension(nactive, nactive) :: matrix, invmatrix

    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhotheta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, Rrect(1), dr)
    jaxis = closest_index(zax, Zrect(1), dz)
    iax = iaxis
    jax = jaxis

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    curref(1:nactive) = 0.
    curnow(1:nactive) = 0.
    curconduc(1:nactive) = 0.
    curdiff = 0.

    raxref = raxp
    zaxref = zaxp
    rbref(1:ntheta) = rbndp(1:ntheta)
    zbref(1:ntheta) = zbndp(1:ntheta)

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.

    do j=1, nactive
        do k=1, ntheta
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:ntheta))/(0. + ntheta)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
    enddo

    do j=1, nactive
        do i=1, nactive
            if (i == j) matrix(i, j) = matrix(i, j) + sigma_energy*sigma_coils(i)*indconduc(i,i)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:ntheta) - G_00c(i))*(G_00(j, 1:ntheta) - G_00c(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j))
       enddo
    enddo

! Calculate inverse
    invmatrix = inv_matrix(matrix, nactive)
! Iteration to find currents
    do j_iter=1, 300000
        if (j_iter > 150) stop
        psiplasrz = get_psiplasrz()
! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nactive
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction ! total flux
            enddo
        enddo

        do j=1, ntheta
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(ntheta + 0.) ! average psi on the boundary

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz

        Ffunc = sigma_B*sum((psicorr - x1)**2) + sigma_axis*(x2**2 + x3**2)
        do i=1, nactive
            Ffunc = Ffunc + 0.5*indconduc(i, i)*sigma_energy*sigma_coils(i)*curdiff(i)**2
        enddo

! Calculate F derivative
        do i=1, nactive
            Fderiv(i) = 2.*(0.5*indconduc(i,i)*sigma_energy*sigma_coils(i)*curdiff(i) + &
                sigma_B*sum((psicorr - x1)*(G_00(i, 1:ntheta) - G_00c(i))) + &
                sigma_axis*(x2*G_00r(i) + x3*G_00z(i)) )
        enddo

! Calculate new currents

        do i=1, nactive
            curnow(i) = curnow(i) - sum(invmatrix(i, 1:nactive)*Fderiv)
        enddo
        curdiff = curnow - curref

! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nactive
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k)
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        call find_new_axis
        call find_psi_boundary
        call new_jrz  ! calculate new right hand side

        temp_err = sum(abs(Fderiv))
        Ffunc_old = Ffunc

        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nactive
        curconduc(i) = curnow(i)
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis
    call find_psi_boundary
    call new_jrz  ! calculate new right hand side

    end subroutine restab_F_function_full_currents_forces

!--------------------------------------------------------------------
    subroutine restab_boundary_with_fourier_wall !not working well

    use errors_params, only: err_find_psistab
    use transport2fbe, only: n_fourier_restab_boundary
    use green_function, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi

    integer :: i, j, k, j_iter, iax, jax, info
    double precision :: temp_err, curr, f_correction, x1, psibt0, psibt1
    double precision, dimension(9) :: bub
    double precision, dimension(npassive) :: anglr
    double precision, dimension(ntheta) :: psicorr
    double precision, dimension(n_fourier_restab_boundary) :: S_00, C_00
    double precision, dimension(4*ntheta*n_fourier_restab_boundary) :: work
    double precision, dimension(ntheta, n_fourier_restab_boundary) :: G_00c, G_00s
    double precision, dimension(ntheta, 2*n_fourier_restab_boundary) :: matrix

! First, initialized initial guess coming from prescribed boundary current density: jrhotheta
    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, Rrect(1), dr)
    jaxis = closest_index(zax, Zrect(1), dz)
    iax = iaxis
    jax = jaxis

    G_00c = 0.
    G_00s = 0.

! Evaluate coils' quantities
    do i=1, npassive
        anglr(i) = ATAN2(z_cond(nactive + i) - zax, r_cond(nactive + i) - rax)
! Find true axis
        do j=1, ntheta
            bub(2) = interp2d_psi(rbndp(j), zbndp(j), Rrect(1:nr), Zrect(1:nz), greeni(1:nr, 1:nz, nactive + i))
            do k=1, n_fourier_restab_boundary
                G_00c(j, k) = G_00c(j, k) + cos(k*anglr(i))*bub(2)
                G_00s(j, k) = G_00s(j, k) + sin(k*anglr(i))*bub(2)
            enddo
        enddo
    enddo

    do j=1, ntheta
        do k=1, n_fourier_restab_boundary
            matrix(j, k) = G_00c(j, k)
            matrix(j, n_fourier_restab_boundary + k) = G_00s(j, k)
        enddo
    enddo

    psibt0 = 1000.
    psibt1 = 1000.

    do j_iter=1, 30
        psiplasrz = get_psiplasrz()
        call compound_psi
        do j=1, ntheta
            psicorr(j) = interp2d_psi(rbndp(j), zbndp(j), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(ntheta + 0.)
        psibt1  = x1
        psicorr = x1 - psicorr

! Least square fit solution

        call dgels('N', ntheta, n_fourier_restab_boundary*2, 1, matrix, ntheta, &
            psicorr, ntheta, WORK, 2*(ntheta)*n_fourier_restab_boundary*2, INFO)

! Construct correction
        call compound_psi
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, n_fourier_restab_boundary
                    f_correction = f_correction +  &
                        psicorr(k)                            *sum(greeni(i, j, nactive + 1:nactive + npassive)*cos(k*anglr(1:npassive))) + &
                        psicorr(n_fourier_restab_boundary + k)*sum(greeni(i, j, nactive + 1:nactive + npassive)*sin(k*anglr(1:npassive)))
                enddo
                psirz(i, j) = psirz(i, j) + f_correction !total flux
            enddo
        enddo

        call find_new_axis
        call find_psi_boundary
        call new_jrz  ! calculate new right hand side

        temp_err = abs(psibt0 - psibt1)
        psibt0 = psibt1
        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, npassive
        do j=1, n_fourier_restab_boundary
            curconduc(nactive + i) = curconduc(nactive + i) + C_00(j)*cos(j*anglr(i)) + S_00(j)*sin(j*anglr(i))
        enddo
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis
    call find_psi_boundary
    call new_jrz  ! calculate new right hand side

    stop
    end subroutine restab_boundary_with_fourier_wall

!--------------------------------------------------------------------
    subroutine restab_axis_with_fourier_wall

    use errors_params, only: err_find_psistab
    use green_function, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, least_square_biquad

    integer :: i, j, j_iter, iax, jax
    double precision :: curr, dum1, psistab1o, psistab2o, delr, delz, &
        S_00r, C_00r, S_00z, C_00z, temp_err
    double precision, dimension(5) :: ddipsi
    double precision, dimension(9) :: bub, xub, yub
    double precision, dimension(npassive) :: anglr, g0_r, g0_z
    double precision, dimension(258, 258) :: C_00, S_00

!first, initialized initial guess coming from prescribed boundary current density: jrhotheta
    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, Rrect(1), dr)
    jaxis = closest_index(zax, Zrect(1), dz)
    iax = iaxis
    jax = jaxis

! Evaluate coils' quantities
    do i=1, npassive
        anglr(i) = ATAN2(z_cond(nactive + i) - zax, r_cond(nactive + i) - rax)
! Find true axis
        xub(1) = Rrect(iax-1)
        xub(2) = Rrect(iax)
        xub(3) = Rrect(iax + 1)
        xub(4) = Rrect(iax)
        xub(5) = Rrect(iax)
        xub(6) = Rrect(iax-1)
        xub(7) = Rrect(iax-1)
        xub(8) = Rrect(iax + 1)
        xub(9) = Rrect(iax + 1)
        yub(1) = Zrect(jax)
        yub(2) = Zrect(jax)
        yub(3) = Zrect(jax)
        yub(4) = Zrect(jax-1)
        yub(5) = Zrect(jax + 1)
        yub(6) = Zrect(jax-1)
        yub(7) = Zrect(jax + 1)
        yub(8) = Zrect(jax-1)
        yub(9) = Zrect(jax + 1)
        bub(1) = greeni(iax - 1, jax    , nactive + i)
        bub(2) = greeni(iax    , jax    , nactive + i)
        bub(3) = greeni(iax + 1, jax    , nactive + i)
        bub(4) = greeni(iax    , jax - 1, nactive + i)
        bub(5) = greeni(iax    , jax + 1, nactive + i)
        bub(6) = greeni(iax - 1, jax - 1, nactive + i)
        bub(7) = greeni(iax - 1, jax + 1, nactive + i)
        bub(8) = greeni(iax + 1, jax - 1, nactive + i)
        bub(9) = greeni(iax + 1, jax + 1, nactive + i)
        ddipsi = least_square_biquad(xub, yub, bub, 9)
        g0_r(i) = ddipsi(1)
        g0_z(i) = ddipsi(2)
    enddo

    do j=1, nz2
        do i=1, nr2
            C_00(i, j) = sum(greeni(i, j, nactive + 1:nactive + npassive)*cos(anglr))
            S_00(i, j) = sum(greeni(i, j, nactive + 1:nactive + npassive)*sin(anglr))
        enddo
    enddo
    C_00r = sum(g0_r*cos(anglr))
    S_00r = sum(g0_r*sin(anglr))
    C_00z = sum(g0_z*cos(anglr))
    S_00z = sum(g0_z*sin(anglr))

    psistab1o = 1000.
    psistab2o = 1000.

    do j_iter=1, 30000
        psiplasrz = get_psiplasrz()
        psistabr = 0.
        psistabz = 0.
        delr = 0.
        delz = 0.
        call compound_psi
        call find_new_axis
        dum1 = C_00r*S_00z - C_00z*S_00r

        bub(1) = interp2d_psi(raxp + dr, zaxp, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxp - dr, zaxp, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxp, zaxp + dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxp, zaxp - dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))

        xub(1) = (bub(1) - bub(2))/(2.*dr)
        yub(1) = (bub(3) - bub(4))/(2.*dz)
        delr = ( S_00z*(-xub(1)) - C_00z*(-yub(1)))/dum1
        delz = (-S_00r*(-xub(1)) + C_00r*(-yub(1)))/dum1
        psistabr = delr
        psistabz = delz

        call compound_psi
        psirz(1:nr2, 1:nz2) = psirz(1:nr2, 1:nz2) +  &
            psistabr*C_00(1:nr2, 1:nz2) + psistabz*S_00(1:nr2, 1:nz2) !total flux

        call find_new_axis
        call find_psi_boundary
        call new_jrz  ! calculate new right hand side

        temp_err = (abs(psistab1o - psistabr) + abs(psistab2o - psistabz))
        psistab1o = psistabr
        psistab2o = psistabz

        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, npassive
        curconduc(nactive + i) = curconduc(nactive + i) + psistabr*cos(anglr(i)) + psistabz*sin(anglr(i))
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis
    call find_psi_boundary
    call new_jrz  ! calculate new right hand side

    end subroutine restab_axis_with_fourier_wall

!--------------------------------------------------------------------
    subroutine estimate_tau_VDE(n_coil_de, n_coil_st, i_coil_de, i_coil_st, tau_LR, tau_VDE, F_stab, F_destab)

    use green_function, only: greeni
    use pi_const, only: GP2, GP2_sq

    integer, intent(IN) :: n_coil_de, n_coil_st, i_coil_de(n_coil_de), i_coil_st(n_coil_st)
    double precision, intent(IN) :: tau_LR
    
    double precision, intent(OUT) :: tau_VDE
    
    integer i,i1,j,k
    double precision :: F_stab, F_destab, f_ratio, dum1, dum2
    
! Calculate the time scale of the VDE based on the simple estimate tau_VDE = ()    

! Compute F_destab
    F_destab = 0.
    do i1=1, n_coil_de
        i = i_coil_de(i1)  
        dum1 = 0.
        do j=2, nz1
            do k=2, nr1
                dum2 = 1./dz**2 * ( greeni(k, j+1, i) - 2.*greeni(k, j, i) + greeni(k, j-1, i) )
                dum1 = dum1 + jrz(k,j) * curconduc(i) * dum2
            enddo
        enddo
        F_destab = F_destab + dum1
    enddo

! Compute F_stab
    F_stab = 0.
    do i1=1, n_coil_st
        i = i_coil_st(i1)  
        dum1 = 0.
        do j=2,nz1
            do k=2,nr1
                dum2 = GP2/indconduc(i, i)*( 1./(2.*dz) * (greeni(k, j+1, i) - greeni(k, j-1, i)) )**2
                dum1 = dum1 + jrz(k,j) * dum2
            enddo
        enddo
        F_stab = F_stab + dum1
    enddo
 
    F_stab = F_stab * iplasma
    f_ratio = F_stab / F_destab     
    tau_VDE = tau_LR * (f_ratio - 1.)

    end subroutine estimate_tau_VDE

!-------------------------------------------------------------------
    subroutine ferro_mag_create

    use ferromagstructure, only: type_ferromag
    use feqis_tools, only: interp2d_psi, green_function, inv_matrix
    use numerical_tools, only: qinterp
    use pi_const, only: GP, muvac

    integer :: i, j, ii, jj, iii, iferro, nval
    double precision :: x1, x2, x3, x4, d
    double precision, dimension(1) :: z1, z2
    type(type_ferromag), dimension(:), allocatable :: ferromag

! Loop over ferromagnetic elements
    do iferro=1, nferromag
        iii  = ferromag(iferro)%position%npoints
        nval = ferromag(iferro)%mhrelation%nvalues
! Construct vacuum field
        do j=1, iii
            x1 = interp2d_psi(ferromag(iferro)%position%r(j) + dr/2, ferromag(iferro)%position%z(j), Rrect, Zrect, psirz)
            x2 = interp2d_psi(ferromag(iferro)%position%r(j) - dr/2, ferromag(iferro)%position%z(j), Rrect, Zrect, psirz)
            x3 = interp2d_psi(ferromag(iferro)%position%r(j), ferromag(iferro)%position%z(j) + dz/2, Rrect, Zrect, psirz)
            x4 = interp2d_psi(ferromag(iferro)%position%r(j), ferromag(iferro)%position%z(j) - dz/2, Rrect, Zrect, psirz)

            z1(1) = -2./(x1 + x2)*(x1 - x2)/dr ! Bz=-1/r dpsi/dr
            z2(1) =  2./(x1 + x2)*(x3 - x4)/dz ! Br=1/r dpsi/dz
            ferromag(iferro)%position%btangfield(j) = z1(1)*sin(ferromag(iferro)%position%tanangl(j)) + &
                                                      z2(1)*cos(ferromag(iferro)%position%tanangl(j)) ! btangent
! Find magnetization chi
            if (ferromag(iferro)%position%btangfield(j) <= ferromag(iferro)%mhrelation%h(1)) then
                ferromag(i)%position%magnetizationchi(j) = ferromag(iferro)%mhrelation%chi(1)
            endif
            if (ferromag(iferro)%position%btangfield(j) >= ferromag(iferro)%mhrelation%h(nval)) then
                ferromag(iferro)%position%magnetizationchi(j)=ferromag(iferro)%mhrelation%chi(nval)
            endif
            if (ferromag(iferro)%position%btangfield(j) < ferromag(iferro)%mhrelation%h(nval) .and. &
                ferromag(iferro)%position%btangfield(j) > ferromag(iferro)%mhrelation%h(1)) then
                z1(1) = ferromag(iferro)%position%btangfield(j)
                call qinterp(ferromag(iferro)%mhrelation%h, ferromag(iferro)%mhrelation%chi, nval, z1, z2, 1)
                ferromag(iferro)%position%magnetizationchi(j) = z2(1)
            endif
        enddo

! Create matrix
        do jj=1, iii
            do ii=1, iii
                matrix_ferro_to_invert(ii, jj) = 1. + ferromag(iferro)%mutual_matrix%Mij(ii, jj)
            enddo
        enddo

! Calculate inverse and currents
        matrix_ferro_inverse(1:iii, 1:iii) = inv_matrix(matrix_ferro_to_invert(1:iii, 1:iii), iii)
        do ii=1, iii
            ferromag(iferro)%position%current(ii) = sum(matrix_ferro_inverse(ii, 1:iii)*ferromag(iferro)%position%btangfield(1:iii)) !ferromag currents in MA
        enddo
    enddo

! Calculate psi from ferromagnetics
    do jj=1, nz2
        do ii=1, nr2
            do iferro=1, nferromag
                iii = ferromag(iferro)%position%npoints
                do j=1, iii
                    d = abs(Rrect(ii) - ferromag(iferro)%position%r(j)) + abs(Zrect(ii) - ferromag(iferro)%position%z(j))
                   if (d > 0) then
                       psiferro(ii, jj) = psiferro(ii, jj) + muvac/GP*green_function(Rrect(ii), Zrect(jj), ferromag(iferro)%position%r(j), ferromag(iferro)%position%z(j))*ferromag(iferro)%position%current(j)
                   endif
               enddo
           enddo
       enddo
    enddo

    end subroutine ferro_mag_create

!---------------------------------------------------------------------
    subroutine interp_j_fromrhotorz

    integer :: i, j

! go from jrhotheta to jrz
    jrz = 0.
    jrhotheta(1:nrho, ntheta+1) = jrhotheta(1:nrho, 1)
    do j=1, nz2
        do i=1, nr2
            jrz(i, j) = curinterp(Rrect(i), Zrect(j), jrhotheta(1:nrho, 1:ntheta+1),  &
                rho(1:nrho, 1:ntheta+1), theta(1:ntheta+1), raxp, zaxp, nrho, ntheta+1)
        enddo
    enddo

    end subroutine interp_j_fromrhotorz

end module circuit
