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

contains
!---------------------------------------------------------------------
    subroutine get_zccurb_feqis(rc_cur, zc_cur, z2c_cur, rgeoc, zgeoc, ahorc)

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
        do j=1, nteta
            avgelem = dator(i, j) !/bpcell(i, j-1)
            perimz = perimz + avgelem
            rgeoc  = rgeoc + rpol(i, j)*avgelem
            zgeoc  = zgeoc + zpol(i, j)*avgelem
        enddo
    enddo
    rgeoc = rgeoc/perimz
    zgeoc = zgeoc/perimz

    do i=1, nrho-1
        do j=1, nteta
            avgelem = dator(i, j) !/bpcell(i, j-1)
            ahorc2  = ahorc2 + (rpol(i, j) - rgeoc)**2 * avgelem
        enddo
    enddo
    ahorc2 = ahorc2/perimz
    ahorc = 2.*sqrt(ahorc2)

    return
    end subroutine get_zccurb_feqis

!---------------------------------------------------------------------
    double precision function psib_ext_feqis
! Returns external flux on plasma boundary

    use feqis_tools, only: interp2d_psi
    
    integer :: i
    double precision :: psiext_out, dlt, dllt, psi_ext_1, psiext1, psiext2

!cycle over boundary
    psiext_out = 0.
    dllt = 0.
    psi_ext_1 = interp2d_psi(rbnd(1), zbnd(1), r(1:nr), z(1:nz), psiextrz(1:nr, 1:nz))
    psiext1 = psi_ext_1
    do i=1, nbnd-1
        psiext2 = interp2d_psi(rbnd(i+1), zbnd(i+1), r(1:nr), z(1:nz), psiextrz(1:nr, 1:nz))
        dlt = sqrt((rbnd(i+1) - rbnd(i))**2 + (zbnd(i+1) - zbnd(i))**2)
        psiext_out = psiext_out + 0.5*(psiext1 + psiext2)*dlt
        dllt = dllt + dlt
        psiext1 = psiext2
    enddo
    psiext1 = psi_ext_1
    dlt = sqrt((rbnd(1) - rbnd(nbnd))**2 + (zbnd(1) - zbnd(nbnd))**2)
    psiext_out = psiext_out + 0.5*(psiext1 + psiext2)*dlt
    dllt = dllt + dlt

    psib_ext_feqis = psiext_out/dllt

    write(*, *) 'psibbb', psibnd, psiext_out

    return
    end function psib_ext_feqis

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
        u001 = interp2d_psi(x1, y1, r(1:nr), z(1:nz), psirz(1: nr, 1:nz))
        u002 = interp2d_psi(x2, y2, r(1:nr), z(1:nz), psirz(1: nr, 1:nz))
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

    return
    end function find_l_gap

!---------------------------------------------------------------------
    function find_demo_gaps_feqis(ngaps, demo_gaps) result(geom1d)

    integer, intent(in) :: ngaps
    double precision, intent(in) :: demo_gaps(ngaps, 4)
    double precision :: geom1d(ngaps)

    integer :: i, onlypos
    double precision :: d_step, gapmin, gapmax, l_gap, l_gap_pos, l_gap_neg

    d_step = 0.1 ! advance in 1 cm steps
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

    return
    end function find_demo_gaps_feqis

!---------------------------------------------------------------------
    subroutine psi_external_calc

    use green_matrix, only: greeni

    integer :: i, j

    do j=1, nz2
        do i=1, nr2
            psiextrz(i, j) = sum(curconduc(1: nconduc)*greeni(i, j, 1: nconduc))
        enddo
    enddo

    return
    end subroutine psi_external_calc

!---------------------------------------------------------------------
    subroutine nine_point_regression(r_in, z_in, pos_xpoint, ddpsi, f00)

    use feqis_tools, only: closest_index, exact_biquad

    integer, parameter :: ndim=9
    double precision, intent(in) :: r_in, z_in
    double precision, intent(out) :: f00
    double precision, intent(out), dimension(2) :: pos_xpoint 
    double precision, intent(out), dimension(ndim-1) :: ddpsi

    integer :: iax, jax, i, j, k
    double precision :: rax, zax
    double precision, dimension(ndim) :: bub

    iax = closest_index(r_in, r(1), dr)
    jax = closest_index(z_in, z(1), dz)

    k = 0
    do j=-1, 1
        do i=-1, 1
            k = k + 1
            bub(k) = psirz(iax+i, jax+j)
        enddo
    enddo
    rax = r(iax)
    zax = z(jax)

    call exact_biquad(rax, zax, bub(1:ndim), ndim,  &
        pos_xpoint(1), pos_xpoint(2), f00, ddpsi, dr, dz)

    return
    end subroutine nine_point_regression

!---------------------------------------------------------------------
    subroutine nine_point_coeffs_only(r_in, z_in, c, rax_out, zax_out)

    use feqis_tools, only: closest_index, A_inv

    integer, parameter :: ndim=9
    double precision, intent(in) :: r_in, z_in
    double precision, intent(out) :: rax_out, zax_out
    double precision, intent(out), dimension(9) :: c

    integer :: iax, jax, i, j, k
    double precision, dimension(90) :: bub

    iax = closest_index(r_in, r(1), dr)
    jax = closest_index(z_in, z(1), dz)

    rax_out = r(iax)
    zax_out = z(jax)

!find true axis
    k = 0
    do j=-1, 1
        do i=-1, 1
            k = k + 1
            bub(k) = psirz(iax+i, jax+j)
        enddo
    enddo

    do k=1, 9
        c(k) = sum(A_inv(k, 1:ndim) * bub(1:ndim))
    enddo

    return
    end subroutine nine_point_coeffs_only

!---------------------------------------------------------------------
    subroutine nine_point_regression_follow(rx, zx, pos_xpoint, ddpsi, f00)

    use feqis_tools, only: interp2d_psi, exact_biquad_regress

    integer, parameter :: ndim=9
    double precision, intent(in) :: rx, zx
    double precision, intent(out) :: f00
    double precision, intent(out), dimension(2) :: pos_xpoint
    double precision, intent(out), dimension(ndim-1) :: ddpsi
    
    double precision, dimension(ndim) :: bub

    bub(1) = interp2d_psi(rx - dr, zx - dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
    bub(2) = interp2d_psi(rx     , zx - dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
    bub(3) = interp2d_psi(rx + dr, zx - dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
    bub(4) = interp2d_psi(rx - dr, zx     , r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
    bub(5) = interp2d_psi(rx     , zx     , r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
    bub(6) = interp2d_psi(rx + dr, zx     , r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
    bub(7) = interp2d_psi(rx - dr, zx + dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
    bub(8) = interp2d_psi(rx     , zx + dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
    bub(9) = interp2d_psi(rx + dr, zx + dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))

    call exact_biquad_regress(rx, zx, bub(1:ndim), ndim,  &
        pos_xpoint(1), pos_xpoint(2), f00, ddpsi, dr, dz)

    return
    end subroutine nine_point_regression_follow

!---------------------------------------------------------------------
    subroutine find_closest_xpoints(rx, zx, ierr, n_add)

!this routine finds the x-points close to the plasma boundary,  irrespective of other x-points

    integer, intent(out) :: ierr, n_add
    double precision, intent(out), dimension(20) :: rx, zx

    integer :: i, jinc, nx
    integer, dimension(250) :: jcycl
    double precision :: bx0, bx1, x1
    double precision, dimension(2) :: posx
    double precision, dimension(8) :: ddipsi

    bx1 = sqrt(dr**2 + dz**2)
    rx = 1000.
    zx = 1000.
    nx    = 0
    ierr  = 1
    n_add = 0
    jinc  = 0

    if (rbnd(1) <= r(1)) return ! boundary doesnt exist yet

    ierr = 0

    do i=1, nteta !cycle over boundary points
!around each boundary point,  do a 3-layer X-point search (25 point search x boundary point)
        call nine_point_regression_follow(rbnd(i), zbnd(i), posx, ddipsi, x1)
        if (.not. isnan(x1)) then
            if (posx(1) <= lim_maxR .and. posx(1) >= lim_minR .and. posx(2) >= lim_minZ .and. posx(2) <= lim_maxZ) then
                jinc = jinc + 1
                jcycl(jinc) = i
            endif
        endif
    enddo

    if (jinc == 0) then
        ierr = 1
        return
    endif

    call nine_point_regression_follow(rbnd(jcycl(1)), zbnd(jcycl(1)), posx, ddipsi, x1)
    n_add = 1
    rx(1) = posx(1)
    zx(1) = posx(2)

    if (jinc == 1) then
        return
    endif

    if (jinc >= 2) then
!remove double counts
        do i=2, jinc
            call nine_point_regression_follow(rbnd(jcycl(i)), zbnd(jcycl(i)), posx, ddipsi, x1)
            bx0 = sqrt((rx(i-1) - posx(1))**2 + (zx(i-1)-posx(2))**2)
            if (bx0 > bx1) then
                n_add = n_add + 1
                rx(n_add) = posx(1)
                zx(n_add) = posx(2)
            endif
        enddo
    endif

    if (n_add == 0) ierr=1

    return
    end subroutine find_closest_xpoints

!-----------------------------------------------------------------------------------
    function boundary(green_in) result(green_out)
! New bc is integral_over_boundary of -Green * dg/dn * dl

    use pi_vars, only: GPI

    double precision, intent(in), dimension(i_dim2, i_dim2) :: green_in
    double precision, dimension(i_dim2, i_dim2) :: green_out

    integer :: i, jcounty
    double precision, dimension(i_dim2, 4) :: integr

    integr = 0.
    jcounty = 0

! lower side
    do i=2, nr1
        integr(i, 1) = bgint(green_in, jcounty)
    enddo

! right side
    do i=2, nz1
        integr(i, 2) = bgint(green_in, jcounty)
    enddo

! upper side
    do i=2, nr1
        integr(i, 3) = bgint(green_in, jcounty)
    enddo

! left side
    do i=2, nz1
        integr(i, 4) = bgint(green_in, jcounty)
    enddo

    green_out(2:nr1,   1) = integr(2:nr1, 1)/GPI
    green_out(nr2, 2:nz1) = integr(2:nz1, 2)/GPI
    green_out(2:nr1, nz2) = integr(2:nr1, 3)/GPI
    green_out(  1, 2:nz1) = integr(2:nz1, 4)/GPI

    return
    end function boundary

!-----------------------------------------------------------------------------------
    double precision function bgint(green_in, jcounty)
! Integral_over_boundary of -Green * dg/dn * dl

    double precision, intent(in), dimension(i_dim2, i_dim2) :: green_in
    integer, intent(inout) :: jcounty

    integer :: j
    double precision :: dgdn(i_dim2), greenf

    bgint = 0.

! lower side
    do j=2, nr1
        jcounty = jcounty + 1
        greenf  = green_bnd_f(jcounty)
        dgdn(j) = -green_in(j, 2)/dz*greenf*dr/r(j)
    enddo
    bgint = bgint - sum(dgdn(2:nr1))

!right side
    do j=2, nz1
        jcounty = jcounty + 1
        greenf  = green_bnd_f(jcounty)
        dgdn(j) = -green_in(nr1, j)/dr*greenf*dz/(r(nr2) + r(nr1))*2.
    enddo
    bgint = bgint - sum(dgdn(2:nz1))

! upper side
    do j=2, nr1
        jcounty = jcounty + 1
        greenf  = green_bnd_f(jcounty)
        dgdn(j) = -green_in(j, nz1)/dz*greenf*dr/r(j)
    enddo
    bgint = bgint - sum(dgdn(2:nr1))

!left side
    do j=2, nz1
        jcounty = jcounty + 1
        greenf  = green_bnd_f(jcounty)
        dgdn(j) = -green_in(2, j)/dr*greenf*dz/(r(1) + r(2))*2.
    enddo
    bgint = bgint - sum(dgdn(2:nz1))

    return
    end function bgint

!---------------------------------------------------------------------
    integer function xpoint_axis_connection(rx, zx, rax, zax, dr, dz)

!this routine checks that going from axis to x point,  the directed gradient of psi never changes sign
!(otherwise it means the x point is not connected to the plasma

    use feqis_tools, only: interp2d_psi, pol_angle

    double precision,  intent(in) :: rx, zx, rax, zax, dr, dz

    integer :: nsteps, i
    double precision :: angl, dbl, dd, t1, t2, t3, t4, t5
    double precision :: z1, z2, psiold, z3

    angl = pol_angle(rax, zax, rx, zx)

    dbl = sqrt((rx - rax)**2 + (zx - zax)**2)
    dd  = sqrt(dr**2 + dz**2)
    nsteps = nint(dbl/dd)
    dd = dbl/nsteps !perfect ratio
  
    xpoint_axis_connection = 1
    psiold = 0.
    do i=2, nsteps
        t1 = rax + dd*(i - 1)*cos(angl)
        t2 = zax + dd*(i - 1)*sin(angl)
        t3 = rax + dd*i*cos(angl)
        t4 = zax + dd*i*sin(angl)
        z1 = interp2d_psi(t1, t2, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        z2 = interp2d_psi(t3, t4, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        z3 = (z2 - z1)*psiold
        psiold = z2 - z1
        if (z3 < 0) then
            xpoint_axis_connection = 0
            EXIT
        endif
        if (dd*i >= dbl) then
            EXIT
        endif
    enddo

    return
    end function xpoint_axis_connection


end module feqis_circuit
