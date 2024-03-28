module feqis_circuit

use feqis_dimensions, only: i_dim1, i_dim2, i_dim5, nrho2d

implicit none

integer :: use_limiter_yesno

!time stepping
double precision :: tau_old, tau_new
double precision, dimension(:), allocatable :: psi_cur_old, dpc

!circuits
integer :: ncoils, nreseqcoil, nlimiter
integer, dimension(:), allocatable :: mequivalence
double precision, dimension(:), allocatable :: limiterR, limiterZ
double precision :: lim_maxR, lim_minR, lim_maxZ, lim_minZ
double precision, dimension(:), allocatable :: Rcoil, Zcoil, drcoil, dzcoil, &
    anglecoil, anglehcoil

integer :: nconduc, nblocks, npassive, nactive, nferromag
double precision, dimension(:), allocatable :: curconduc, voltage, voltage_old, &
    cur_con_old, r_cond, z_cond
double precision, dimension(:, :), allocatable :: resconduc, indconduc
double precision, dimension(:), allocatable :: psiplasmatoconduc !plasma --> conduc at t

! coordinates:
! r, z --> rectangular grid in meters
! rho --> distance of point from magnetic axis in meters
! teta --> angle from LFS midplane which is 0, coincides with Zmag plane
! psigrid --> in poloidal flux equispaced

!grids
integer :: nr, nz, nrho, nteta, nr2, nz2, nr1, nz1, nbnd, ngbnd, &
    redo_bnd
integer, dimension(:, :), allocatable :: zlimpotential

double precision :: rmin, rmax, zmin, zmax, dr, dz, dteta, &
    zbot, ztop, raus, rinner
double precision, dimension(:), allocatable :: r, z, rcomp, zcomp
double precision, dimension(:), allocatable :: teta, psigrid
double precision, dimension(:, :), allocatable :: rho, area_eff, &
    rpol, zpol, rpul, zpul, psirhoteta

double precision, dimension(:, :), allocatable :: u_n, omega_pl, &
    psirz, psiextrz, psiplasrz, psiferro

! r(z)pol: R, Z in polar coordinates half radial grid
! r(z)pul: R, Z in polar coordinates full radial grid
double precision, dimension(nrho2d) :: psia_2d, ffp_2d, ppp_2d
double precision, dimension(:), allocatable :: psia_1d, ffp_1d, ppp_1d
double precision, dimension(8) :: derivpsi

! boundary and axis FBE, PBE
integer, parameter :: max_xpoints=500
integer :: i_plasmatype !(0-limited, 1-single null, 2-double null)
integer :: iaxis, jaxis, n_of_xpoints
double precision :: psibnd, psiaxis, psibndp, psiaxisp, &
    rax, zax, trax, tzax, raxp, zaxp, &
    alpsep, psistabR, psistabZ, dr_factor_init, dz_factor_init
double precision, dimension(:), allocatable :: rbnd, zbnd
double precision, dimension(max_xpoints) :: r_xpoint, z_xpoint
double precision :: deriv_x(5, max_xpoints)
double precision, dimension(:), allocatable :: green_bnd_f
double precision, dimension(:), allocatable :: rbndp, zbndp, rexp, zexp, tetaexp

! plasma parameters
double precision :: iplasma, btor0, rgeom0, psplex, li3, li_aug, betapol
double precision, dimension(:), allocatable :: pprime, ffprime, pressure, psigrida, ipol
double precision, dimension(:, :), allocatable :: jrz
double precision, dimension(:, :), allocatable :: jrhoteta

!ferromag
double precision, dimension(100,100) :: matrix_ferro_to_invert, matrix_ferro_inverse

contains

!---------------------------------------------------------------------
    subroutine psi_mutual_effect_conductors(psi_to_conductors)

    use green_matrix, only: greeni

    implicit none

    double precision, intent(out) :: psi_to_conductors(nconduc)
    integer :: i

    do i=1, nconduc
        psi_to_conductors(i) = sum(jrz(1: nr2, 1: nz2) * area_eff(1: nr2, 1: nz2) * greeni(1: nr2, 1: nz2, i))
    enddo

    return
    end subroutine psi_mutual_effect_conductors

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

!    write(*, *) 'psibbb', psibnd, psiext_out

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

    call exact_biquad(rax, zax, bub(1:ndim), ndim, &
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
    double precision, intent(out), dimension(200) :: rx, zx

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

    double precision, intent(in), dimension(nr2, nz2) :: green_in
    double precision, dimension(nr2, nz2) :: green_out

    integer :: i, jcounty
    double precision, dimension(nr2+nz2, 4) :: integr

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

    double precision, intent(in), dimension(nr2, nz2) :: green_in
    integer, intent(inout) :: jcounty

    integer :: j
    double precision :: dgdn(1000), greenf

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

!--------------------------------------------------------------------
    subroutine restab_F_function_full_fonfit
! Refits all currents

    use errors_params, only: err_find_psistab
    use astra2fbe, only: sigma_coils, sigma_b, sigma_axis
    use green_matrix, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix

    integer :: i, j, k, j_iter, iax, jax
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, raxref, zaxref
    double precision, dimension(9) :: bub
    double precision, dimension(500) :: rbref, zbref
    double precision, dimension(nconduc) :: G_00c, G_00r, G_00z, &
        Fderiv, curref, curnow, curdiff
    double precision, dimension(nteta) :: psicorr
    double precision, dimension(nr2, nz2) :: g
    double precision, dimension(nconduc, nteta) :: G_00
    double precision, dimension(nconduc, nconduc) :: matrix, invmatrix

    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhoteta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, r(1), dr)
    jaxis = closest_index(zax, z(1), dz)
    iax = iaxis
    jax = jaxis

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    sigma_coils(nactive + 1:nconduc) = sigma_coils(nactive + 1)

    curref(1:nconduc) = curconduc(1:nconduc)
    curnow(1:nconduc) = curconduc(1:nconduc)
    curdiff = 0.

    raxref = raxp
    zaxref = zaxp
    rbref(1:nteta) = rbndp(1:nteta)
    zbref(1:nteta) = zbndp(1:nteta)

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.
    do j=1, nconduc
        do k=1, nteta
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:nteta))/(0. + nteta)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
    enddo

    do j=1, nconduc
        do i=1, nconduc
            if (i == j) matrix(i, j) = matrix(i, j) + 2.*sigma_coils(i)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:nteta) - G_00c(i))*(G_00(j, 1:nteta) - G_00c(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j))
       enddo
    enddo

! Calculate inverse
    invmatrix = inv_matrix(matrix, nconduc)

    do j_iter=1, 300000 !iterations to find currents
        if (j_iter > 150) stop
        g = 0.
        call solve_gs2d(g) !jrz as right hand side
        g = boundary(g)  ! gbound = integral (Green*dg/dn) over the boundary

        call solve_gs2d(g) ! again jrz as right hand side
        psiplasrz(1:nr2, 1:nz2) = g(1:nr2, 1:nz2)

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

        do j=1, nteta
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(nteta + 0.) !average psi on the boundary

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz

        Ffunc = sigma_B*sum((psicorr - x1)**2) + sum(sigma_coils(1:nconduc)*curdiff**2) + sigma_axis*(x2**2 + x3**2)

! Calculate F derivative
        do i=1, nconduc
            Fderiv(i) = 2.*(sigma_coils(i)*curdiff(i) + &
                sigma_B*sum((psicorr - x1)*(G_00(i, 1:nteta) - G_00c(i))) + &
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

        call find_new_axis_part1
        call find_psi_boundary
        call new_jrz_feqis  ! calculate new right hand side

        temp_err = abs(Ffunc - Ffunc_old)
        Ffunc_old = Ffunc
        write(*,*) 'iteration ',j_iter, curnow(1:nconduc),rax,zax,temp_err, Ffunc,sum(abs(Fderiv(1:nconduc)))
        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nconduc
        curconduc(i) = curnow(i)
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis_part1
    call find_psi_boundary
    call new_jrz_feqis  ! calculate new right hand side

    write(*, *) curconduc(1:nconduc), rax, zax
    write(*, *) 'full fonfit eddy currents converged',temp_err,err_find_psistab

    return
    end subroutine restab_F_function_full_fonfit

!--------------------------------------------------------------------

!--------------------------------------------------------------------
    subroutine restab_F_function_full_fonfit_xpoints ! valid only if n_xpoint_fit > 0
! Refits all currents and x points

    use errors_params, only: err_find_psistab
    use astra2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_xpoint, &
        r_xpoint_fit, z_xpoint_fit, n_xpoint_fit
    use green_matrix, only: greeni
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

    double precision, dimension(nteta) :: psicorr
    double precision, dimension(nr2, nz2) :: g
    double precision, dimension(nconduc, nteta) :: G_00
    double precision, dimension(nconduc, nconduc) :: matrix, invmatrix

    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhoteta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, r(1), dr)
    jaxis = closest_index(zax, z(1), dz)
    iax = iaxis
    jax = jaxis

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    sigma_coils(nactive + 1:nconduc) = sigma_coils(nactive + 1)

    curref(1:nconduc) = curconduc(1:nconduc)
    curnow(1:nconduc) = curconduc(1:nconduc)
    curdiff = 0.

    raxref = raxp
    zaxref = zaxp
    rbref(1:nteta) = rbndp(1:nteta)
    zbref(1:nteta) = zbndp(1:nteta)

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.
    do j=1, nconduc
        do k=1, nteta
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:nteta))/(0. + nteta)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
        do k=1, n_xpoint_fit
            bub(1) = interp2d_psi(r_xpoint_fit(k) - dr/2., z_xpoint_fit(k), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            bub(2) = interp2d_psi(r_xpoint_fit(k) + dr/2., z_xpoint_fit(k), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(k) - dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(k) + dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            G_00xr(j,k) = (bub(2) - bub(1))/dr
            G_00xz(j,k) = (bub(4) - bub(3))/dz
        enddo
    enddo

    do j=1, nconduc
        do i=1, nconduc
            if (i == j) matrix(i, j) = matrix(i, j) + 2.*sigma_coils(i)
            dummyx(1:n_xpoint_fit) = G_00xr(i,1:n_xpoint_fit)*G_00xr(j,1:n_xpoint_fit) + G_00xz(i,1:n_xpoint_fit)*G_00xz(j,1:n_xpoint_fit)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:nteta) - G_00c(i))*(G_00(j, 1:nteta) - G_00c(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j)) + &
                2.*sigma_xpoint*sum(dummyx)
       enddo
    enddo

! Calculate inverse
    invmatrix = inv_matrix(matrix, nconduc)

    do j_iter=1, 300000 !iterations to find currents
        if (j_iter > 150) stop
        g = 0.
        call solve_gs2d(g) !jrz as right hand side
        g = boundary(g)  ! gbound = integral (Green*dg/dn) over the boundary

        call solve_gs2d(g) ! again jrz as right hand side
        psiplasrz(1:nr2, 1:nz2) = g(1:nr2, 1:nz2)

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

        do j=1, nteta
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(nteta + 0.) !average psi on the boundary

        Ffunc = 0.
        Ffunc = sigma_B*sum((psicorr - x1)**2) + sum(sigma_coils(1:nconduc)*curdiff**2)

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz
        Ffunc = Ffunc +  sigma_axis*(x2**2 + x3**2)
        do k=1, n_xpoint_fit
            bub(1) = interp2d_psi(r_xpoint_fit(k) - 0.5*dr, z_xpoint_fit(K), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(2) = interp2d_psi(r_xpoint_fit(k) + 0.5*dr, z_xpoint_fit(K), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            x2 = (bub(2) - bub(1))/dr ! dPsi/dr
            x3 = (bub(4) - bub(3))/dz ! dPsi/dz
            Ffunc = Ffunc +  sigma_xpoint*(x2**2 + x3**2)
        enddo



! Calculate F derivative
        do i=1, nconduc
            Fderiv(i) = 2.*(sigma_coils(i)*curdiff(i) + sigma_B*sum((psicorr - x1)*(G_00(i, 1:nteta) - G_00c(i))))
            bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            x2 = (bub(2) - bub(1))/dr ! dPsi/dr
            x3 = (bub(4) - bub(3))/dz ! dPsi/dz
            Fderiv(i) = Fderiv(i) + 2.*sigma_axis*(x2*G_00r(i) + x3*G_00z(i))
            do k=1, n_xpoint_fit
                bub(1) = interp2d_psi(r_xpoint_fit(k) - 0.5*dr, z_xpoint_fit(K), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
                bub(2) = interp2d_psi(r_xpoint_fit(k) + 0.5*dr, z_xpoint_fit(K), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
                bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
                bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
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

        call find_new_axis_part1
        call find_psi_boundary
        call new_jrz_feqis  ! calculate new right hand side

        temp_err = abs(Ffunc - Ffunc_old)
        Ffunc_old = Ffunc
        write(*,*) 'iteration ',j_iter, curnow(1:nconduc), rax, zax, temp_err, Ffunc, sum(abs(Fderiv(1:nconduc)))
        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nconduc
        curconduc(i) = curnow(i)
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis_part1
    call find_psi_boundary
    call new_jrz_feqis  ! calculate new right hand side

    write(*, *) curconduc(1:nconduc), rax, zax
    write(*, *) 'full fonfit eddy currents converged', temp_err, err_find_psistab

    return
    end subroutine restab_F_function_full_fonfit_xpoints

!--------------------------------------------------------------------
    subroutine restab_2_timepoints_evolution

! Finds active currents from scratch including evolution from time point t1 to time point t2

    use errors_params, only: err_find_psistab
    use astra2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy
    use green_matrix, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix
    use pi_vars, only: GPI, GPI2, GPI4, muvac, mu0
    use numerical_tools, only: linterp

    integer, parameter :: n_evol = 2
    integer :: i, j, k, j_iter, iax, jax, jt
    integer, dimension(n_evol) :: iax_ev, jax_ev, nxp_ev
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, lambda, deltapsiext !lagrange multplier lambda
    double precision, dimension(9) :: bub
    double precision, dimension(n_evol) :: rax_ev,zax_ev,ip_ev, V_loop, L_ext, delta_t, psi_ext_ev
    double precision, dimension(nteta) :: psicorr
    double precision, dimension(:), allocatable :: raxref_ev, zaxref_ev, &
        ffunc_ev, result_vector, Fderiv
    double precision, dimension(nrho, n_evol) :: psia_ev,pprim_ev,ffprim_ev
    double precision, dimension(nr2, nz2) :: g
    double precision, dimension(nactive, n_evol) :: i_totev
    double precision, dimension(1000, n_evol) :: rxp_ev, zxp_ev
    double precision, dimension(:, :), allocatable :: rbref_ev, zbref_ev, G_00c, G_00r, G_00z, &
        curref, curnow, curdiff, matrix, invmatrix
    double precision, dimension(:, :, :), allocatable :: j_ev, psi_ev, G_00
    character(len=80) :: file_time

    allocate(rbref_ev(nteta, n_evol))
    allocate(zbref_ev(nteta, n_evol))
    allocate(raxref_ev(n_evol))
    allocate(zaxref_ev(n_evol))
    allocate(ffunc_ev(n_evol))
    allocate(j_ev(nr2, nz2, n_evol))
    allocate(psi_ev(nr2, nz2, n_evol))
    allocate(G_00(nactive, nteta, n_evol))
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
        if (jt == 1) file_time = 'dat/feqis_time1_file.dat' !contains the prescribed boundary data for time t1: jrhoteta, rho, teta, rb, zb, rho augmented to nteta+1, teta augmented to nteta+1, raxp, zaxp, iplasma, psia,pprim,ffprim
        if (jt == 2) file_time = 'dat/feqis_time2_file.dat' !contains the prescribed boundary data for time t2: jrhoteta, rho, teta, rho augmented to nteta+1, teta augmented to nteta+1, raxp, zaxp, iplasma, psia2,pprim2,ffprim2, Lext, deltaT, Vloop
        open(32, file=TRIM(file_time))
            read(32, *) jrhoteta(1:nrho, 1:nteta)
            read(32, *) rho(1:nrho, 1:nteta)
            read(32, *) teta(1:nteta)
            read(32, *) rbndp(1:nteta)
            read(32, *) zbndp(1:nteta)
            read(32, *) raxp, zaxp, ip_ev(jt)
            read(32, *) psia_ev(1:nrho, jt), pprim_ev(1:nrho, jt), ffprim_ev(1:nrho, jt) !pprime is Pascal / grad(FP), ffprime is F dF/dFP
            read(32, *) L_ext(jt), delta_t(jt), V_loop(jt)
        close(32)
        rho(1:nrho, nteta+1) = rho(1:nrho, 1)
        teta(nteta+1) = teta(1) + GPI2
        call interp_j_fromrhotorz
        curr = SUM(jrz)*dr*dz
        j_ev(:, :, jt) = jrz/curr*ip_ev(jt)
        psia_ev(:, jt) = (psia_ev(:, jt) - psia_ev(1, jt))/(psia_ev(nrho, jt) - psia_ev(1, jt))
        rax_ev(jt) = raxp
        iax_ev(jt) = closest_index(rax_ev(jt), rmin, dr)
        zax_ev(jt) = zaxp
        jax_ev(jt) = closest_index(zax_ev(jt), zmin, dz)
        rbref_ev(1: nteta, jt) = rbndp(1: nteta)
        zbref_ev(1: nteta, jt) = zbndp(1: nteta)
    enddo

    deltapsiext = -(0.5*sum(V_loop)*(delta_t(2) - delta_t(1)) + 0.5*sum(L_ext)*(ip_ev(2) - ip_ev(1)))

    write(*, *) 'deltapsi', deltapsiext, rax_ev, zax_ev, j_ev(30, 30, :)

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
            do k=1, nteta
                G_00(j, k, jt) = interp2d_psi(rbref_ev(k,jt), zbref_ev(k,jt), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            enddo
            G_00c(j, jt) = sum(G_00(j, 1:nteta, jt))/(0. + nteta)
            bub(1) = interp2d_psi(raxref_ev(jt) - dr/2., zaxref_ev(jt), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            bub(2) = interp2d_psi(raxref_ev(jt) + dr/2., zaxref_ev(jt), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            bub(3) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) - dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            bub(4) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) + dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            G_00r(j,jt) = (bub(2) - bub(1))/dr
            G_00z(j,jt) = (bub(4) - bub(3))/dz
        enddo

        do j=1, nactive
            do i=1, nactive
                if (i == j) matrix((jt-1)*nactive+i, (jt-1)*nactive+j) = matrix((jt-1)*nactive+i, (jt-1)*nactive+j) + sigma_coils(i)*sigma_energy*indconduc(i, i)
                matrix((jt-1)*nactive+i, (jt-1)*nactive+j) = matrix((jt-1)*nactive+i, (jt-1)*nactive+j) +  &
                    2.*sigma_B*sum((G_00(i, 1:nteta, jt) - G_00c(i, jt))*(G_00(j, 1:nteta, jt) - G_00c(j, jt))) +  &
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
    write(*, *) 'invmatrix', invmatrix(10, 10)

    do j_iter=1, 300000 !iterations to find currents
        if (j_iter > 150) stop
        do jt=1, n_evol
            jrz = j_ev(:, :, jt)   ! assign previous current density
            g = 0.
            call solve_gs2d(g) ! jrz as right hand side
            g = boundary(g)    ! gbound = integral (Green*dg/dn) over the boundary
            call solve_gs2d(g) ! again jrz as right hand side
            psiplasrz = g ! solution for pure plasma
            psi_ev(:, :, jt) = g
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
           do j=1, nteta
               psicorr(j) = interp2d_psi(rbref_ev(j,jt), zbref_ev(j,jt), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
           enddo
           x1 = sum(psicorr)/(nteta + 0.) !average psi on the boundary
! Derivative at ref axis
           bub(1) = interp2d_psi(raxref_ev(jt) - 0.5*dr, zaxref_ev(jt), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
           bub(2) = interp2d_psi(raxref_ev(jt) + 0.5*dr, zaxref_ev(jt), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
           bub(3) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
           bub(4) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
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
                   sigma_B*sum((psicorr - x1)*(G_00(i, 1:nteta, jt) - G_00c(i, jt))) + &
                   sigma_axis*(x2*G_00r(i, jt) + x3*G_00z(i, jt)) )
           enddo
       enddo !end time loop

       write(*, *) 'psi ext', psi_ext_ev, lambda

! add lambda contributions
       Ffunc = sum(Ffunc_ev) + lambda*(deltapsiext - (psi_ext_ev(2) - psi_ext_ev(1)))
       do i=1, nactive
           Fderiv(i) = Fderiv(i) + lambda*G_00c(i, 1)
       enddo
       do i=1, nactive
           Fderiv(nactive+i) = Fderiv(nactive+i) - lambda*G_00c(i, 2)
       enddo
       Fderiv(2*nactive+1)=deltapsiext-(psi_ext_ev(2)-psi_ext_ev(1)) ! dF/dlambda

       write(*, *) 'fderiv', fderiv

! Calculate new currents
       do i=1, n_evol*nactive+1
           result_vector(i) = result_vector(i) - sum(invmatrix(i, 1:n_evol*nactive+1)*Fderiv)
       enddo
       curdiff(:, 1) = result_vector(1:nactive)
       curdiff(:, 2) = result_vector(nactive+1:n_evol*nactive)
       lambda = result_vector(n_evol*nactive+1)

       write(*, *) 'result_vector', result_vector
       write(*, *) 'stop jere in fonfit times'

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
           call find_new_axis_part1
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
           ffp_2d = -GPI2/mu0*ffp_2d
           ppp_2d = -GPI2*1.e-6*ppp_2d
           call new_jrz_feqis  ! calculate new right hand side
           j_ev(:, :, jt) = jrz
           i_totev(1: nactive, jt) = curdiff(1: nactive, jt)
       enddo

       temp_err = sum(abs(Fderiv)) !error
       Ffunc_old = Ffunc
       write(*, *) 'iteration ', j_iter, temp_err, sum(abs(Fderiv(1: 2*nactive))), rax, zax
       if (temp_err <= err_find_psistab) EXIT
    enddo

    write(*, *) ''
    write(*, *) 'full fonfit active currents converged',temp_err,err_find_psistab
    write(*, *) ''
    write(*, *) 'closing program but saving data in output_timefit.dat'
    open(32, file='dat/output_timefit.dat')
        write(32,*) 't1 currents: ', i_totev(1: nactive, 1)
        write(32,*) 't2 currents: ', i_totev(1: nactive, 2)
    close(32)
    pause

    return
    end subroutine restab_2_timepoints_evolution

!--------------------------------------------------------------------
    subroutine restab_2_timepoints_evolution_limits

! Finds active currents from scratch including evolution from time point t1 to time point t2

    use errors_params, only: err_find_psistab
    use astra2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy, &
        current_limit_feqis
    use green_matrix, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix
    use pi_vars, only: GPI, GPI2, GPI4, muvac, mu0
    use numerical_tools, only: linterp

    integer, parameter :: n_evol = 2
    integer :: i, j, k, j_iter, iax, jax, jt, iactive, jactive
    integer, dimension(n_evol) :: iax_ev, jax_ev, nxp_ev
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, lambda, deltapsiext !lagrange multplier lambda
    double precision, dimension(9) :: bub
    double precision, dimension(nteta) :: psicorr
    double precision, dimension(n_evol) :: rax_ev, zax_ev, ip_ev, &
         V_loop, L_ext, delta_t, psi_ext_ev
    double precision, dimension(:), allocatable :: raxref_ev, zaxref_ev, &
         ffunc_ev, result_vector, Fderiv
    double precision, dimension(nrho, n_evol) :: psia_ev,pprim_ev,ffprim_ev
    double precision, dimension(1000, n_evol) :: rxp_ev, zxp_ev
    double precision, dimension(nr2, nz2) :: g
    double precision, dimension(nactive, n_evol) :: i_totev
    double precision, dimension(:, :), allocatable :: rbref_ev, zbref_ev, &
        G_00c, G_00r, G_00z, curref, curnow, curdiff, matrix, invmatrix
    double precision, dimension(:, :, :), allocatable :: j_ev, psi_ev, G_00
    character(len=80) :: file_time

    allocate(rbref_ev(nteta, n_evol))
    allocate(zbref_ev(nteta, n_evol))
    allocate(raxref_ev(n_evol))
    allocate(zaxref_ev(n_evol))
    allocate(ffunc_ev(n_evol))
    allocate(j_ev(nr2, nz2, n_evol))
    allocate(psi_ev(nr2, nz2, n_evol))
    allocate(G_00(nactive, nteta, n_evol))
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
        if (jt==1) file_time='dat/feqis_time1_file.dat'   !contains the prescribed boundary data for time t1: jrhoteta, rho, teta, rb, zb, rho augmented to nteta+1, teta augmented to nteta+1, raxp, zaxp, iplasma, psia,pprim,ffprim
        if (jt==2) file_time='dat/feqis_time2_file.dat'   !contains the prescribed boundary data for time t2: jrhoteta, rho, teta, rho augmented to nteta+1, teta augmented to nteta+1, raxp, zaxp, iplasma, psia2,pprim2,ffprim2, Lext, deltaT, Vloop
        open(32, file=TRIM(file_time))
            read(32, *) jrhoteta(1: nrho, 1: nteta)
            read(32, *) rho(1: nrho, 1: nteta)
            read(32, *) teta(1: nteta)
            read(32, *) rbndp(1: nteta)
            read(32, *) zbndp(1: nteta)
            read(32, *) raxp, zaxp, ip_ev(jt)
            read(32, *) psia_ev(1: nrho,jt), pprim_ev(1: nrho, jt), ffprim_ev(1: nrho, jt) !pprime is Pascal / grad(FP), ffprime is F dF/dFP
            read(32,*) L_ext(jt), delta_t(jt), V_loop(jt)
        close(32)
        rho(1: nrho, nteta+1) = rho(1: nrho, 1)
        teta(nteta+1) = teta(1) + GPI2
        call interp_j_fromrhotorz
        curr = SUM(jrz) *dr*dz
        j_ev(:, :, jt) = jrz/curr*ip_ev(jt)
        psia_ev(:, jt) = (psia_ev(:, jt) - psia_ev(1, jt))/(psia_ev(nrho, jt) - psia_ev(1, jt))
        rax_ev(jt) = raxp
        iax_ev(jt) = closest_index(rax_ev(jt), rmin, dr)
        zax_ev(jt) = zaxp
        jax_ev(jt) = closest_index(zax_ev(jt), zmin, dz)
        rbref_ev(1: nteta, jt) = rbndp(1: nteta)
        zbref_ev(1: nteta, jt) = zbndp(1: nteta)
    enddo

    deltapsiext = -(0.5*sum(V_loop)*(delta_t(2) - delta_t(1)) + 0.5*sum(L_ext) * (ip_ev(2) - ip_ev(1)))
    write(*, *) 'deltapsi', deltapsiext, rax_ev, zax_ev, j_ev(30, 30, :)

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
            do k=1, nteta
                G_00(j, k, jt) = interp2d_psi(rbref_ev(k, jt), zbref_ev(k, jt), r(1: nr), z(1: nz), greeni(1: nr, 1: nz, j))
            enddo
            G_00c(j, jt) = sum(G_00(j, 1:nteta, jt))/(0. + nteta)
            bub(1) = interp2d_psi(raxref_ev(jt) - dr/2., zaxref_ev(jt), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            bub(2) = interp2d_psi(raxref_ev(jt) + dr/2., zaxref_ev(jt), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            bub(3) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) - dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
            bub(4) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) + dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
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
                   2.*sigma_B*sum((G_00(i, 1:nteta, jt) - G_00c(i, jt))*(G_00(j, 1:nteta, jt) - G_00c(j, jt))) +  &
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
    write(*, *) 'invmatrix', invmatrix(10, 10)

    do j_iter=1, 300000 ! iterations to find currents
        if (j_iter > 150) stop
        do jt=1, n_evol
            jrz = j_ev(:, :, jt)    ! assign previous current density
            g = 0.
            call solve_gs2d(g) ! jrz as right hand side
            g = boundary(g)    ! gbound = integral (Green*dg/dn) over the boundary
            call solve_gs2d(g) ! again jrz as right hand side
            psiplasrz = g      ! solution for pure plasma
            psi_ev(:, :, jt) = g
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
            do j=1, nteta
                psicorr(j) = interp2d_psi(rbref_ev(j, jt), zbref_ev(j, jt), r(1: nr), z(1: nz), psirz(1: nr, 1: nz))
            enddo
            x1 = sum(psicorr)/(nteta + 0.) !average psi on the boundary
! Derivative at ref axis
            bub(1) = interp2d_psi(raxref_ev(jt) - 0.5*dr, zaxref_ev(jt), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(2) = interp2d_psi(raxref_ev(jt) + 0.5*dr, zaxref_ev(jt), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(3) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(4) = interp2d_psi(raxref_ev(jt), zaxref_ev(jt) + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
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
                    sigma_B*sum((psicorr - x1)*(G_00(i, 1:nteta, jt) - G_00c(i, jt))) + &
                    sigma_axis*(x2*G_00r(i, jt) + x3*G_00z(i, jt)) )
            enddo
        enddo !end time loop

        write(*, *) 'psi ext', psi_ext_ev, lambda

! Add lambda contributions
        Ffunc = sum(Ffunc_ev) + lambda*(deltapsiext - (psi_ext_ev(2) - psi_ext_ev(1)))
        do i=1, nactive
            Fderiv(i) = Fderiv(i) + lambda*G_00c(i, 1)
        enddo
        do i=1, nactive
            Fderiv(nactive+i) = Fderiv(nactive+i) - lambda*G_00c(i, 2)
        enddo
        Fderiv(2*nactive+1) = deltapsiext - (psi_ext_ev(2) - psi_ext_ev(1)) ! dF/dlambda
        write(*, *) 'fderiv', fderiv

! Calculate new currents
        do i=1, n_evol*nactive+1
            result_vector(i) = result_vector(i) - sum(invmatrix(i, 1: n_evol*nactive+1)*Fderiv)
        enddo
        curdiff(:, 1) = result_vector(1: nactive)
        curdiff(:, 2) = result_vector(nactive+1: n_evol*nactive)
        lambda = result_vector(n_evol*nactive+1)
        write(*, *) 'result_vector', result_vector
        write(*, *) 'stop jere in fonfit times'
! Apply limits
        do jt=1, n_evol
            do i=1, nactive
                curdiff(i, jt) = min(curdiff(i, jt), current_limit_feqis(i, 1))
                curdiff(i, jt) = max(curdiff(i, jt), current_limit_feqis(i, 2))
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
            call find_new_axis_part1
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
            ffp_2d = -GPI2/mu0*ffp_2d
            ppp_2d = -GPI2*1.e-6*ppp_2d
            call new_jrz_feqis  ! calculate new right hand side
            j_ev(:, :, jt) = jrz

            i_totev(1: nactive, jt) = curdiff(1: nactive, jt)
        enddo

        temp_err = sum(abs(Fderiv)) !error
        Ffunc_old = Ffunc
        write(*, *) 'iteration ', j_iter, temp_err, sum(abs(Fderiv(1: 2*nactive))), rax, zax
        if (temp_err <= err_find_psistab) EXIT
    enddo

    write(*, *) ''
    write(*, *) 'full fonfit active currents converged',temp_err,err_find_psistab
    write(*, *) ''
    write(*, *) 'closing program but saving data in output_timefit.dat'
    open(32, file='dat/output_timefit.dat')
        write(32, *) 't1 currents: ', i_totev(1: nactive, 1)
        write(32, *) 't2 currents: ', i_totev(1: nactive, 2)
    close(32)

    pause

    return
    end subroutine restab_2_timepoints_evolution_limits

!--------------------------------------------------------------------
    subroutine restab_j_timepoints_evolution_limits_xpoints

! Finds active currents from scratch including evolution from time point j-1 to point j

! to calculate forces using this mode, call: call coil_forces_feqis(ncoil_blocks, force_R, force_Z, 0)   !the 0 at the end means no plasma contribution. also eddy currents are ignored.

    use errors_params, only: err_find_psistab
    use astra2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy, &
        current_limit_feqis, sigma_xpoint, r_xpoint_fit, z_xpoint_fit, &
       n_xpoint_fit, sigma_coils_psiext, vloop_avg, L_ext, &
       dIp_dt, tau_gseq_feqis, time_astra
    use green_matrix, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix
    use pi_vars, only: GPI, GPI2, GPI4, muvac, mu0
    use numerical_tools, only: linterp

    integer :: i, j, k, j_iter, iax, jax, jt, j_time
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, lambda, deltapsiext, psibexto, psibext !lagrange multplier lambda
    double precision, dimension(9) :: bub
    double precision, dimension(nteta) :: psicorr, rbref, zbref
    double precision, dimension(nactive+1) :: result_vector, Fderiv
    double precision, dimension(nactive, nteta) :: G_00
    double precision, dimension(nr2, nz2) :: g
    double precision, dimension(nactive) :: G_00c, G_00r, G_00z, curref, curnow, curdiff
    double precision, dimension(nactive+1,nactive+1) :: matrix, invmatrix
    double precision, dimension(:, :), allocatable :: G_00xr, G_00xz
    double precision, dimension(:), allocatable :: dummyx
    double precision, dimension(nblocks-npassive) :: force_r, force_z
    double precision :: raxref, zaxref
    character(len=80) :: file_time
    logical :: file_exists
    data j_time/0/
    save j_time

    deltapsiext = tau_gseq_feqis*(vloop_avg+L_ext*dIp_dt)   ! jump in psiext
    write(*, *) 'deltapsi', deltapsiext, tau_gseq_feqis, vloop_avg, L_ext, dIp_dt, &
   curconduc(1:nactive), j_time

    if (n_xpoint_fit > 0) then
        allocate(G_00xr(nactive, n_xpoint_fit))
        allocate(G_00xz(nactive, n_xpoint_fit))
        allocate(dummyx(n_xpoint_fit))
    endif

    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhoteta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, r(1), dr)
    jaxis = closest_index(zax, z(1), dz)
    iax = iaxis
    jax = jaxis
    raxref = raxp
    zaxref = zaxp
    rbref(1:nteta) = rbndp(1:nteta)
    zbref(1:nteta) = zbndp(1:nteta)
    nbnd = nteta
    rbnd(1:nbnd) = rbref(1:nteta)
    zbnd(1:nbnd) = zbref(1:nteta)

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)
    curref(1:nactive) = curconduc(1:nactive)
    curnow(1:nactive) = curconduc(1:nactive)
    curdiff = 0.

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.
    lambda = 0.

    write(*, *) 'psiext', psiextrz(30,30), rax, zax, jrz(30, 30)

    do j=1, nactive
        do k=1, nteta
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:nteta))/(0. + nteta)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
        if (n_xpoint_fit>0) then
            do k=1, n_xpoint_fit
                bub(1) = interp2d_psi(r_xpoint_fit(k) - dr/2., z_xpoint_fit(k), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
                bub(2) = interp2d_psi(r_xpoint_fit(k) + dr/2., z_xpoint_fit(k), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
                bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(k) - dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
                bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(k) + dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
                G_00xr(j,k) = (bub(2) - bub(1))/dr
                G_00xz(j,k) = (bub(4) - bub(3))/dz
            enddo
        endif
    enddo

    psibexto = sum(G_00c*curnow)

    do j=1, nactive
        do i=1, nactive
            if (i == j) matrix(i, j) = matrix(i, j) + sigma_energy*sigma_coils(i)*indconduc(i,i) + 2.*sigma_coils(i)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:nteta) - G_00c(i))*(G_00(j, 1:nteta) - G_00c(j))) +  &
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
        matrix(nactive+1, i) = G_00c(i)
        matrix(i, nactive+1) = G_00c(i)
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
        g = 0.
        call solve_gs2d(g) !jrz as right hand side
        g = boundary(g)  ! gbound = integral (Green*dg/dn) over the boundary

        call solve_gs2d(g) ! again jrz as right hand side
        psiplasrz(1:nr2, 1:nz2) = g(1:nr2, 1:nz2)

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

        do j=1, nteta
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(nteta + 0.) !average psi on the boundary

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz

        Ffunc = sigma_B*sum((psicorr - x1)**2) + sigma_axis*(x2**2 + x3**2)
        do i=1, nactive
            Ffunc = Ffunc+0.5*indconduc(i, i)*sigma_energy*sigma_coils(i)*curnow(i)**2 + sigma_coils(i)*curdiff(i)**2
        enddo

        if (n_xpoint_fit > 0) then
            do k=1, n_xpoint_fit
                bub(1) = interp2d_psi(r_xpoint_fit(k) - 0.5*dr, z_xpoint_fit(K), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
                bub(2) = interp2d_psi(r_xpoint_fit(k) + 0.5*dr, z_xpoint_fit(K), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
                bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
                bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
                x2 = (bub(2) - bub(1))/dr ! dPsi/dr
                x3 = (bub(4) - bub(3))/dz ! dPsi/dz
                Ffunc = Ffunc +  sigma_xpoint*(x2**2 + x3**2)
             enddo
        endif

        psibext = sum(G_00c*curnow)

        write(*,*) 'loop', deltapsiext/GPI2, psibext, psibexto, psibext - psibexto, -deltapsiext/GPI2
        if (j_time > 0) Ffunc = Ffunc + lambda*(psibext - psibexto + deltapsiext/GPI2)

! Calculate F derivative
        Fderiv = 0.
        do i=1, nactive
            Fderiv(i) = 2.*(0.5*indconduc(i,i)*sigma_energy*sigma_coils(i)*curnow(i) + &
                sigma_coils(i)*curdiff(i) + sigma_B*sum((psicorr - x1)*(G_00(i, 1:nteta) - G_00c(i))))
            bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            x2 = (bub(2) - bub(1))/dr ! dPsi/dr
            x3 = (bub(4) - bub(3))/dz ! dPsi/dz
            Fderiv(i) = Fderiv(i) + 2.*sigma_axis*(x2*G_00r(i) + x3*G_00z(i))
            if (n_xpoint_fit > 0) then
                do k=1, n_xpoint_fit
                    bub(1) = interp2d_psi(r_xpoint_fit(k) - 0.5*dr, z_xpoint_fit(K), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
                    bub(2) = interp2d_psi(r_xpoint_fit(k) + 0.5*dr, z_xpoint_fit(K), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
                    bub(3) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
                    bub(4) = interp2d_psi(r_xpoint_fit(k), z_xpoint_fit(K) + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
                    x2 = (bub(2) - bub(1))/dr ! dPsi/dr
                    x3 = (bub(4) - bub(3))/dz ! dPsi/dz
                    Fderiv(i) = Fderiv(i) +  2.*sigma_xpoint*(x2*G_00xr(i,k) + x3*G_00xz(i,k))
                enddo
            endif
        enddo

        if (j_time > 0) Fderiv(1:nactive) = Fderiv(1:nactive) + lambda*G_00c(1:nactive)
        if (j_time > 0) Fderiv(nactive+1) = (psibext-psibexto+deltapsiext/GPI2)

        write(*, *) 'fderiv', fderiv

! Calculate new currents

        do i=1, nactive
           if (j_time == 0) then
               curnow(i) = curnow(i) - sum(invmatrix(i, 1:nactive)*Fderiv(1:nactive))
           else
               curnow(i) = curnow(i) - sum(invmatrix(i, 1:nactive+1)*Fderiv(1:nactive+1))
           endif
        enddo
        lambda = lambda - sum(invmatrix(nactive+1, 1:nactive+1)*Fderiv(1:nactive+1))

! cut new currents to limits
        do i=1, nactive
            curnow(i) = min(curnow(i), current_limit_feqis(i, 1))
            curnow(i) = max(curnow(i), current_limit_feqis(i, 2))
        enddo
        curdiff = curnow - curref

        write(*, *) 'result vector', curnow(1:nactive), lambda

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

        call find_new_axis_part1
        call find_psi_boundary
        call new_jrz_feqis  ! calculate new right hand side

        temp_err = sum(abs(Fderiv))
        Ffunc_old = Ffunc
        write(*,*) 'iteration ',j_iter, curnow(1:nactive),rax,zax,temp_err,sum(abs(Fderiv(1:nactive)))

        if (j_iter > 15000) stop
        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nactive
        curconduc(i) = curnow(i)
    enddo

    write(*, *) ' '
    write(*, *) curconduc(1:nconduc)
    write(*, *) 'full fonfit active currents converged', temp_err, err_find_psistab

    write(file_time,'(1A, I0)') 'dat/file_fit_currents_output_', j_time

    inquire(FILE=file_time, EXIST=file_exists)
    if (file_exists) call execute_command_line('rm '//trim(file_time), wait=.true.)

! Evaluate forces
    call coil_forces_feqis(nblocks - npassive, force_R, force_Z, 0)

    call psi_mutual_effect_conductors(psiplasmatoconduc)

    write(*, *) 'Saving new currents in ', file_time
    open(32, file = file_time)
    write(32, *) nactive, curconduc(1:nactive)*1e3, -GPI2*psibext, rax, zax, &
        tau_gseq_feqis, vloop_avg, L_ext*dIp_dt, &
        nr2, nz2, psirz(1:nr2, 1:nz2), nblocks - npassive, force_R, force_Z, &
        psibnd, psiaxis, indconduc(1:nactive, 1:nactive), resconduc(1:nactive, 1:nactive), &
        nlimiter, limiterR(1:nlimiter), limiterZ(1:nlimiter), r(1), r(nr2), z(1), z(nz2), nteta, &
        rbref(1:nteta), zbref(1:nteta), time_astra, tau_gseq_feqis, GPI2*psiplasmatoconduc(1:nactive)  !saved in kA
    close(32)

!update time
    j_time = j_time +1

    if (n_xpoint_fit > 0) then
        deallocate(G_00xr)
        deallocate(G_00xz)
        deallocate(dummyx)
    endif

    return
    end subroutine restab_j_timepoints_evolution_limits_xpoints
  
!--------------------------------------------------------------------
    subroutine diagnose_feqis(filename)  !call it in sbr/assign_geom.f90 in astra for example, it has access to this module.

    use metric_coefficients_pbe, only: dator

    implicit none

    integer, parameter :: unit=32
    character(*), intent(in) :: filename

    open(unit, file=trim(filename))
!first PBE stuff
        write(unit, *) nrho, nteta
        write(unit, *) pprime, ffprime, psigrida  ! pprim, ffprim, psigrid
        write(unit, *) rpol, zpol, jrhoteta, psirhoteta  !R, Z, jrhoteta, PSI
        write(unit, *) dator !area elements
!now free boundary
        if (allocated(r)) then
            write(unit, *) nr2, nz2, r, z !grid
            write(unit, *) psiextrz, psiplasrz, psirz !psivacuum, psiplasma, psitotal maps [radiants]
            write(unit, *) jrz, area_eff !current densiy, area elements
            write(unit, *) nconduc, r_cond, z_cond !conductors positions
            write(unit, *) curconduc, voltage, psiplasmatoconduc !currents, voltages, psiplasma mutual induct
            write(unit, *) nlimiter, limiterr, limiterz !limiter
        else
            write(unit, *) '-1'
        endif
    close(unit)

    return
    end subroutine diagnose_feqis
  
!--------------------------------------------------------------------
    subroutine restab_F_function_full_currents
! Finds active currents from scratch. Passive currents are given.

    use errors_params, only: err_find_psistab
    use astra2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy
    use green_matrix, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix

    integer :: i, j, k, j_iter, iax, jax
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, raxref, zaxref
    double precision, dimension(9) :: bub
    double precision, dimension(500) :: rbref, zbref
    double precision, dimension(nactive) :: G_00c, G_00r, G_00z, &
        Fderiv, curref, curnow, curdiff
    double precision, dimension(nteta) :: psicorr
    double precision, dimension(nr2, nz2) :: g
    double precision, dimension(nactive, nteta) :: G_00
    double precision, dimension(nactive, nactive) :: matrix, invmatrix

    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhoteta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, r(1), dr)
    jaxis = closest_index(zax, z(1), dz)
    iax = iaxis
    jax = jaxis

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    curref(1:nactive) = 0.
    curnow(1:nactive) = 0.
    curconduc(1:nactive) = 0.
    curdiff = 0.

    raxref = raxp
    zaxref = zaxp
    rbref(1:nteta) = rbndp(1:nteta)
    zbref(1:nteta) = zbndp(1:nteta)

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.

    write(*,*) 'psiext',psiextrz(30,30),rax,zax,jrz(30,30)

    do j=1, nactive
        do k=1, nteta
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:nteta))/(0. + nteta)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
    enddo

    do j=1, nactive
        do i=1, nactive
            if (i == j) matrix(i, j) = matrix(i, j) + sigma_energy*sigma_coils(i)*indconduc(i,i)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:nteta) - G_00c(i))*(G_00(j, 1:nteta) - G_00c(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j))
        enddo
    enddo

! Calculate inverse
    invmatrix = inv_matrix(matrix, nactive)

! Iteration to find currents
    do j_iter=1, 300000 
        if (j_iter > 150) stop
        g = 0.
        call solve_gs2d(g) ! jrz as right hand side
        g = boundary(g)    ! gbound = integral (Green*dg/dn) over the boundary

        call solve_gs2d(g) ! again jrz as right hand side
        psiplasrz(1:nr2, 1:nz2) = g(1:nr2, 1:nz2)

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

        do j=1, nteta
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(nteta + 0.) !average psi on the boundary

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz

        Ffunc = sigma_B*sum((psicorr - x1)**2) + sigma_axis*(x2**2 + x3**2)
        do i=1,nactive
            Ffunc = Ffunc + 0.5*indconduc(i, i)*sigma_energy*sigma_coils(i)*curdiff(i)**2
        enddo

! Calculate F derivative
        do i=1, nactive
            Fderiv(i) = 2.*(0.5*indconduc(i,i)*sigma_energy*sigma_coils(i)*curdiff(i) + &
                sigma_B*sum((psicorr - x1)*(G_00(i, 1:nteta) - G_00c(i))) + &
                sigma_axis*(x2*G_00r(i) + x3*G_00z(i)) )
        enddo

        write(*, *) 'fderiv', fderiv

! Calculate new currents

        do i=1, nactive
            curnow(i) = curnow(i) - sum(invmatrix(i, 1:nactive)*Fderiv)
        enddo
        curdiff = curnow - curref
        write(*,*) 'result vector',curdiff(1:nactive)

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

        call find_new_axis_part1
        call find_psi_boundary
        call new_jrz_feqis  ! calculate new right hand side

        temp_err = sum(abs(Fderiv))
        Ffunc_old = Ffunc
        write(*,*) 'iteration ', j_iter, curnow(1:nactive), rax, zax, temp_err, sum(abs(Fderiv(1:nactive)))

        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nactive
        curconduc(i) = curnow(i)
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis_part1
    call find_psi_boundary
    call new_jrz_feqis  ! calculate new right hand side

    write(*, *) ' '
    write(*, *) curconduc(1:nconduc), rax, zax
    write(*, *) 'full fonfit active currents converged', temp_err, err_find_psistab

    return
    end subroutine restab_F_function_full_currents

!--------------------------------------------------------------------
    subroutine restab_F_function_full_currents_limits
! Finds active currents from scratch. Passive currents are given.

    use errors_params, only: err_find_psistab
    use astra2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy, &
      current_limit_feqis
    use green_matrix, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix

    integer :: i, j, k, j_iter, iax, jax
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, raxref, zaxref
    double precision, dimension(9) :: bub
    double precision, dimension(500) :: rbref, zbref
    double precision, dimension(nactive) :: G_00c, G_00r, G_00z, &
        Fderiv, curref, curnow, curdiff
    double precision, dimension(nteta) :: psicorr
    double precision, dimension(nr2, nz2) :: g
    double precision, dimension(nactive, nteta) :: G_00
    double precision, dimension(nactive, nactive) :: matrix, invmatrix

    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhoteta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, r(1), dr)
    jaxis = closest_index(zax, z(1), dz)
    iax = iaxis
    jax = jaxis

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    curref(1:nactive) = 0.
    curnow(1:nactive) = 0.
    curconduc(1:nactive) = 0.
    curdiff = 0.

    raxref = raxp
    zaxref = zaxp
    rbref(1:nteta) = rbndp(1:nteta)
    zbref(1:nteta) = zbndp(1:nteta)

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.

    write(*,*) 'psiext',psiextrz(30,30),rax,zax,jrz(30,30)

    do j=1, nactive
        do k=1, nteta
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:nteta))/(0. + nteta)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
    enddo

    do j=1, nactive
        do i=1, nactive
            if (i == j) matrix(i, j) = matrix(i, j) + sigma_energy*sigma_coils(i)*indconduc(i,i)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:nteta) - G_00c(i))*(G_00(j, 1:nteta) - G_00c(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j))
       enddo
    enddo

! Calculate inverse
    invmatrix = inv_matrix(matrix, nactive)
! Iteration to find currents
    do j_iter=1, 300000
        if (j_iter > 150) stop
        g = 0.
        call solve_gs2d(g) ! jrz as right hand side
        g = boundary(g)    ! gbound = integral (Green*dg/dn) over the boundary
        call solve_gs2d(g) ! again jrz as right hand side
        psiplasrz(1:nr2, 1:nz2) = g(1:nr2, 1:nz2)

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

        do j=1, nteta
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(nteta + 0.) !average psi on the boundary

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz

        Ffunc=sigma_B*sum((psicorr - x1)**2) + sigma_axis*(x2**2 + x3**2)
        do i=1, nactive
            Ffunc = Ffunc + 0.5*indconduc(i,i)*sigma_energy*sigma_coils(i)*curdiff(i)**2
        enddo

! Calculate F derivative
        do i=1, nactive
            Fderiv(i) = 2.*(0.5*indconduc(i,i)*sigma_energy*sigma_coils(i)*curdiff(i) + &
                sigma_B*sum((psicorr - x1)*(G_00(i, 1:nteta) - G_00c(i))) + &
                sigma_axis*(x2*G_00r(i) + x3*G_00z(i)) )
        enddo

        write(*, *) 'fderiv', fderiv

! Calculate new currents

        do i=1, nactive
            curnow(i) = curnow(i) - sum(invmatrix(i, 1:nactive)*Fderiv)
        enddo
        curdiff = curnow - curref
        write(*, *) 'result vector', curdiff(1:nactive)

! cut new currents to limits
        do i=1, nactive
            curdiff(i) = min(curdiff(i), current_limit_feqis(i, 1))
            curdiff(i) = max(curdiff(i), current_limit_feqis(i, 2))
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

        call find_new_axis_part1
        call find_psi_boundary
        call new_jrz_feqis  ! calculate new right hand side

        temp_err = sum(abs(Fderiv))
        Ffunc_old = Ffunc
        write(*,*) 'iteration ', j_iter, curnow(1:nactive), rax, zax, temp_err, sum(abs(Fderiv(1:nactive)))

        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nactive
        curconduc(i) = curnow(i)
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis_part1
    call find_psi_boundary
    call new_jrz_feqis  ! calculate new right hand side

    write(*, *) ' '
    write(*, *) curconduc(1:nconduc), rax, zax
    write(*, *) 'full fonfit active currents converged',temp_err, err_find_psistab

    return
    end subroutine restab_F_function_full_currents_limits

!--------------------------------------------------------------------
    subroutine restab_F_function_full_currents_forces
! Finds active currents from scratch. Passive currents are given. Forces get minimized too.

! how to include forces???

    use errors_params, only: err_find_psistab
    use astra2fbe, only: sigma_coils, sigma_b, sigma_axis, sigma_energy, sigma_forces
    use green_matrix, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, inv_matrix

    integer :: i, j, k, j_iter, iax, jax
    double precision :: temp_err, curr, f_correction, x1, x2, x3, &
        Ffunc, Ffunc_old, raxref, zaxref
    double precision, dimension(9) :: bub
    double precision, dimension(500) :: rbref, zbref
    double precision, dimension(nactive) :: G_00c, G_00r, G_00z, &
        Fderiv, curref, curnow, curdiff
    double precision, dimension(nteta) :: psicorr
    double precision, dimension(nr2, nz2) :: g
    double precision, dimension(nactive, nteta) :: G_00
    double precision, dimension(nactive, nactive) :: matrix, invmatrix

    psicorr   = 0.
    Ffunc_old = 1.e6

! First, initialized initial guess coming from prescribed boundary current density: jrhoteta

    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, r(1), dr)
    jaxis = closest_index(zax, z(1), dz)
    iax = iaxis
    jax = jaxis

! Axis is given by raxp, zaxp; boundary by rbndp, zbndp; coilref by curconduc(passive)

    curref(1:nactive) = 0.
    curnow(1:nactive) = 0.
    curconduc(1:nactive) = 0.
    curdiff = 0.

    raxref = raxp
    zaxref = zaxp
    rbref(1:nteta) = rbndp(1:nteta)
    zbref(1:nteta) = zbndp(1:nteta)

! calculate the matrix F_li of the F function, including the green function terms
    matrix    = 0.
    invmatrix = 0.

    write(*,*) 'psiext', psiextrz(30, 30), rax, zax, jrz(30, 30)

    do j=1, nactive
        do k=1, nteta
            G_00(j, k) = interp2d_psi(rbref(k), zbref(k), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        enddo
        G_00c(j) = sum(G_00(j, 1:nteta))/(0. + nteta)
        bub(1) = interp2d_psi(raxref - dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(2) = interp2d_psi(raxref + dr/2., zaxref, r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(3) = interp2d_psi(raxref, zaxref - dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        bub(4) = interp2d_psi(raxref, zaxref + dz/2., r(1:nr), z(1:nz), greeni(1:nr, 1:nz, j))
        G_00r(j) = (bub(2) - bub(1))/dr
        G_00z(j) = (bub(4) - bub(3))/dz
    enddo

    do j=1, nactive
        do i=1, nactive
            if (i == j) matrix(i, j) = matrix(i, j) + sigma_energy*sigma_coils(i)*indconduc(i,i)
            matrix(i, j) = matrix(i, j) +  &
                2.*sigma_B*sum((G_00(i, 1:nteta) - G_00c(i))*(G_00(j, 1:nteta) - G_00c(j))) +  &
                2.*sigma_axis*(G_00r(i)*G_00r(j) + G_00z(i)*G_00z(j))
       enddo
    enddo

! Calculate inverse
    invmatrix = inv_matrix(matrix, nactive)
! Iteration to find currents
    do j_iter=1, 300000
        if (j_iter > 150) stop
        g = 0.
        call solve_gs2d(g) !jrz as right hand side
        g = boundary(g)  ! gbound = integral (Green*dg/dn) over the boundary

        call solve_gs2d(g) ! again jrz as right hand side
        psiplasrz(1:nr2, 1:nz2) = g(1:nr2, 1:nz2)

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

        do j=1, nteta
            psicorr(j) = interp2d_psi(rbref(j), zbref(j), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(nteta + 0.) ! average psi on the boundary

! Derivative at ref axis
        bub(1) = interp2d_psi(raxref - 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxref + 0.5*dr, zaxref, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxref, zaxref - 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxref, zaxref + 0.5*dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        x2 = (bub(2) - bub(1))/dr ! dPsi/dr
        x3 = (bub(4) - bub(3))/dz ! dPsi/dz

        Ffunc = sigma_B*sum((psicorr - x1)**2) + sigma_axis*(x2**2 + x3**2)
        do i=1, nactive
            Ffunc = Ffunc + 0.5*indconduc(i, i)*sigma_energy*sigma_coils(i)*curdiff(i)**2
        enddo

! Calculate F derivative
        do i=1, nactive
            Fderiv(i) = 2.*(0.5*indconduc(i,i)*sigma_energy*sigma_coils(i)*curdiff(i) + &
                sigma_B*sum((psicorr - x1)*(G_00(i, 1:nteta) - G_00c(i))) + &
                sigma_axis*(x2*G_00r(i) + x3*G_00z(i)) )
        enddo

        write(*, *) 'fderiv', fderiv

! Calculate new currents

        do i=1, nactive
            curnow(i) = curnow(i) - sum(invmatrix(i, 1:nactive)*Fderiv)
        enddo
        curdiff = curnow - curref
        write(*, *) 'result vector', curdiff(1:nactive)

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

        call find_new_axis_part1
        call find_psi_boundary
        call new_jrz_feqis  ! calculate new right hand side

        temp_err = sum(abs(Fderiv))
        Ffunc_old = Ffunc
        write(*,*) 'iteration ',j_iter, curnow(1:nactive),rax,zax,temp_err,sum(abs(Fderiv(1:nactive)))

        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nactive
        curconduc(i) = curnow(i)
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis_part1
    call find_psi_boundary
    call new_jrz_feqis  ! calculate new right hand side

    write(*, *) ' '
    write(*, *) curconduc(1:nconduc), rax, zax
    write(*, *) 'full fonfit active currents converged', temp_err, err_find_psistab

    return
    end subroutine restab_F_function_full_currents_forces

!--------------------------------------------------------------------
    subroutine restab_boundary_with_fourier_wall !not working well

    use errors_params, only: err_find_psistab
    use astra2fbe, only: n_fourier_restab_boundary
    use green_matrix, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, pol_angle

    integer :: i, j, k, j_iter, iax, jax, info
    double precision :: temp_err, curr, f_correction, x1, psibt0, psibt1
    double precision, dimension(9) :: bub
    double precision, dimension(npassive) :: anglr
    double precision, dimension(nteta) :: psicorr
    double precision, dimension(n_fourier_restab_boundary) :: S_00, C_00
    double precision, dimension(4*nteta*n_fourier_restab_boundary) :: work
    double precision, dimension(nr2, nz2) :: g
    double precision, dimension(nteta, n_fourier_restab_boundary) :: G_00c, G_00s
    double precision, dimension(nteta, 2*n_fourier_restab_boundary) :: matrix

! First, initialized initial guess coming from prescribed boundary current density: jrhoteta
    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, r(1), dr)
    jaxis = closest_index(zax, z(1), dz)
    iax = iaxis
    jax = jaxis

    G_00c = 0.
    G_00s = 0.

! Evaluate coils' quantities
    do i=1, npassive
        anglr(i) = pol_angle(rax, zax, r_cond(nactive + i), z_cond(nactive + i))
! Find true axis
        do j=1, nteta
            bub(2) = interp2d_psi(rbndp(j), zbndp(j), r(1:nr), z(1:nz), greeni(1:nr, 1:nz, nactive + i))
            do k=1, n_fourier_restab_boundary
                G_00c(j, k) = G_00c(j, k) + cos(k*anglr(i))*bub(2)
                G_00s(j, k) = G_00s(j, k) + sin(k*anglr(i))*bub(2)
            enddo
        enddo
    enddo

    do j=1, nteta
        do k=1, n_fourier_restab_boundary
            matrix(j, k) = G_00c(j, k)
            matrix(j, n_fourier_restab_boundary + k) = G_00s(j, k)
        enddo
    enddo

    psibt0 = 1000.
    psibt1 = 1000.

    do j_iter=1, 30
        g = 0.
        call solve_gs2d(g) !jrz as right hand side
        g = boundary(g)   ! gbound = integral (Green*dg/dn) over the boundary
        call solve_gs2d(g) ! again jrz as right hand side
        psiplasrz(1:nr2, 1:nz2) = g(1:nr2, 1:nz2)
        call compound_psi
        do j=1, nteta
            psicorr(j) = interp2d_psi(rbndp(j), zbndp(j), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        enddo

        x1 = sum(psicorr)/(nteta + 0.)
        psibt1  = x1
        psicorr = x1 - psicorr

! Least square fit solution

        call dgels('N', nteta, n_fourier_restab_boundary*2, 1, matrix, nteta, &
            psicorr, nteta, WORK, 2*(nteta)*n_fourier_restab_boundary*2, INFO)

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

        call find_new_axis_part1
        call find_psi_boundary
        call new_jrz_feqis  ! calculate new right hand side

        temp_err = abs(psibt0 - psibt1)
        psibt0 = psibt1
        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, npassive
        do j=1, n_fourier_restab_boundary
            curconduc(nactive + i) = curconduc(nactive + i) + &
                C_00(j)*cos(j*anglr(i)) + S_00(j)*sin(j*anglr(i))
        enddo
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis_part1
    call find_psi_boundary
    call new_jrz_feqis  ! calculate new right hand side

    write(*, *) curconduc(1:nconduc), rax, zax

    stop
    return
    end subroutine restab_boundary_with_fourier_wall

!--------------------------------------------------------------------
    subroutine restab_axis_with_fourier_wall

    use errors_params, only: err_find_psistab
    use green_matrix, only: greeni
    use feqis_tools, only: closest_index, interp2d_psi, &
        pol_angle, least_square_biquad

    integer :: i, j, j_iter, iax, jax
    double precision :: curr, dum1, dum2, zum1, psistab1o, psistab2o, delr, delz, &
        S_00r, C_00r, S_00z, C_00z, temp_err
    double precision, dimension(6) :: ccc
    double precision, dimension(8) :: ddipsi
    double precision, dimension(9) :: bub, xub, yub
    double precision, dimension(npassive) :: anglr, g0_r, g0_z
    double precision, dimension(258, 258) :: C_00, S_00
    double precision, dimension(nr2, nz2) :: g

!first, initialized initial guess coming from prescribed boundary current density: jrhoteta
    write(*, *) 'reinterp curr, restab'
    call interp_j_fromrhotorz

! Rescale current density
    curr = SUM(jrz(1:nr2, 1:nz2)) *dr*dz
    jrz = jrz/curr*iplasma

    rax = raxp
    zax = zaxp
    iaxis = closest_index(rax, r(1), dr)
    jaxis = closest_index(zax, z(1), dz)
    iax = iaxis
    jax = jaxis

! Evaluate coils' quantities
    do i=1, npassive
        anglr(i) = pol_angle(rax, zax, r_cond(nactive + i), z_cond(nactive + i))
! Find true axis
        xub(1) = r(iax-1)
        xub(2) = r(iax)
        xub(3) = r(iax + 1)
        xub(4) = r(iax)
        xub(5) = r(iax)
        xub(6) = r(iax-1)
        xub(7) = r(iax-1)
        xub(8) = r(iax + 1)
        xub(9) = r(iax + 1)
        yub(1) = z(jax)
        yub(2) = z(jax)
        yub(3) = z(jax)
        yub(4) = z(jax-1)
        yub(5) = z(jax + 1)
        yub(6) = z(jax-1)
        yub(7) = z(jax + 1)
        yub(8) = z(jax-1)
        yub(9) = z(jax + 1)
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
        g = 0.
        call solve_gs2d(g) ! jrz as right hand side
        g = boundary(g)   ! gbound = integral (Green*dg/dn) over the boundary
        call solve_gs2d(g) ! again jrz as right hand side

        psiplasrz(1:nr2, 1:nz2) = g(1:nr2, 1:nz2)
        psistabr = 0.
        psistabz = 0.
        delr = 0.
        delz = 0.
        call compound_psi
        call find_new_axis_part1
        dum1 = C_00r*S_00z - C_00z*S_00r

        bub(1) = interp2d_psi(raxp + dr, zaxp, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(2) = interp2d_psi(raxp - dr, zaxp, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(3) = interp2d_psi(raxp, zaxp + dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
        bub(4) = interp2d_psi(raxp, zaxp - dz, r(1:nr), z(1:nz), psirz(1:nr, 1:nz))

        xub(1) = (bub(1) - bub(2))/(2.*dr)
        yub(1) = (bub(3) - bub(4))/(2.*dz)
        delr = ( S_00z*(-xub(1)) - C_00z*(-yub(1)))/dum1
        delz = (-S_00r*(-xub(1)) + C_00r*(-yub(1)))/dum1
        psistabr = delr
        psistabz = delz

        call compound_psi
        psirz(1:nr2, 1:nz2) = psirz(1:nr2, 1:nz2) +  &
            psistabr*C_00(1:nr2, 1:nz2) + psistabz*S_00(1:nr2, 1:nz2) !total flux

        call find_new_axis_part1
        call find_psi_boundary
        call new_jrz_feqis  ! calculate new right hand side

        temp_err = (abs(psistab1o - psistabr) + abs(psistab2o - psistabz))
        psistab1o = psistabr
        psistab2o = psistabz

        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, npassive
        curconduc(nactive + i) = curconduc(nactive + i) + &
            psistabr*cos(anglr(i)) + psistabz*sin(anglr(i))
    enddo

    call psi_external_calc
    call compound_psi
    call find_new_axis_part1
    call find_psi_boundary
    call new_jrz_feqis  ! calculate new right hand side

    write(*, *) 'Code restab-axis with Fourier wall converged'

    return
    end subroutine restab_axis_with_fourier_wall

!--------------------------------------------------------------------
    subroutine find_new_axis_part1

    integer :: j, iax, jax, i1, i2, j1, j2
    double precision :: errtol, tolerr, raxm, zaxm
    double precision, dimension(2) :: ppx

! 1) find new magnetic axis

    iax = iaxis
    jax = jaxis
    i1 = 10000
    j1 = 10000
    i2 = 10000
    j2 = 10000
    errtol = 1.e6
    tolerr = 10.
    j = 1
    raxm = rax
    zaxm = zax

    do while(errtol > tolerr)
        i1 = i2
        j1 = j2
        i2 = iax
        j2 = jax
        if (psirz(iax + 1, jax) > psirz(iax, jax)) then
            iax = iax + 1
            jax = jax
        endif
        if (psirz(iax-1, jax) > psirz(iax, jax)) then
            iax = iax - 1
            jax = jax
        endif
        if (psirz(iax, jax + 1) > psirz(iax, jax)) then
            iax = iax
            jax = jax + 1
        endif
        if (psirz(iax, jax-1) > psirz(iax, jax)) then
            iax = iax
            jax = jax - 1
        endif
        if (psirz(iax + 1, jax-1) > psirz(iax, jax)) then
            iax = iax + 1
            jax = jax - 1
        endif
        if (psirz(iax-1, jax + 1) > psirz(iax, jax)) then
            iax = iax - 1
            jax = jax + 1
        endif
        if (psirz(iax-1, jax-1) > psirz(iax, jax)) then
            iax = iax - 1
            jax = jax - 1
        endif
        if (psirz(iax + 1, jax + 1) > psirz(iax, jax)) then
            iax = iax + 1
            jax = jax + 1
        endif

        j = j + 1

        if ((iax == i1) .and. (jax == j1)) EXIT
        if (j >= 100000) EXIT
    enddo

    iaxis = iax
    jaxis = jax
    call nine_point_regression(r(iax), z(jax), ppx, derivpsi, psiaxis)

    rax = ppx(1)  !r(iaxis)
    zax = ppx(2) !z(jaxis)
    trax = rax
    tzax = zax

    if (isnan(rax)) then
        write(*, *) 'rax is nan in fbe find axis'
        stop
    endif
    if (rax < 0.1 .or. rax > 100) then
        write(*, *) 'rax is nan in fbe find axis', rax, zax, ppx, derivpsi, psiaxis
        stop
    endif

    return
    end subroutine find_new_axis_part1

!-------------------------------------------------------------------
    subroutine compound_psi

    implicit none

    psirz = psiplasrz + psiextrz

    if (nferromag >= 1) then
        psiferro = 0.
        call ferro_mag_create
        psirz = psirz + psiferro
    endif

    return
    end subroutine compound_psi
  
!-------------------------------------------------------------------
    subroutine ferro_mag_create

    use ferromagstructure, only: type_ferromag
    use feqis_tools, only: interp2d_psi, green_function, inv_matrix
    use numerical_tools, only: qinterp
    use pi_vars, only: GPI, GPI2, GPI4, muvac

    implicit none

    integer :: i, j, ii, jj, iii, jjj, iferro, nval
    double precision :: x1, x2, x3, x4, d
    double precision, dimension(1) :: z1, z2
    type(type_ferromag), dimension(:), allocatable :: ferromag

! Loop over ferromagnetic elements
    do iferro=1, nferromag
        iii  = ferromag(iferro)%position%npoints
        nval = ferromag(iferro)%mhrelation%nvalues
! Construct vacuum field
        do j=1, iii
            x1 = interp2d_psi(ferromag(iferro)%position%r(j) + dr/2, ferromag(iferro)%position%z(j), r, z, psirz)
            x2 = interp2d_psi(ferromag(iferro)%position%r(j) - dr/2, ferromag(iferro)%position%z(j), r, z, psirz)
            x3 = interp2d_psi(ferromag(iferro)%position%r(j), ferromag(iferro)%position%z(j) + dz/2, r, z, psirz)
            x4 = interp2d_psi(ferromag(iferro)%position%r(j), ferromag(iferro)%position%z(j) - dz/2, r, z, psirz)

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
                    d = abs(r(ii) - ferromag(iferro)%position%r(j)) + abs(z(ii) - ferromag(iferro)%position%z(j))
                   if (d > 0) then
                       psiferro(ii, jj) = psiferro(ii, jj) + muvac/GPI*green_function(r(ii), z(jj), ferromag(iferro)%position%r(j), ferromag(iferro)%position%z(j))*ferromag(iferro)%position%current(j)
                   endif
               enddo
           enddo
       enddo
    enddo

    return
    end subroutine ferro_mag_create

!--------------------------------------------------------------------
    subroutine find_psi_boundary

    use pi_vars, only: GPI
    use astra2fbe, only: x_point_save, plasma_config
    use errors_params, only: err_find_oxpoints_derivs
    use feqis_tools, only: closest_index, pol_angle, &
        interp2d_psi

    integer :: iaold, niter, i, j, k, i1, i4, i5, i9, n_adding, & ! oldpointnum,
        i_county
    double precision :: x1, x2, x5
    double precision, dimension(2) :: pos_xpoint(2)
    double precision, dimension(8) :: ddipsi
    double precision, dimension(200) :: rx_add, zx_add
    double precision, dimension(500) :: psi_limp
    double precision, dimension(max_xpoints) :: psi_xpoint

    data i_county/0/
    save i_county  !, oldpointnum

    i_plasmatype = 0

    if (n_of_xpoints == 0) i_county = 0 !reset to full search if there are no x points!

    if (n_of_xpoints >= 1) then
        x_point_save(1:n_of_xpoints, 1) = r_xpoint(1:n_of_xpoints)
        x_point_save(1:n_of_xpoints, 2) = z_xpoint(1:n_of_xpoints)
    endif

    call find_closest_xpoints(rx_add, zx_add, i9, n_adding)
    if (i9 == 0) then
        r_xpoint(n_of_xpoints + 1:n_of_xpoints + n_adding) = rx_add(1:n_adding)
        z_xpoint(n_of_xpoints + 1:n_of_xpoints + n_adding) = zx_add(1:n_adding)
        n_of_xpoints = n_of_xpoints + n_adding
    endif

    if (i_county == 1) then
!go through old x-points and see where they end up
        if (n_of_xpoints >= 1) then
            iaold = n_of_xpoints
            do i=1, iaold
                do niter=1, 101
                    call nine_point_regression_follow(r_xpoint(i), z_xpoint(i), pos_xpoint, ddipsi, x1)
                    r_xpoint(i) = pos_xpoint(1)
                    z_xpoint(i) = pos_xpoint(2)
                    if ( (abs(ddipsi(1)) + abs(ddipsi(2))) <= err_find_oxpoints_derivs) then
! Check if point outside of domain
                        if ((pos_xpoint(1) > r(nr2) - dr) .or. (pos_xpoint(1) < r(1) + dr) .or.  &
                            (pos_xpoint(2) > z(nz2) - dz) .or. (pos_xpoint(2) < z(1) + dz)) then
! xpoint doesnt exist anymore
                            r_xpoint(i) = 1.e6
                            z_xpoint(i) = 0.
                        else if ((pos_xpoint(1) > rax - dr) .and. (pos_xpoint(1) < rax + dr) .and.  &
                                 (pos_xpoint(2) > zax - dz) .and. (pos_xpoint(2) < zax + dz)) then
! xpoint doesn't exist anymore
                            r_xpoint(i) = 1.e6
                            z_xpoint(i) = 0.
                        else
                            r_xpoint(i) = pos_xpoint(1)
                            z_xpoint(i) = pos_xpoint(2)
                        endif
                        EXIT
                    endif
                    if (niter > 100) then
                        r_xpoint(i) = 1.e6
                        z_xpoint(i) = 0.
                        EXIT
                    endif
                    if ((pos_xpoint(1) > r(nr2) - dr) .or. (pos_xpoint(1) < r(1) + dr) .or.  &
                        (pos_xpoint(2) > z(nz2) - dz) .or. (pos_xpoint(2) < z(1) + dz)) then
! xpoint doesn't exist anymore
                        r_xpoint(i) = 1.e6
                        z_xpoint(i) = 0.
                        EXIT
                    endif
                enddo
            enddo
        endif

! Scan the boundary to find new x-points
        do j=2, nz1, nz1-2
            do i=2, nr1
                call nine_point_regression(r(i), z(j), pos_xpoint, ddipsi, x1)
                x5 = (ddipsi(5)**2 - ddipsi(3)*ddipsi(4))
                if ((pos_xpoint(1) >= r(i) - dr) .and.  &
                    (pos_xpoint(1) <= r(i) + dr) .and.  &
                    (pos_xpoint(2) >= z(j) - dz) .and.  &
                    (pos_xpoint(2) <= z(j) + dz) .and.  &
                    (x5 >= 0.)) then

                    n_of_xpoints = min(max_xpoints, n_of_xpoints + 1)
                    r_xpoint(n_of_xpoints) = pos_xpoint(1)
                    z_xpoint(n_of_xpoints) = pos_xpoint(2)
                endif
            enddo
        enddo

        do i=2, nr1, nr1-2
            do j=2, nz1
                call nine_point_regression(r(i), z(j), pos_xpoint, ddipsi, x1)
                x5 = (ddipsi(5)**2 - ddipsi(3)*ddipsi(4))
                if ((pos_xpoint(1) >= r(i) - dr) .and.  &
                    (pos_xpoint(1) <= r(i) + dr) .and.  &
                    (pos_xpoint(2) >= z(j) - dz) .and.  &
                    (pos_xpoint(2) <= z(j) + dz) .and.  &
                    (x5 >= 0.)) then

                    n_of_xpoints = min(max_xpoints, n_of_xpoints + 1)
                    r_xpoint(n_of_xpoints) = pos_xpoint(1)
                    z_xpoint(n_of_xpoints) = pos_xpoint(2)
                endif
            enddo
        enddo

    endif

!-------------------
    if (i_county == 0) then ! do a full pass to find all X-points
        n_of_xpoints = 0
        do j=2, nz1
            do i=2, nr1
                call nine_point_regression(r(i), z(j), pos_xpoint, ddipsi, x1)
                x5 = (ddipsi(5)**2 - ddipsi(3)*ddipsi(4))

                if ((pos_xpoint(1) >= r(i) - dr) .and.  &
                    (pos_xpoint(1) <= r(i) + dr) .and.  &
                    (pos_xpoint(2) >= z(j) - dz) .and.  &
                    (pos_xpoint(2) <= z(j) + dz) .and.  &
                    (x5 >= 0.)) then

                    if (abs(pos_xpoint(1) - rax) > 2.*dr .or. abs(pos_xpoint(2) - zax) > 2.*dz) then
                        n_of_xpoints = min(max_xpoints, n_of_xpoints + 1)
                        r_xpoint(n_of_xpoints) = pos_xpoint(1)
                        z_xpoint(n_of_xpoints) = pos_xpoint(2)
                    endif
                endif
            enddo
        enddo
        i_county = 1
    endif

! Remove disappeared x-points and double counts

    i = 0
    xpoints_loop: do
        i = i + 1
        if (i > n_of_xpoints) EXIT xpoints_loop
        if (r_xpoint(i) >= 1.e5) then
            if (n_of_xpoints == 1) then
                n_of_xpoints = 0
                EXIT xpoints_loop
            endif
            r_xpoint(i:n_of_xpoints-1) = r_xpoint(i+1:n_of_xpoints)
            z_xpoint(i:n_of_xpoints-1) = z_xpoint(i+1:n_of_xpoints)
            n_of_xpoints = n_of_xpoints - 1
            i = i - 1
            CYCLE xpoints_loop
        endif

        do k=1, i-1 ! check if double counted
            if ((abs(r_xpoint(i) - r_xpoint(k)) <= 2.*dr) .and.  &
                (abs(z_xpoint(i) - z_xpoint(k)) <= 2.*dz)) then

                r_xpoint(k) = 0.5*(r_xpoint(i) + r_xpoint(k))
                z_xpoint(k) = 0.5*(z_xpoint(i) + z_xpoint(k))
                r_xpoint(i:n_of_xpoints-1) = r_xpoint(i+1:n_of_xpoints)
                z_xpoint(i:n_of_xpoints-1) = z_xpoint(i+1:n_of_xpoints)
                n_of_xpoints = n_of_xpoints - 1
                i = i - 1
                CYCLE xpoints_loop
            endif
        enddo
    enddo xpoints_loop


! Ignore limiter if use_limiter_astra is 0, da trasferirsi in init
    if (use_limiter_yesno == 0) then
        limiterR = r(nr1)
        limiterZ = z(nz1)
    endif

! Calculate limiter flux
    do i=1, nlimiter
        psi_limp(i) = interp2d_psi(limiterR(i), limiterZ(i), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
    enddo

    if (n_of_xpoints == 0) then ! No x-points, take largest limiter flux
        psibnd = maxval(psi_limp(1:nlimiter), 1)
        i_plasmatype = 0
    endif

! now, remove limiters that are in the shadow of xpoints
    ztop   =  1.e6
    zbot   = -1.e6
    raus   =  1.e6
    rinner =  0.

    if (n_of_xpoints >= 1) then
        i_plasmatype = 1

! First pass, remove X-points behind the limiter area
        do i=1, n_of_xpoints
            j = closest_index(r_xpoint(i), r(1), dr)
            k = closest_index(z_xpoint(i), z(1), dz)
            if (zlimpotential(j, k) == 0) then
                psi_xpoint(i) = -1.e6
            else
                psi_xpoint(i) = interp2d_psi(r_xpoint(i), z_xpoint(i), r(1:nr), z(1:nz), psirz(1:nr, 1:nz))
            endif
        enddo

! Second pass, remove limiter points that are in x-points shadow, simple "straight line method" --> to be refined later on
        do i=1, n_of_xpoints
            if (psi_xpoint(i) > -1.e5) then
                x1 = pol_angle(rax, zax, r_xpoint(i), z_xpoint(i))
                if ((r_xpoint(i) > rax) .and. (x1 >= 7./4.*GPI .or.  x1 <= GPI/4.   )) raus   = min(raus, r_xpoint(i))
                if ((z_xpoint(i) > zax) .and. (x1 >=   GPI/4.  .and. x1 <= 3./4.*GPI)) ztop   = min(ztop, z_xpoint(i))
                if ((r_xpoint(i) < rax) .and. (x1 >= 3./4.*GPI .and. x1 <= 5./4.*GPI)) rinner = max(rinner, r_xpoint(i))
                if ((z_xpoint(i) < zax) .and. (x1 >= 5./4.*GPI .and. x1 <= 7./4.*GPI)) zbot   = max(zbot, z_xpoint(i))
            endif
        enddo
        do j=1, nlimiter
            if (limiterr(j) <= rinner) psi_limp(j) = -1.e6
            if (limiterr(j) >= raus) psi_limp(j) = -1.e6
            if (limiterz(j) <= zbot) psi_limp(j) = -1.e6
            if (limiterz(j) >= ztop) psi_limp(j) = -1.e6
        enddo

! Third pass, remove x-points that are non-monotonically connected to the plasma.
        do i=1, n_of_xpoints
            if (psi_xpoint(i) > -1.e5) then
                i1 = xpoint_axis_connection(r_xpoint(i), z_xpoint(i), rax, zax, dr, dz)
                if (i1 == 0) psi_xpoint(i) = -1.e6
            endif
        enddo

! Psi boundary is the highest of all the values
        i4 = maxloc(psi_xpoint(1:n_of_xpoints), 1)
        i5 = maxloc(psi_limp(1:nlimiter), 1)
        x1 = psi_xpoint(i4)
        x2 = psi_limp(i5)
        psibnd = max(x1, x2)
        if (x2 > x1) then
            i_plasmatype = 0
        else
            i_plasmatype = 1
        endif
    endif

    plasma_config = i_plasmatype
    x_point_save(20, 1) = r_xpoint(i4)
    x_point_save(20, 2) = z_xpoint(i4)

! if (i_plasmatype == 1)
    psibnd = psiaxis + (psibnd - psiaxis)*alpsep

! Normalized flux
    u_n(1:nr2, 1:nz2) = (psirz(1:nr2, 1:nz2) - psiaxis)/(psibnd - psiaxis)

    return
    end subroutine find_psi_boundary

!--------------------------------------------------------------------
    subroutine new_jrz_feqis ! calculate new right hand side given new boundary!

    use rcurr_zcurr_2def, only: R_curr_2d, Z_curr_2D
    use feqis_tools, only: fill_in_current, floor_index

    integer :: i, j, i1, i2, j1, quadrant, ipluz, jpluz, &
        ilast, totpoints, istart, j_griddo_j, i_griddo_j
    integer, dimension(nr2*nz2, 2) :: external_griddo_j, internal_griddo
    double precision :: curr, darea, t1, t2, t3, t4, je1, je2, je3, je4, &
        z11, z12, z13, z14
    double precision, dimension(nr2, nz2) :: iconvex
    double precision, dimension(nr2, nz2) :: dumc

! in entry: rbnd, zbnd, nbnd, psiaxis, psibnd, u_n

    iconvex = 0.
    dumc    = 0.
    area_eff = dr*dz
    darea = dr*dz
    i_griddo_j = 0
    j_griddo_j = 0

    totpoints = 0
    quad_loop: do quadrant=1, 4
! sweep from axis to exterior and fill in the current
!start from axis position

        i1 = floor_index(rax, r(1), dr)
        j1 = floor_index(zax, z(1), dz)

        SELECT CASE(quadrant)
        CASE(1)
            i = i1 + 1
            j = j1 + 1
            ipluz = 1
            jpluz = 1
        CASE(2)
            i = i1
            j = j1 + 1
            ipluz = -1
            jpluz =  1
        CASE(3)
            i = i1
            j = j1
            ipluz = -1
            jpluz = -1
        CASE(4)
            i = i1 + 1
            j = j1
            ipluz =  1
            jpluz = -1
        END SELECT

        istart = i

        inner_loop: do

            totpoints = totpoints + 1
            if (totpoints > nz2*nr2) then
                write(*, *) 'Error in find new boundary (totpoints > nz2*nr2)', totpoints, nr2, nz2
                stop
            endif

            dumc(i, j) = fill_in_current(r(i), nrho2d, ppp_2d, ffp_2d, u_n(i, j))
            iconvex(i, j) = 1.
            i_griddo_j = i_griddo_j + 1
            internal_griddo(i_griddo_j, 1) = i
            internal_griddo(i_griddo_j, 2) = j
            i1 = i
            i2 = j

            t1 = t_find_u_n(i1 + ipluz, i2, i1, i2)
            t2 = t_find_u_n(i1, i2 + jpluz, i1, i2)
            if (t1 < 0 .or. t1 > 1.) t1 = 1.e6
            if (t2 < 0 .or. t2 > 1.) t2 = 1.e6
            if (t2 > 0. .and. t2 <= 1.) then
                j_griddo_j = j_griddo_j + 1
                external_griddo_j(j_griddo_j, 1) = i
                external_griddo_j(j_griddo_j, 2) = j + jpluz
                iconvex(i, j + jpluz) = 1.
            endif

            if (t1 > 1.e5) then
! Move horizontally to the right
                i = i + ipluz
            else
                ilast = i + ipluz
                j_griddo_j = j_griddo_j + 1
                external_griddo_j(j_griddo_j, 1) = ilast
                external_griddo_j(j_griddo_j, 2) = j
                iconvex(ilast, j) = 1.
! Found boundary, go back, check vertically
                i = istart
                do
                    t2 = t_find_u_n(i, j + jpluz, i, j)
                    if (t2 < 0 .or. t2 > 1.) then
                        j = j + jpluz
                        CYCLE inner_loop
                    endif

! Found boundary on Z, need to advance 1 more
                    j_griddo_j = j_griddo_j + 1
                    external_griddo_j(j_griddo_j, 1) = i
                    external_griddo_j(j_griddo_j, 2) = j + jpluz
                    iconvex(i, j + jpluz) = 1.

                    i = i + ipluz
                    istart = i
                    if (i == ilast) CYCLE quad_loop
                enddo
                EXIT inner_loop
            endif
        enddo inner_loop
    enddo quad_loop

! Trick at boundary for ciurrent

! fill current in external griddo
    do i=1, j_griddo_j
        t1 = 0.
        t2 = 0.
        t3 = 0.
        t4 = 0.
        i1 = external_griddo_j(i, 1)
        i2 = external_griddo_j(i, 2)
        t1 = t_find_u_n(i1 - 1, i2, i1, i2)
        t2 = t_find_u_n(i1, i2 + 1, i1, i2)
        t3 = t_find_u_n(i1 + 1, i2, i1, i2)
        t4 = t_find_u_n(i1, i2 - 1, i1, i2)
        if (t1 < 0 .or. t1 > 1.) t1 = 0.
        if (t2 < 0 .or. t2 > 1.) t2 = 0.
        if (t3 < 0 .or. t3 > 1.) t3 = 0.
        if (t4 < 0 .or. t4 > 1.) t4 = 0.

        je1 = 0.
        je2 = 0.
        je3 = 0.
        je4 = 0.
        z11 = r(i1-1)*(1 - t1) + r(i1)*t1
        z12 = r(i1)
        z13 = r(i1 + 1)*(1 - t3) + r(i1)*t3
        z14 = r(i1)
        if (t1 > 0.) je1 = fill_in_current(z11, nrho2d, ppp_2d, ffp_2d, 1.d0)
        if (t2 > 0.) je2 = fill_in_current(z12, nrho2d, ppp_2d, ffp_2d, 1.d0)
        if (t3 > 0.) je3 = fill_in_current(z13, nrho2d, ppp_2d, ffp_2d, 1.d0)
        if (t4 > 0.) je4 = fill_in_current(z14, nrho2d, ppp_2d, ffp_2d, 1.d0)

! defining S1 = dR - d1, S2 = dZ - d2, S3 = dZ - d3, S4 = dR - d4
! Jvacuum = C*Jb
! if 1 point only: C = 1 - S1/dR
! if 2 points: C = 1 - S1*S2/(dR*dZ)
! if 3 points: C = 1 - (S1*S2 + S1*S3)/(2*dR*dZ)
! if 4 points: C = 1 - (S1*S2 + S1*S3 + S3*S4 + S4*S2)/(4*dR*dZ)
! t_j = d_j /(dR or dZ depending on direction)

        dumc(i1, i2) = t1*je1 + t2*je2 + t3*je3 + t4*je4 -  &
            (t1*t2*(je1 + je2)/2. +  &
             t1*t3*(je1 + je3)/2. +  &
             t1*t4*(je1 + je4)/2. +  &
             t2*t3*(je3 + je2)/2. +  &
             t2*t4*(je4 + je2)/2. +  &
             t3*t4*(je3 + je4)/2.)

        iconvex(i1, i2) = t1 + t2 + t3 + t4 -  &
            (t1*t2 + t1*t3 + t1*t4 + t2*t3 + t2*t4 + t3*t4)

    enddo

    jrz = 0.
    do j=2, nz1
        do i=2, nr1
            jrz(i, j) = iconvex(i, j)*0.5*(dumc(i, j) + 0.25*(dumc(i+1, j) + dumc(i-1, j) + dumc(i, j-1) + dumc(i, j+1))) ! this is a first order Shapiro filter with alpha = 0.5, with in addition the iconvex multiplier which reduces the ghost currents even more to avoid overshooting. Seems to work well vs spider
        enddo
    enddo

!uncomment below for consistent current
!jrz(1:nr2, 1:nz2)=dumc(1:nr2, 1:nz2)
    curr = sum(jrz)*darea
    jrz = jrz/curr*iplasma

    if (isnan(curr)) then
        write(*, *) 'Total current is nan in fbe current rescaling'
        stop
    endif

    t1 = 0.
    t2 = 0.
    t3 = 0.
    do j=1, nz2
         do i=1, nr2
            t1 = t1 + r(i)**2*jrz(i, j)*darea
            t2 = t2 + z(j)*jrz(i, j)*darea
            t3 = t3 + jrz(i, j)*darea
        enddo
    enddo

    R_curr_2D =  sqrt(t1/t3)
    Z_curr_2D =  t2/t3

    return
    end subroutine new_jrz_feqis

!---------------------------------------------------------------------
    subroutine interp_j_fromrhotorz

    use feqis_tools, only: curinterp

    integer :: i, j, k, k1, k2
    double precision :: t1, t2, t3, t4

! go from jrhoteta to jrz
    jrz = 0.
    jrhoteta(1:nrho, nteta+1) = jrhoteta(1:nrho, 1)
    do j=1, nz2
        do i=1, nr2
            jrz(i, j) = curinterp(r(i), z(j), jrhoteta(1:nrho, 1:nteta+1),  &
                rho(1:nrho, 1:nteta+1), teta(1:nteta+1), raxp, zaxp, nrho, nteta+1)
        enddo
    enddo

    return
    end subroutine interp_j_fromrhotorz

!---------------------------------------------------------------------
    double precision function t_find_u_n(i1, j1, i2, j2)

    integer, intent(in) :: i1, i2, j1, j2

    t_find_u_n = (1. - u_n(i1, j1))/(u_n(i2, j2) - u_n(i1, j1))

    return
    end function t_find_u_n


end module feqis_circuit
