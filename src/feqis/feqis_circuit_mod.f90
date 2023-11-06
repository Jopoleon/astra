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
    anglecoil, anglehcoil

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
    double precision, dimension(300, 300) :: g
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
        if (j_iter > 50) stop
        g = 0.
        call solve_gs2d(g) !jrz as right hand side
        g = boundary(g)  ! gbound = integral (Green*dg/dn) over the boundary

        call solve_gs2d(g) ! again jrz as right hand side
        psiplasrz(1:nr2, 1:nz2) = g(1:nr2, 1:nz2)

! Construct correction
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nconduc
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k) 
                enddo
                psirz(i, j) = psiplasrz(i, j) + psiextrz(i, j) + f_correction !total flux
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
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, nconduc
                    f_correction = f_correction + curdiff(k)*greeni(i, j, k) 
                enddo
                psirz(i, j) = psiplasrz(i, j) + psiextrz(i, j) + f_correction !total flux
            enddo
        enddo

        call find_new_axis_part1
        call find_psi_boundary
        call new_jrz_feqis  ! calculate new right hand side

        temp_err = abs(Ffunc - Ffunc_old)
        Ffunc_old = Ffunc
        if (temp_err <= err_find_psistab) EXIT

    enddo

    do i=1, nconduc
        curconduc(i) = curnow(i)
    enddo

    call psi_external_calc
    psirz(1:nr2, 1:nz2) = psiplasrz(1:nr2, 1:nz2) + psiextrz(1:nr2, 1:nz2)
    call find_new_axis_part1
    call find_psi_boundary
    call new_jrz_feqis  ! calculate new right hand side

    write(*, *) curconduc(1:nconduc), rax, zax
    write(*, *) 'full fonfit eddy currents converged'

    return
    end subroutine restab_F_function_full_fonfit

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
    double precision, dimension(300, 300) :: g
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
        psirz(1:nr2, 1:nz2) = psiplasrz(1:nr2, 1:nz2) + psiextrz(1:nr2, 1:nz2) !total flux
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
        do j=1, nz2
            do i=1, nr2
                f_correction = 0.
                do k=1, n_fourier_restab_boundary
                    f_correction = f_correction +  & 
                        psicorr(k)                            *sum(greeni(i, j, nactive + 1:nactive + npassive)*cos(k*anglr(1:npassive))) + & 
                        psicorr(n_fourier_restab_boundary + k)*sum(greeni(i, j, nactive + 1:nactive + npassive)*sin(k*anglr(1:npassive)))
                enddo
                psirz(i, j) = psiplasrz(i, j) + psiextrz(i, j) + f_correction !total flux
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
    psirz(1:nr2, 1:nz2) = psiplasrz(1:nr2, 1:nz2) + psiextrz(1:nr2, 1:nz2)
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
    double precision, dimension(300, 300) :: g

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
        psirz(1:nr2, 1:nz2) = psiplasrz(1:nr2, 1:nz2) + psiextrz(1:nr2, 1:nz2) !total flux
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

        psirz(1:nr2, 1:nz2) = psiplasrz(1:nr2, 1:nz2) + psiextrz(1:nr2, 1:nz2) +  &
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
    psirz(1:nr2, 1:nz2) = psiplasrz(1:nr2, 1:nz2) + psiextrz(1:nr2, 1:nz2)
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

!--------------------------------------------------------------------
    subroutine find_psi_boundary

    use pi_vars, only: GPI
    use astra2fbe, only: x_point_save, plasma_config
    use errors_params, only: err_find_oxpoints_derivs
    use feqis_tools, only: closest_index, pol_angle, &
        interp2d_psi

    integer :: iaold, niter, i, j, k, oldpointnum, i1, i4, i5, i9, n_adding, &
        i_county
    double precision :: x1, x2, x5
    double precision, dimension(2) :: pos_xpoint(2)
    double precision, dimension(8) :: ddipsi
    double precision, dimension(20) :: rx_add, zx_add
    double precision, dimension(500) :: psi_limp
    double precision, dimension(max_xpoints) :: psi_xpoint

    data i_county/0/
    save i_county, oldpointnum

    i_plasmatype = 0

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
        oldpointnum = n_of_xpoints
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

    oldpointnum = n_of_xpoints

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
            if (zlimpotential(j, k) < 0.5) then
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
    integer, dimension(90000, 2) :: external_griddo_j, internal_griddo
    double precision :: curr, darea, t1, t2, t3, t4, je1, je2, je3, je4, &
        z11, z12, z13, z14
    double precision, dimension(300, 300) :: iconvex
    double precision, dimension(i_dim2, i_dim2) :: dumc

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
            jrz(i, j) = iconvex(i, j)*0.5*(dumc(i, j) + 0.25*(dumc(i+1, j) + dumc(i-1, j) + dumc(i, j-1) + dumc(i, j+1)))
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
