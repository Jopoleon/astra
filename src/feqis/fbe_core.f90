module fbe_core

implicit none

integer, parameter :: nrho2d=5999
integer :: use_limiter

integer :: nlimiter
double precision, dimension(:), allocatable :: limiterR, limiterZ
double precision :: lim_maxR, lim_minR, lim_maxZ, lim_minZ

!grids
integer :: nr, nz, nr2, nz2, nr1, nz1, nbnd, ngbnd, &
    redo_bnd
integer, dimension(:, :), allocatable :: zlimpotential

double precision :: rmin, rmax, zmin, zmax, dr, dz, dteta, &
    zbot, ztop, raus, rinner
double precision, dimension(:), allocatable :: Rrect, Zrect, rcomp, zcomp
double precision, dimension(:, :), allocatable :: area_eff

double precision, dimension(:, :), allocatable :: u_n, &
    psirz, psiextrz, psiplasrz, psiferro

double precision, dimension(nrho2d) :: psia_2d, ffp_2d, ppp_2d
double precision, dimension(8) :: derivpsi

!conductors
integer :: nconduc
double precision, dimension(:), allocatable :: curconduc

! boundary and axis FBE, PBE
integer, parameter :: max_xpoints=500
integer :: i_plasmatype !(0-limited, 1-single null, 2-double null)
integer :: iaxis, jaxis, n_of_xpoints, active_x_point
double precision :: psibnd, psiaxis, &
    rax, zax, trax, tzax, &
    alpsep, psistabR, psistabZ, dr_factor_init, dz_factor_init
double precision, dimension(max_xpoints) :: r_xpoint, z_xpoint
double precision :: deriv_x(5, max_xpoints)
double precision, dimension(:), allocatable :: green_bnd_f

! plasma parameters
double precision, dimension(:, :), allocatable :: jrz
double precision, dimension(:), allocatable :: rbnd, zbnd

contains

!--------------------------------------------------------------------
    subroutine solve_gs2d(g)

    use pi_vars, only: mu0
    use fft_mod_eff, only: costable
    use feqis_tools, only: discrete_sine_transform, solve_tridiag_fbe

    double precision, intent(inout), dimension(nr2, nz2) :: g

    integer :: i, j, k, j_init
    double precision :: x1, x2, r1m_1, r2m_1
    double precision, dimension(1000) :: A, B, C, z_fourier
    double precision, dimension(500, 500) :: gt, rhs, wrhs

    data j_init/0/
    save A, B, C, j_init, z_fourier

    do j=1, nz2
        do i=1, nr2
            rhs(i, j) = -mu0*Rrect(i)*jrz(i, j)
        enddo
    enddo

    r1m_1 = Rrect(  2)/dr**2/((Rrect(  1) + Rrect(  2))/2.)
    r2m_1 = Rrect(nr1)/dr**2/((Rrect(nr1) + Rrect(nr2))/2.)

    rhs(2:nr1,   2) = rhs(2:nr1,   2) - g(2:nr1,   1)/dz**2
    rhs(2:nr1, nz1) = rhs(2:nr1, nz1) - g(2:nr1, nz2)/dz**2

    rhs(  2, 2:nz1) = rhs(  2, 2:nz1) - g(  1, 2:nz1)*r1m_1
    rhs(nr1, 2:nz1) = rhs(nr1, 2:nz1) - g(nr2, 2:nz1)*r2m_1

    do i=2, nr1
        wrhs(i, 2:nz1) = discrete_sine_transform(nz, rhs(i, 2:nz1))
    enddo

! Create inverse matrix for gs2d
    if (j_init == 0) then
        A = 0.
        B = 0.
        C = 0.
        do i=2, nr-1
            x1 = 0.5*(rcomp(  i) + rcomp(i-1))
            x2 = 0.5*(rcomp(i+1) + rcomp(i  ))
            B(i) = -rcomp(i)/dr**2*(1./x2 + 1./x1)
            A(i) =  rcomp(i)/x2/dr**2
            C(i) =  rcomp(i)/x1/dr**2
        enddo
        i = 1
        x1 = 0.5*(rcomp(  i) + Rrect(1))
        x2 = 0.5*(rcomp(i+1) + rcomp(i))
        B(i) = -rcomp(i)/dr**2 * (1./x2 + 1./x1)
        A(i) =  rcomp(i)/x2/dr**2
        i = nr
        x1 = 0.5*(rcomp(i)   + rcomp(i-1))
        x2 = 0.5*(Rrect(nr2) + rcomp(i)  )
        B(i) = -rcomp(i)/dr**2 * (1./x2 + 1./x1)
        C(i) =  rcomp(i)/x1/dr**2
        j_init = 1
        do k=2, nz1
            z_fourier(k) = 2./dz**2 * (costable(k-1) - 1.)
        enddo
    endif

    gt = 0.
! Solve matrix
    do k=2, nz1
        gt(2:nr1, k) = solve_tridiag_fbe(C(1:nr), B(1:nr) + z_fourier(k), A(1:nr), wrhs(2:nr1, k), nr)
    enddo

! Invert fourier from gt(1:nr, 1:kfourier) to g(2:nr1, 2:nz1)
!  gt(i, k)=sum(invMM_gs2d(i-1, 1:nr, k-1)*wrhs(2:nr1, k))

    do i=2, nr1
        gt(i, 2:nz1) = discrete_sine_transform(nz, gt(i, 2:nz1))
    enddo
    g(2:nr1, 2:nz1) = 2./(nz + 1)*gt(2:nr1, 2:nz1)

    return
    end subroutine solve_gs2d

!---------------------------------------------------------------------
    subroutine psi_external_calc

    use green_function, only: greeni

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

    iax = closest_index(r_in, Rrect(1), dr)
    jax = closest_index(z_in, Zrect(1), dz)

    k = 0
    do j=-1, 1
        do i=-1, 1
            k = k + 1
            bub(k) = psirz(iax+i, jax+j)
        enddo
    enddo
    rax = Rrect(iax)
    zax = Zrect(jax)

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

    iax = closest_index(r_in, Rrect(1), dr)
    jax = closest_index(z_in, Zrect(1), dz)

    rax_out = Rrect(iax)
    zax_out = Zrect(jax)

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

    bub(1) = interp2d_psi(rx - dr, zx - dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    bub(2) = interp2d_psi(rx     , zx - dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    bub(3) = interp2d_psi(rx + dr, zx - dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    bub(4) = interp2d_psi(rx - dr, zx     , Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    bub(5) = interp2d_psi(rx     , zx     , Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    bub(6) = interp2d_psi(rx + dr, zx     , Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    bub(7) = interp2d_psi(rx - dr, zx + dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    bub(8) = interp2d_psi(rx     , zx + dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    bub(9) = interp2d_psi(rx + dr, zx + dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))

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

    if (rbnd(1) <= Rrect(1)) return ! boundary doesnt exist yet

    ierr = 0

    do i=1, nbnd !cycle over boundary points
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
        dgdn(j) = -green_in(j, 2)/dz*greenf*dr/Rrect(j)
    enddo
    bgint = bgint - sum(dgdn(2:nr1))

!right side
    do j=2, nz1
        jcounty = jcounty + 1
        greenf  = green_bnd_f(jcounty)
        dgdn(j) = -green_in(nr1, j)/dr*greenf*dz/(Rrect(nr2) + Rrect(nr1))*2.
    enddo
    bgint = bgint - sum(dgdn(2:nz1))

! upper side
    do j=2, nr1
        jcounty = jcounty + 1
        greenf  = green_bnd_f(jcounty)
        dgdn(j) = -green_in(j, nz1)/dz*greenf*dr/Rrect(j)
    enddo
    bgint = bgint - sum(dgdn(2:nr1))

!left side
    do j=2, nz1
        jcounty = jcounty + 1
        greenf  = green_bnd_f(jcounty)
        dgdn(j) = -green_in(2, j)/dr*greenf*dz/(Rrect(1) + Rrect(2))*2.
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
        z1 = interp2d_psi(t1, t2, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        z2 = interp2d_psi(t3, t4, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
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
    subroutine find_new_axis

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
    call nine_point_regression(Rrect(iax), Zrect(jax), ppx, derivpsi, psiaxis)

    rax = ppx(1) !r(iaxis)
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
    end subroutine find_new_axis

!-------------------------------------------------------------------
    subroutine compound_psi   ! to think about ferromags...

    psirz = psiplasrz + psiextrz

!    if (nferromag >= 1) then
!        psiferro = 0.
!        call ferro_mag_create
!        psirz = psirz + psiferro
!    endif

    return
    end subroutine compound_psi

!--------------------------------------------------------------------
    subroutine find_psi_boundary

    use pi_vars, only: GPI
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
                        if ((pos_xpoint(1) > Rrect(nr2) - dr) .or. (pos_xpoint(1) < Rrect(1) + dr) .or.  &
                            (pos_xpoint(2) > Zrect(nz2) - dz) .or. (pos_xpoint(2) < Zrect(1) + dz)) then
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
                    if ((pos_xpoint(1) > Rrect(nr2) - dr) .or. (pos_xpoint(1) < Rrect(1) + dr) .or.  &
                        (pos_xpoint(2) > Zrect(nz2) - dz) .or. (pos_xpoint(2) < Zrect(1) + dz)) then
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
                call nine_point_regression(Rrect(i), Zrect(j), pos_xpoint, ddipsi, x1)
                x5 = (ddipsi(5)**2 - ddipsi(3)*ddipsi(4))
                if ((pos_xpoint(1) >= Rrect(i) - dr) .and.  &
                    (pos_xpoint(1) <= Rrect(i) + dr) .and.  &
                    (pos_xpoint(2) >= Zrect(j) - dz) .and.  &
                    (pos_xpoint(2) <= Zrect(j) + dz) .and.  &
                    (x5 >= 0.)) then

                    n_of_xpoints = min(max_xpoints, n_of_xpoints + 1)
                    r_xpoint(n_of_xpoints) = pos_xpoint(1)
                    z_xpoint(n_of_xpoints) = pos_xpoint(2)
                endif
            enddo
        enddo

        do i=2, nr1, nr1-2
            do j=2, nz1
                call nine_point_regression(Rrect(i), Zrect(j), pos_xpoint, ddipsi, x1)
                x5 = (ddipsi(5)**2 - ddipsi(3)*ddipsi(4))
                if ((pos_xpoint(1) >= Rrect(i) - dr) .and.  &
                    (pos_xpoint(1) <= Rrect(i) + dr) .and.  &
                    (pos_xpoint(2) >= Zrect(j) - dz) .and.  &
                    (pos_xpoint(2) <= Zrect(j) + dz) .and.  &
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
                call nine_point_regression(Rrect(i), Zrect(j), pos_xpoint, ddipsi, x1)
                x5 = (ddipsi(5)**2 - ddipsi(3)*ddipsi(4))

                if ((pos_xpoint(1) >= Rrect(i) - dr) .and.  &
                    (pos_xpoint(1) <= Rrect(i) + dr) .and.  &
                    (pos_xpoint(2) >= Zrect(j) - dz) .and.  &
                    (pos_xpoint(2) <= Zrect(j) + dz) .and.  &
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
    if (use_limiter == 0) then
        limiterR = Rrect(nr1)
        limiterZ = Zrect(nz1)
    endif

! Calculate limiter flux
    do i=1, nlimiter
        psi_limp(i) = interp2d_psi(limiterR(i), limiterZ(i), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
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
            j = closest_index(r_xpoint(i), Rrect(1), dr)
            k = closest_index(z_xpoint(i), Zrect(1), dz)
            if (zlimpotential(j, k) == 0) then
                psi_xpoint(i) = -1.e6
            else
                psi_xpoint(i) = interp2d_psi(r_xpoint(i), z_xpoint(i), Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
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
            i_plasmatype = 0  ! limiter
        else
            i_plasmatype = 1  ! X-point
        endif
        active_x_point = i4
    endif

!    plasma_config = i_plasmatype

! if (i_plasmatype == 1)
    psibnd = psiaxis + (psibnd - psiaxis)*alpsep

! Normalized flux
    u_n(1:nr2, 1:nz2) = (psirz(1:nr2, 1:nz2) - psiaxis)/(psibnd - psiaxis)

    return
    end subroutine find_psi_boundary

!--------------------------------------------------------------------
    subroutine new_jrz ! calculate new right hand side given new boundary!

    use feqis_tools, only: fill_in_current, floor_index
    use global_params, only: iplasma

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

        i1 = floor_index(rax, Rrect(1), dr)
        j1 = floor_index(zax, Zrect(1), dz)

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

            dumc(i, j) = fill_in_current(Rrect(i), nrho2d, ppp_2d, ffp_2d, u_n(i, j))
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
        z11 = Rrect(i1-1)*(1 - t1) + Rrect(i1)*t1
        z12 = Rrect(i1)
        z13 = Rrect(i1 + 1)*(1 - t3) + Rrect(i1)*t3
        z14 = Rrect(i1)
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

    return
    end subroutine new_jrz

!---------------------------------------------------------------------
    double precision function t_find_u_n(i1, j1, i2, j2)

    integer, intent(in) :: i1, i2, j1, j2

    t_find_u_n = (1. - u_n(i1, j1))/(u_n(i2, j2) - u_n(i1, j1))

    return
    end function t_find_u_n

!--------------------------------------------------------------------
    subroutine solve_fbe_instantaneous(j_init, j_stab, raxold, zaxold)

    use feqis_tools, only: closest_index

    integer, intent(in) :: j_init, j_stab
    double precision, intent(in) :: raxold, zaxold

    integer :: i, j
    double precision :: curr, dum1, dum2, zum1, zum2, delr, delz
    double precision, dimension(9) :: c
    double precision, dimension(nr2, nz2) :: g

    g = 0.
    call solve_gs2d(g) !jrz as right hand side
    g = boundary(g)   ! gbound = integral (Green*dg/dn) over the boundary
    call solve_gs2d(g) ! again jrz as right hand side
    psiplasrz(1:nr2, 1:nz2) = g(1:nr2, 1:nz2)

    if (j_stab == 1) then
        psistabr = 0.
        psistabz = 0.
        delr = 0.
        delz = 0.
        call compound_psi
        call find_new_axis
        call nine_point_coeffs_only(raxold, zaxold, c, zum1, zum2)

! dpsidr
        zum1 = (raxold - zum1)/dr
        zum2 = (zaxold - zum2)/dz
        dum1 = 2.*c(2)*(zum1*zum2**2 + 2.*c(2)*zum1*zum2) +  &
            c(3)*zum2**2 + c(4)*zum2 + 2*c(5)*zum1 + c(7)

! dpsidz
       dum2 = 2.*c(1)*zum1**2*zum2 + c(2)*zum1**2 +  &
           2.*c(3)*zum1*zum2    + c(4)*zum1 + 2.*c(6)*zum2 + c(8)

        psistabr = -1./(2.*raxold)*dum1/dr
        psistabz = -dum2/dz

        call compound_psi
        do i=1, nr2
            do j=1, nz2
                psirz(i, j) = psirz(i, j) + psistabr*Rrect(i)**2 + psistabz*Zrect(j) ! Total flux
            enddo
        enddo
        call find_new_axis
    else
        call compound_psi ! Total flux
    endif

    call find_new_axis
    call find_psi_boundary
    call new_jrz

    return
    end subroutine solve_fbe_instantaneous

!--------------------------------------------------------------------
    subroutine solve_fbe_static_iterations_curgiven(j_init, raxp, zaxp, n_of_newton_iterations)

    use feqis_tools, only: closest_index
    use errors_params, only: err_find_psistab

    double precision, intent(in):: raxp, zaxp
    integer, intent(in):: j_init
    integer, intent(in):: n_of_newton_iterations

    integer :: j_iter, j_iter2, j_cyclo, jeppa
    double precision :: temp_err, raxold, zaxold, temp_err2, raxoldo, zaxoldo, &
        raxtmp, zaxtmp, det, psistab1o, psistab2o, psro, pszo, dist1, dist2, &
        cibapr, cibazr, rleft, rright, zup, zdown, dcrdr, dcrdz, dczdr, dczdz

! Start iterations to find self-consistent solution
    iaxis = closest_index(raxp, Rrect(1), dr)
    jaxis = closest_index(zaxp, Zrect(1), dz)
    rax = Rrect(iaxis)
    zax = Zrect(jaxis)
    raxold = rax
    zaxold = zax
    raxoldo = rax
    zaxoldo = zax
    psistabR = 0.
    psistabZ = 0.
    psistab1o = 0.
    psistab2o = 0.
    temp_err = 100.
    temp_err2 = 100.
    rleft = raxold
    rright = raxold
    zup = zaxold
    zdown = zaxold
! outer cycle, calculate new axis
    j_cyclo = 0
    jeppa = 1
    dist1 = dr*dr_factor_init
    dist2 = dz*dz_factor_init

    do j_iter2=1, 10000000

        SELECT CASE(j_cyclo)
        CASE(1)
            raxold = raxoldo + dist1
            zaxold = zaxoldo
            psistab1o = psistabr
            psistab2o = psistabz

        CASE(2)
            raxold = raxoldo
            zaxold = zaxoldo + dist2
            dcrdr = (psistabr - psistab1o)/dist1
            dczdr = (psistabz - psistab2o)/dist1

        CASE(3) ! inverse of (dcrdr dcrdz  ; dczdr dczdz ) = ( dczdz -dcrdz ; -dczdr dcrdr)
            dcrdz = (psistabr - psistab1o)/dist2
            dczdz = (psistabz - psistab2o)/dist2
            det = (dcrdr*dczdz - dcrdz*dczdr)
            cibapr = ( dczdz*psistab1o - dcrdz*psistab2o)/det
            cibazr = (-dczdr*psistab1o + dcrdr*psistab2o)/det
            raxold = raxoldo - cibapr
            zaxold = zaxoldo - cibazr
            raxoldo = raxold
            zaxoldo = zaxold
        END SELECT

! Inner cycle, calculate psi1 and psi2

        psro = 1000.
        pszo = 1000.
        redo_bnd = 1

        do j_iter=1, 10000
            raxtmp = trax
            zaxtmp = tzax
            call solve_fbe_instantaneous(j_iter - 1 + j_iter2 - 1, 1, raxold, zaxold)
            temp_err = (abs(psro - psistabr) + abs(pszo - psistabz))
            if (temp_err <= err_find_psistab) EXIT
            psro = psistabr
            pszo = psistabz
        enddo

        if (j_iter >= 10000) then
            write(*, *) 'total iterations passed, stopping!'
            stop
        endif

        if (j_cyclo == 3) then
            temp_err2 = abs(psistabr) + abs(psistabZ)
            j_cyclo = 1
            if (nint((0. + j_iter2)/3.) >= n_of_newton_iterations) EXIT
        else
            j_cyclo = j_cyclo + 1
        endif

        if (j_iter2 >= 400000) EXIT
        if (temp_err2 <= err_find_psistab) EXIT

    enddo

    return
    end subroutine solve_fbe_static_iterations_curgiven

end module fbe_core
