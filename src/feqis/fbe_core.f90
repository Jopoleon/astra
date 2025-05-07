module fbe_core

implicit none

integer, parameter :: nrho2d=251
integer :: use_limiter

integer :: nlimiter
double precision, dimension(:), allocatable :: limiterR, limiterZ
double precision :: lim_maxR, lim_minR, lim_maxZ, lim_minZ

!grids
integer :: nr, nz, nr2, nz2, nr1, nz1, nbnd
integer, dimension(:, :), allocatable :: zlimpotential

double precision :: dr, dz, zbot, ztop, raus, rinner
double precision, dimension(:), allocatable :: Rrect, Zrect
double precision, dimension(:, :), allocatable :: psi_n, psirz, psiextrz, psiplasrz, psiferro
double precision, dimension(nrho2d) :: ffp_2d, ppp_2d, psia_2d

!conductors
integer :: nconduc
double precision, dimension(:), allocatable :: curconduc

! boundary and axis FBE, PBE
integer, parameter :: max_xpoints=500
integer :: i_plasmatype !(0-limited, 1-single null, 2-double null)
integer :: iaxis, jaxis, n_of_xpoints, active_x_point
double precision :: psibnd, psiaxis, rax, zax, &
    alpsep, psistabR, psistabZ, dr_factor_init, dz_factor_init
double precision, dimension(max_xpoints) :: r_xpoint, z_xpoint, psi_xpoint
double precision, dimension(:), allocatable :: green_bnd_f

! plasma parameters
double precision, dimension(:, :), allocatable :: jrz
double precision, dimension(:), allocatable :: rbnd, zbnd

! reshape does transpose by default, order=(/2, 1/) restores the natural ordering

double precision, dimension(9, 9), parameter :: A_inv = reshape( (/ &
     0.25, -0.5,  0.25, -0.5,  1., -0.5,  0.25, -0.5, 0.25, &
    -0.25,  0.5, -0.25,  0. ,  0.,  0. ,  0.25, -0.5, 0.25, &
    -0.25,  0. ,  0.25,  0.5,  0., -0.5, -0.25,  0. , 0.25, &
     0.25,  0. , -0.25,  0. ,  0.,  0. , -0.25,  0. , 0.25, &
     0.  ,  0. ,  0.  ,  0.5, -1.,  0.5,  0.  ,  0. , 0.  , &
     0.  ,  0.5,  0.  ,  0. , -1.,  0. ,  0.  ,  0.5, 0.  , &
     0.  ,  0. ,  0.  , -0.5,  0.,  0.5,  0.  ,  0. , 0.  , &
     0.  , -0.5,  0.  ,  0. ,  0.,  0. ,  0.  ,  0.5, 0.  , &
     0.  ,  0. ,  0.  ,  0. ,  1.,  0. ,  0.  ,  0. , 0.  /), (/9, 9/), order=(/2, 1/) )

contains

!--------------------------------------------------------------------
    function solve_gs2d(greenBnd_in) result(green_out)

    use pi_vars, only: mu0
    use fft_mod_eff, only: costable
    use feqis_tools, only: discrete_sine_transform, solve_tridiag_fbe

    double precision, intent(in), dimension(2*nr+2*nz) :: greenBnd_in
    double precision, dimension(nr2, nz2) :: green_out

    integer :: i, j, k, j_init
    double precision :: x1, x2, r1m_1, r2m_1
    double precision, dimension(500) :: A, B, C, z_fourier
    double precision, dimension(nr2, nz2) :: gt1, gt2, rhs, wrhs

    data j_init/0/
    save A, B, C, j_init, z_fourier

    do j=1, nz2
        do i=1, nr2
            rhs(i, j) = -mu0*Rrect(i)*jrz(i, j)
        enddo
    enddo

    r1m_1 = Rrect(  2)/dr**2/((Rrect(  1) + Rrect(  2))/2.)
    r2m_1 = Rrect(nr1)/dr**2/((Rrect(nr1) + Rrect(nr2))/2.)

    green_out = 0.
    green_out(2:nr1,   1) = greenBnd_in(1:nr)
    green_out(nr2, 2:nz1) = greenBnd_in(nr+1: nr+nz)
    green_out(2:nr1, nz2) = greenBnd_in(nr+nz+1: 2*nr+nz)
    green_out(  1, 2:nz1) = greenBnd_in(2*nr+nz+1: 2*nr+2*nz)
    rhs(2:nr1,   2) = rhs(2:nr1,   2) - green_out(2:nr1,   1)/dz**2
    rhs(2:nr1, nz1) = rhs(2:nr1, nz1) - green_out(2:nr1, nz2)/dz**2
    rhs(  2, 2:nz1) = rhs(  2, 2:nz1) - green_out(  1, 2:nz1)*r1m_1
    rhs(nr1, 2:nz1) = rhs(nr1, 2:nz1) - green_out(nr2, 2:nz1)*r2m_1

    do i=2, nr1
        wrhs(i, 2:nz1) = discrete_sine_transform(nz, rhs(i, 2:nz1))
    enddo

! Create inverse matrix for gs2d
    if (j_init == 0) then
        A = 0.
        B = 0.
        C = 0.
        do i=1, nr
            x1 = 0.5*(Rrect(i+1) + Rrect(i  ))
            x2 = 0.5*(Rrect(i+2) + Rrect(i+1))
            B(i) = -Rrect(i+1)/dr**2*(1./x2 + 1./x1)
            A(i) =  Rrect(i+1)/x2/dr**2
            C(i) =  Rrect(i+1)/x1/dr**2
        enddo
        j_init = 1
        do k=2, nz1
            z_fourier(k) = 2./dz**2 * (costable(k-1) - 1.)
        enddo
    endif
     
    gt1 = 0.
! Solve matrix
    do k=2, nz1
        gt1(2:nr1, k) = solve_tridiag_fbe(C(1:nr), B(1:nr) + z_fourier(k), A(1:nr), wrhs(2:nr1, k), nr)
    enddo

! Invert fourier from gt(1:nr, 1:kfourier) to g(2:nr1, 2:nz1)
    do i=2, nr1
        gt2(i, 2:nz1) = discrete_sine_transform(nz, gt1(i, 2:nz1))
    enddo
    green_out(2:nr1, 2:nz1) = 2./(nz + 1)*gt2(2:nr1, 2:nz1)

    return
    end function solve_gs2d

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
    function getCoeffs(i_in, j_in) result(coeff)

    integer, intent(in) :: i_in, j_in

    integer :: i, j, k
    double precision, dimension(9) :: psi9, coeff

    k = 0
    do j=-1, 1
        do i=-1, 1
            k = k + 1
            psi9(k) = psirz(i_in+i, j_in+j)
        enddo
    enddo
    do k=1, 9
        coeff(k) = sum(A_inv(k, :) * psi9)
    enddo
      
    return
    end function getCoeffs

!---------------------------------------------------------------------
    subroutine nine_point_regression(i_in, j_in, r_out, z_out, psi_loc, hessian)

    use feqis_tools, only: exact_biquad

    integer, intent(in) :: i_in, j_in
    double precision, intent(out) :: r_out, z_out, psi_loc, hessian

    double precision :: dr_out, dz_out
    double precision, dimension(9) :: coeff

    coeff = getCoeffs(i_in, j_in)
    call exact_biquad(coeff, dr_out, dz_out, psi_loc, hessian, dr, dz)
    r_out = Rrect(i_in) + dr_out
    z_out = Zrect(j_in) + dz_out
    
    return
    end subroutine nine_point_regression

!---------------------------------------------------------------------
    subroutine nine_point_regression_follow(r_in, z_in, r_out, z_out, dpsi, psi_loc)

    use feqis_tools, only: interp2d_psi, exact_biquad_regress

    double precision, intent(in) :: r_in, z_in
    double precision, intent(out) :: psi_loc
    double precision, intent(out) :: r_out, z_out
    double precision, intent(out), dimension(5) :: dpsi

    integer :: k
    double precision :: dr_out, dz_out
    double precision, dimension(9) :: psi9, coeff

    psi9(1) = interp2d_psi(r_in - dr, z_in - dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    psi9(2) = interp2d_psi(r_in     , z_in - dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    psi9(3) = interp2d_psi(r_in + dr, z_in - dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    psi9(4) = interp2d_psi(r_in - dr, z_in     , Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    psi9(5) = interp2d_psi(r_in     , z_in     , Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    psi9(6) = interp2d_psi(r_in + dr, z_in     , Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    psi9(7) = interp2d_psi(r_in - dr, z_in + dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    psi9(8) = interp2d_psi(r_in     , z_in + dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    psi9(9) = interp2d_psi(r_in + dr, z_in + dz, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
    
    do k=1, 9
        coeff(k) = sum(A_inv(k, :) * psi9)
    enddo
    call exact_biquad_regress(coeff, dr_out, dz_out, psi_loc, dpsi, dr, dz)
    r_out = r_in + dr_out
    z_out = z_in + dz_out

    return
    end subroutine nine_point_regression_follow

!---------------------------------------------------------------------
    subroutine find_closest_xpoints(rx, zx, ierr, n_add)

!this routine finds the x-points close to the plasma boundary,  irrespective of other x-points

    integer, intent(out) :: ierr, n_add
    double precision, intent(out), dimension(200) :: rx, zx

    integer :: i, jinc, nx
    integer, dimension(250) :: jcycl
    double precision :: bx0, bx1, x1, posxR, posxZ
    double precision, dimension(5) :: dpsi

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
        call nine_point_regression_follow(rbnd(i), zbnd(i), posxR, posxZ, dpsi, x1)
        if (.not. isnan(x1)) then
            if (posxR <= lim_maxR .and. posxR >= lim_minR .and. posxZ >= lim_minZ .and. posxZ <= lim_maxZ) then
                jinc = jinc + 1
                jcycl(jinc) = i
            endif
        endif
    enddo

    if (jinc == 0) then
        ierr = 1
        return
    endif

    call nine_point_regression_follow(rbnd(jcycl(1)), zbnd(jcycl(1)), rx(1), zx(1), dpsi, x1)
    n_add = 1

    if (jinc == 1) then
        return
    endif

    if (jinc >= 2) then
!remove double counts
        do i=2, jinc
            call nine_point_regression_follow(rbnd(jcycl(i)), zbnd(jcycl(i)), posxR, posxZ, dpsi, x1)
            bx0 = sqrt((rx(i-1) - posxR)**2 + (zx(i-1) - posxZ)**2)
            if (bx0 > bx1) then
                n_add = n_add + 1
                rx(n_add) = posxR
                zx(n_add) = posxZ
            endif
        enddo
    endif

    if (n_add == 0) ierr=1

    return
    end subroutine find_closest_xpoints

!-----------------------------------------------------------------------------------
    function boundary(green_in) result(green_bnd)
! New bc is integral_over_boundary of -Green * dg/dn * dl

    double precision, intent(in), dimension(nr2, nz2) :: green_in

    integer :: j, jcount_in
    double precision, dimension(2*nr+2*nz) :: green_bnd

    jcount_in = 0
! lower, right, upper, left
    do j=1, 2*nr+2*nz
        green_bnd(j) = bgint(green_in, jcount_in)
        jcount_in = jcount_in + 2*nr + 2*nz
    enddo

    return
    end function boundary

!-----------------------------------------------------------------------------------
    double precision function bgint(green_in, jcount_in)
! Integral_over_boundary of -Green * dg/dn * dl

    use pi_vars, only: GPI

    double precision, intent(in), dimension(nr2, nz2) :: green_in
    integer, intent(in)  :: jcount_in

    integer :: j, jcount
    double precision, dimension(2*nr+2*nz) :: dgdn

    jcount = jcount_in
    do j=2, nr1 ! bottom
        jcount = jcount + 1
        dgdn(jcount-jcount_in) = green_in(j, 2) * green_bnd_f(jcount) * dr/dz * 1./Rrect(j)
    enddo      ! right
    do j=2, nz1
        jcount = jcount + 1
        dgdn(jcount-jcount_in) = green_in(nr1, j) * green_bnd_f(jcount) * dz/dr * 2./(Rrect(nr2) + Rrect(nr1))
    enddo      ! top
    do j=2, nr1
        jcount = jcount + 1
        dgdn(jcount-jcount_in) = green_in(j, nz1) * green_bnd_f(jcount) * dr/dz * 1./Rrect(j)
    enddo      ! left
    do j=2, nz1
        jcount = jcount + 1
        dgdn(jcount-jcount_in) = green_in(2, j) * green_bnd_f(jcount) * dz/dr * 2./(Rrect(1) + Rrect(2))
    enddo

    bgint = sum(dgdn)/GPI

    return
    end function bgint

!---------------------------------------------------------------------
    integer function xpoint_axis_connection(rx, zx, r0_in, z0_in, dr, dz)

!this routine checks that going from axis to x point,  the directed gradient of psi never changes sign
!(otherwise it means the x point is not connected to the plasma

    use feqis_tools, only: interp2d_psi

    double precision,  intent(in) :: rx, zx, r0_in, z0_in, dr, dz

    integer :: nsteps, i
    double precision :: angl, cos_ang, sin_ang, norm, dgrid, dd, &
        t1, t2, t3, t4, z1, z2, z3, psiold

    angl = ATAN2(rx - r0_in, zx - z0_in)
    cos_ang = COS(angl)
    sin_ang = SIN(angl)
    norm = sqrt((rx - r0_in)**2 + (zx - z0_in)**2)
    dgrid = sqrt(dr**2 + dz**2)
    nsteps = nint(norm/dgrid)
    dd = norm/nsteps !perfect ratio

    xpoint_axis_connection = 1
    psiold = 0.
    do i=2, nsteps
        t1 = r0_in + dd*(i - 1)*cos_ang
        t2 = z0_in + dd*(i - 1)*sin_ang
        t3 = r0_in + dd*i*cos_ang
        t4 = z0_in + dd*i*sin_ang
        z1 = interp2d_psi(t1, t2, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        z2 = interp2d_psi(t3, t4, Rrect(1:nr), Zrect(1:nz), psirz(1:nr, 1:nz))
        z3 = (z2 - z1)*psiold
        psiold = z2 - z1
        if (z3 < 0) then
            xpoint_axis_connection = 0
            EXIT
        endif
        if (dd*i >= norm) then
            EXIT
        endif
    enddo

    return
    end function xpoint_axis_connection

!--------------------------------------------------------------------
    subroutine find_new_axis

    integer :: j, iax, jax, i_old, j_old, ijmax(2)
    double precision :: hessian

! 1) find new magnetic axis

    iax = iaxis
    jax = jaxis
    i_old = 10000
    j_old = 10000
    do j=1, 100000
        ijmax = MAXLOC(psirz(iax-1:iax+1, jax-1:jax+1))
        iax = iax - 2 + ijmax(1)
        jax = jax - 2 + ijmax(2)
        if ((iax == i_old) .and. (jax == j_old)) EXIT
        i_old = iax
        j_old = jax
    enddo
    iaxis = iax
    jaxis = jax
    call nine_point_regression(iax, jax, rax, zax, psiaxis, hessian)

    if (isnan(rax)) then
        write(*, *) 'rax is nan in fbe find axis'
        stop
    endif
    if (rax < 0.1 .or. rax > 100) then
        write(*, *) 'rax is out of boundaries in fbe find axis', rax, zax, psiaxis, hessian
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
    subroutine add_xpoint(i_in, j_in)

    integer, intent(in) :: i_in, j_in
    double precision :: pos_xpointR, pos_xpointZ, x1, hessian

    call nine_point_regression(i_in, j_in, pos_xpointR, pos_xpointZ, x1, hessian)

    if ((pos_xpointR >= Rrect(i_in) - dr) .and. (pos_xpointR <= Rrect(i_in) + dr) .and.  &
        (pos_xpointZ >= Zrect(j_in) - dz) .and. (pos_xpointZ <= Zrect(j_in) + dz) .and.  &
        (hessian <= 0.)) then
        if (abs(pos_xpointR - rax) > 2.*dr .or. abs(pos_xpointZ - zax) > 2.*dz) then
            n_of_xpoints = min(max_xpoints, n_of_xpoints + 1)
            r_xpoint(n_of_xpoints) = pos_xpointR
            z_xpoint(n_of_xpoints) = pos_xpointZ
        endif
    endif

    return
    end subroutine add_xpoint

!--------------------------------------------------------------------
    subroutine find_psi_boundary

    use pi_vars, only: GPI
    use errors_params, only: err_find_oxpoints_derivs
    use feqis_tools, only: closest_index, pol_angle, interp2d_psi

    logical :: from_scratch
    integer :: niter, i, j, k, i1, i4, i5, i9, n_adding
    double precision :: x1, x2, x5, pos_xpointR, pos_xpointZ
    double precision, dimension(5) :: dpsi
    double precision, dimension(200) :: rx_add, zx_add
    double precision, dimension(500) :: psi_limp

    data from_scratch/.TRUE./
    save from_scratch

    i_plasmatype = 0

    if (n_of_xpoints == 0) from_scratch = .TRUE. !reset to full search if there are no x points!

    call find_closest_xpoints(rx_add, zx_add, i9, n_adding)
    if (i9 == 0) then
        r_xpoint(n_of_xpoints + 1:n_of_xpoints + n_adding) = rx_add(1:n_adding)
        z_xpoint(n_of_xpoints + 1:n_of_xpoints + n_adding) = zx_add(1:n_adding)
        n_of_xpoints = n_of_xpoints + n_adding
    endif

    if (from_scratch) then ! do a full pass to find all X-points
        n_of_xpoints = 0
        do j=2, nz1
            do i=2, nr1
                call add_xpoint(i, j)
            enddo
        enddo
    else
! go through old x-points and see where they end up
        if (n_of_xpoints >= 1) then
           do i=1, n_of_xpoints           
                do niter=1, 101
                    call nine_point_regression_follow(r_xpoint(i), z_xpoint(i), pos_xpointR, pos_xpointZ, dpsi, x1)
                    r_xpoint(i) = pos_xpointR
                    z_xpoint(i) = pos_xpointZ
                    if ( (abs(dpsi(1)) + abs(dpsi(2))) <= err_find_oxpoints_derivs) then
! Check if point outside of domain
                        if ((pos_xpointR > Rrect(nr2) - dr) .or. (pos_xpointR < Rrect(1) + dr) .or.  &
                            (pos_xpointZ > Zrect(nz2) - dz) .or. (pos_xpointZ < Zrect(1) + dz)) then
! xpoint doesnt exist anymore
                            r_xpoint(i) = 1.e6
                            z_xpoint(i) = 0.
                        else if ((pos_xpointR > rax - dr) .and. (pos_xpointR < rax + dr) .and.  &
                                 (pos_xpointZ > zax - dz) .and. (pos_xpointZ < zax + dz)) then
! xpoint doesn't exist anymore
                            r_xpoint(i) = 1.e6
                            z_xpoint(i) = 0.
                        else
                            r_xpoint(i) = pos_xpointR
                            z_xpoint(i) = pos_xpointZ
                        endif
                        EXIT
                    endif
                    if (niter > 100) then
                        r_xpoint(i) = 1.e6
                        z_xpoint(i) = 0.
                        EXIT
                    endif
                    if ((pos_xpointR > Rrect(nr2) - dr) .or. (pos_xpointR < Rrect(1) + dr) .or.  &
                        (pos_xpointZ > Zrect(nz2) - dz) .or. (pos_xpointZ < Zrect(1) + dz)) then
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
                call add_xpoint(i, j)
            enddo
        enddo
        do i=2, nr1, nr1-2
            do j=2, nz1
                call add_xpoint(i, j)
            enddo
        enddo
    endif
    from_scratch = .FALSE.

! Remove disappeared x-points and double counts

    i = 0
    xpoints_loop: do
        i = i + 1
        if (i > n_of_xpoints) EXIT xpoints_loop
        if (r_xpoint(i) >= 1.e5) then
            n_of_xpoints = n_of_xpoints - 1
            if (n_of_xpoints == 0) then
                EXIT xpoints_loop
            endif
            r_xpoint(i:n_of_xpoints) = r_xpoint(i+1:n_of_xpoints+1)
            z_xpoint(i:n_of_xpoints) = z_xpoint(i+1:n_of_xpoints+1)
            i = i - 1
            CYCLE xpoints_loop
        endif

        do k=1, i-1 ! check if double counted
            if ((abs(r_xpoint(i) - r_xpoint(k)) <= 2.*dr) .and.  &
                (abs(z_xpoint(i) - z_xpoint(k)) <= 2.*dz)) then

                r_xpoint(k) = 0.5*(r_xpoint(i) + r_xpoint(k))
                z_xpoint(k) = 0.5*(z_xpoint(i) + z_xpoint(k))
                n_of_xpoints = n_of_xpoints - 1
                r_xpoint(i:n_of_xpoints) = r_xpoint(i+1:n_of_xpoints+1)
                z_xpoint(i:n_of_xpoints) = z_xpoint(i+1:n_of_xpoints+1)
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
            if (limiterr(j) >= raus  ) psi_limp(j) = -1.e6
            if (limiterz(j) <= zbot  ) psi_limp(j) = -1.e6
            if (limiterz(j) >= ztop  ) psi_limp(j) = -1.e6
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
    psi_n(1:nr2, 1:nz2) = (psirz(1:nr2, 1:nz2) - psiaxis)/(psibnd - psiaxis)

    return
    end subroutine find_psi_boundary

!--------------------------------------------------------------------
    subroutine new_jrz ! calculate new right hand side given new boundary!

    use feqis_tools, only: fill_in_current, floor_index
    use global_params, only: iplasma

    integer :: i, j, i1, i2, j1, quadrant, ipluz, jpluz, &
        ilast, totpoints, istart, j_griddo_j
    integer, dimension(nr2*nz2, 2) :: external_griddo_j
    double precision :: curr, t1, t2, t3, t4, je1, je2, je3, je4, &
        z11, z12, z13
    double precision, dimension(nr2, nz2) :: iconvex
    double precision, dimension(nr2, nz2) :: dumc

! in entry: rbnd, zbnd, nbnd, psiaxis, psibnd, psi_n

    iconvex = 0.
    dumc    = 0.
    j_griddo_j = 0

    i1 = floor_index(rax, Rrect(1), dr)
    j1 = floor_index(zax, Zrect(1), dz)

    totpoints = 0
    do quadrant=1, 4
! sweep from axis to exterior and fill in the current
!start from axis position

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

            dumc(i, j) = fill_in_current(Rrect(i), nrho2d, ppp_2d, ffp_2d, psi_n(i, j))
            iconvex(i, j) = 1.
            t1 = (1. - psi_n(i + ipluz, j))/(psi_n(i, j) - psi_n(i + ipluz, j))
            t2 = (1. - psi_n(i, j + jpluz))/(psi_n(i, j) - psi_n(i, j + jpluz))
            if (t2 > 0. .and. t2 <= 1.) then
                j_griddo_j = j_griddo_j + 1
                external_griddo_j(j_griddo_j, 1) = i
                external_griddo_j(j_griddo_j, 2) = j + jpluz
            endif

            if (t1 < 0 .or. t1 > 1.) then
! Move horizontally to the right
                i = i + ipluz
            else
                ilast = i + ipluz
                j_griddo_j = j_griddo_j + 1
                external_griddo_j(j_griddo_j, 1) = ilast
                external_griddo_j(j_griddo_j, 2) = j
! Found boundary, go back, check vertically
                i = istart
                do
                    t2 = (1. - psi_n(i, j + jpluz))/(psi_n(i, j) - psi_n(i, j + jpluz))
                    if (t2 < 0 .or. t2 > 1.) then
                        j = j + jpluz
                        EXIT
                    endif

! Found boundary on Z, need to advance 1 more
                    j_griddo_j = j_griddo_j + 1
                    external_griddo_j(j_griddo_j, 1) = i
                    external_griddo_j(j_griddo_j, 2) = j + jpluz

                    i = i + ipluz
                    istart = i
                    if (i == ilast) EXIT inner_loop
                enddo
            endif
        enddo inner_loop
    enddo

! fill current in external griddo
    do i=1, j_griddo_j
        i1 = external_griddo_j(i, 1)
        i2 = external_griddo_j(i, 2)
        t1 = (1. - psi_n(i1-1, i2))/(psi_n(i1, i2) - psi_n(i1-1, i2))
        t2 = (1. - psi_n(i1, i2+1))/(psi_n(i1, i2) - psi_n(i1, i2+1))
        t3 = (1. - psi_n(i1+1, i2))/(psi_n(i1, i2) - psi_n(i1+1, i2))
        t4 = (1. - psi_n(i1, i2-1))/(psi_n(i1, i2) - psi_n(i1, i2-1))
        if (t1 < 0 .or. t1 > 1.) t1 = 0.
        if (t2 < 0 .or. t2 > 1.) t2 = 0.
        if (t3 < 0 .or. t3 > 1.) t3 = 0.
        if (t4 < 0 .or. t4 > 1.) t4 = 0.
        z11 = Rrect(i1-1) + t1*dr
        z12 = Rrect(i1)
        z13 = Rrect(i1+1) - t3*dr
        je1 = ppp_2d(nrho2d)*z11 + ffp_2d(nrho2d)/z11
        je2 = ppp_2d(nrho2d)*z12 + ffp_2d(nrho2d)/z12
        je3 = ppp_2d(nrho2d)*z13 + ffp_2d(nrho2d)/z13
        je4 = je2

! defining S1 = dR - d1, S2 = dZ - d2, S3 = dZ - d3, S4 = dR - d4
! Jvacuum = C*Jb
! if 1 point only: C = 1 - S1/dR
! if 2 points: C = 1 - S1*S2/(dR*dZ)
! if 3 points: C = 1 - (S1*S2 + S1*S3)/(2*dR*dZ)
! if 4 points: C = 1 - (S1*S2 + S1*S3 + S3*S4 + S4*S2)/(4*dR*dZ)
! t_j = d_j /(dR or dZ depending on direction)

        dumc(i1, i2) = t1*je1 + t2*je2 + t3*je3 + t4*je4 - 0.5* &
            (t1*t2*(je1 + je2) +  &
             t1*t3*(je1 + je3) +  &
             t1*t4*(je1 + je4) +  &
             t2*t3*(je3 + je2) +  &
             t2*t4*(je4 + je2) +  &
             t3*t4*(je3 + je4) )

        iconvex(i1, i2) = t1 + t2 + t3 + t4 - (t1*t2 + t1*t3 + t1*t4 + t2*t3 + t2*t4 + t3*t4)

    enddo

    jrz = 0.
    do j=2, nz1
        do i=2, nr1
            jrz(i, j) = iconvex(i, j)*0.5*(dumc(i, j) + 0.25*(dumc(i+1, j) + dumc(i-1, j) + dumc(i, j-1) + dumc(i, j+1))) ! this is a first order Shapiro filter with alpha = 0.5, with in addition the iconvex multiplier which reduces the ghost currents even more to avoid overshooting. Seems to work well vs spider
        enddo
    enddo

!uncomment below for consistent current
!jrz(1:nr2, 1:nz2)=dumc(1:nr2, 1:nz2)
    curr = sum(jrz)*dr*dz
    jrz = jrz/curr*iplasma

    if (isnan(curr)) then
        write(*, *) 'Total current is nan in fbe current rescaling'
        stop
    endif

    return
    end subroutine new_jrz

!--------------------------------------------------------------------
    function get_psiplasrz result(psi_plas)

    double precision, dimension(nr2, nz2) :: green, psi_plas
    double precision, dimension(2*nr + 2*nz) :: greenBnd

    greenBnd = 0.
    green = solve_gs2d(greenBnd)    ! jrz as right hand side
    greenBnd = boundary(green)      ! gbound = integral (Green*dg/dn) over the boundary
    psi_plas = solve_gs2d(greenBnd) ! again jrz as right hand side

    return
    end function get_psiplasrz

!--------------------------------------------------------------------
    subroutine solve_fbe_instantaneous(j_stab, raxold, zaxold)

    use feqis_tools, only: closest_index, expandCoeffs

    integer, intent(in) :: j_stab
    double precision, intent(in) :: raxold, zaxold

    integer :: i, j, iloc, jloc
    double precision :: u_loc, zum1, zum2
    double precision, dimension(5) :: dpsi
    double precision, dimension(9) :: coeff

    psiplasrz = get_psiplasrz()
    call compound_psi
    if (j_stab == 1) then
        iloc = closest_index(raxold, Rrect(1), dr)
        jloc = closest_index(zaxold, Zrect(1), dz)
        zum1 = (raxold - Rrect(iloc))/dr
        zum2 = (zaxold - Zrect(jloc))/dz
        coeff = getCoeffs(iloc, jloc)
        call expandCoeffs(coeff, zum1, zum2, u_loc, dpsi)
        psistabr = -dpsi(1)/(2.*raxold*dr)
        psistabz = -dpsi(2)/dz
        do i=1, nr2
            do j=1, nz2
                psirz(i, j) = psirz(i, j) + psistabr*Rrect(i)**2 + psistabz*Zrect(j)
            enddo
        enddo
    endif
    call find_new_axis
    call find_psi_boundary
    call new_jrz

    return
    end subroutine solve_fbe_instantaneous

!--------------------------------------------------------------------
    subroutine solve_fbe_static_iterations_curgiven(raxp, zaxp, n_of_newton_iterations)

    use feqis_tools, only: closest_index
    use errors_params, only: err_find_psistab

    double precision, intent(in):: raxp, zaxp
    integer, intent(in):: n_of_newton_iterations

    integer :: j_iter, j_iter2, j_cyclo
    double precision :: temp_err, raxold, zaxold, temp_err2, raxoldo, zaxoldo, &
        det, psistab1o, psistab2o, psro, pszo, dist1, dist2, &
        cibapr, cibazr, dcrdr, dcrdz, dczdr, dczdz

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
! outer cycle, calculate new axis
    j_cyclo = 0
    dist1 = dr*dr_factor_init
    dist2 = dz*dz_factor_init

    do j_iter2=1, 400000

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

        do j_iter=1, 10000
            call solve_fbe_instantaneous(1, raxold, zaxold)
            temp_err = abs(psro - psistabr) + abs(pszo - psistabz)
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
        if (temp_err2 <= err_find_psistab) EXIT
    enddo

    return
    end subroutine solve_fbe_static_iterations_curgiven

end module fbe_core
