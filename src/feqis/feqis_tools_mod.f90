module feqis_tools

implicit none

! reshape does transpose!

double precision, dimension(9, 9), parameter :: A_inv = reshape( (/ &
     0.25, -0.25, -0.25,  0.25,  0.00,  0.00,  0.00,  0.00,  0.00, &
    -0.50,  0.50,  0.00,  0.00,  0.00,  0.50,  0.00, -0.50,  0.00, &
     0.25, -0.25,  0.25, -0.25,  0.00,  0.00,  0.00,  0.00,  0.00, &
    -0.50,  0.00,  0.50,  0.00,  0.50,  0.00, -0.50,  0.00,  0.00, &
     1.00,  0.00,  0.00,  0.00, -1.00, -1.00,  0.00,  0.00,  1.00, &
    -0.50,  0.00, -0.50,  0.00,  0.50,  0.00,  0.50,  0.00,  0.00, &
     0.25,  0.25, -0.25, -0.25,  0.00,  0.00,  0.00,  0.00,  0.00, &
    -0.50, -0.50,  0.00,  0.00,  0.00,  0.50,  0.00,  0.50,  0.00, &
     0.25,  0.25,  0.25,  0.25,  0.00,  0.00,  0.00,  0.00,  0.00 /), (/9, 9/) )

contains

!---------------------------------------------------------------------
    double precision function pol_angle(r0, z0, r, z)

    use pi_vars, only: GPI2

    double precision, intent(in) :: r0, z0, r, z

    pol_angle = ATAN2(z - z0, r - r0)
    if (pol_angle < 0) pol_angle = pol_angle + GPI2

    return
    end function pol_angle

!---------------------------------------------------------------------
    function inv_matrix(a_in, ndim) result(ainv_out)

    integer, intent(in) :: ndim
    double precision, intent(in), dimension(ndim, ndim) ::  a_in
    double precision, dimension(ndim, ndim) :: ainv_out

    integer :: i, j, k
    double precision :: coeff
    double precision, dimension(ndim) :: b, d, x
    double precision, dimension(ndim, ndim) :: a, L, U

    L = 0.0
    U = 0.0
    b = 0.0
    a = a_in

! step 1: forward elimination
    do k=1, ndim-1
        do i=k+1, ndim
            coeff = a(i, k)/a(k, k)
            L(i, k) = coeff
            do j=k+1, ndim
                 a(i, j) = a(i, j) - coeff*a(k, j)
            enddo
        enddo
    enddo

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
    do i=1, ndim
        L(i, i) = 1.0
    enddo
! U matrix is the upper triangular part of A
    do j=1, ndim
        do i=1, j
            U(i, j) = a(i, j)
        enddo
    enddo

! Step 3: compute columns of the inverse matrix C
    do k=1, ndim
        b(k) = 1.0
        d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
        do i=2, ndim
            d(i) = b(i)
            do j=1, i-1
                d(i) = d(i) - L(i, j)*d(j)
            enddo
        enddo
! Step 3b: Solve Ux=d using the back substitution
        x(ndim) = d(ndim)/U(ndim, ndim)
        do i = ndim-1, 1, -1
            x(i) = d(i)
            do j=ndim, i+1, -1
                x(i) = x(i) - U(i, j)*x(j)
            enddo
            x(i) = x(i)/u(i, i)
        enddo
! Step 3c: fill the solutions x(n) into column k of ainv_out
        do i=1, ndim
            ainv_out(i, k) = x(i)
        enddo
        b(k) = 0.0
    enddo
    
    return
    end function inv_matrix

!---------------------------------------------------------------------
    double precision function curinterp(r_in, z_in, jrho, rho, teta, rax, zax, nrho, nteta)

    use pi_vars, only: GPI2
    use numerical_tools, only: linterp

    integer,  intent(in) :: nrho, nteta
    double precision, intent(in) :: r_in, z_in, rax, zax
    double precision, intent(in), dimension(nteta) :: teta
    double precision, intent(in), dimension(nrho, nteta) :: jrho, rho

    integer :: i, k, k1, k2, k3, k4, j1, j2
    double precision :: anglr, rho0, r1, r2, r3, r4, z1, z2, z3, z4, a1, a2, a3, a4, d1, d2, d3, d4
    double precision, dimension(2) :: x1, y1, jt
    double precision, dimension(4) :: jj1, coef
    double precision, dimension(4, 4) :: matrix, imatrix

    anglr = pol_angle(rax, zax, r_in, z_in)
    rho0 = sqrt((r_in - rax)**2 + (z_in - zax)**2)

    if (anglr < teta(1)) anglr = anglr + GPI2

    j1 = 1
    do i=1, nteta
        if (anglr >= teta(i)) j1 = i
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

    r1 = rax + rho(k1, j1)*cos(teta(j1))
    r2 = rax + rho(k2, j1)*cos(teta(j1))
    r3 = rax + rho(k4, j2)*cos(teta(j2))
    r4 = rax + rho(k3, j2)*cos(teta(j2))
    z1 = zax + rho(k1, j1)*sin(teta(j1))
    z2 = zax + rho(k2, j1)*sin(teta(j1))
    z3 = zax + rho(k4, j2)*sin(teta(j2))
    z4 = zax + rho(k3, j2)*sin(teta(j2))
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

    return
    end function curinterp

!---------------------------------------------------------------------
    function discrete_sine_transform(ndim, f_in) result(f_out)

    use fft_mod_eff, only: dp, sintable

    integer, intent(in) :: ndim
    double precision,  intent(in), dimension(ndim) :: f_in
    double precision, dimension(ndim) :: f_out

    integer :: i, j, imethod1, icall, k
    double precision :: z(ndim)
    complex(kind=dp) :: d(2*(ndim+1))

    imethod1 = 1
    f_out = f_in
    if (imethod1 == 1) then
        z = 0.
        do i=1, ndim
            do j=1, ndim
                z(i) = z(i) + f_out(j)*sintable(i, j)
            enddo
        enddo
        f_out = z
        return
    else if (imethod1 == 2) then 
! fast sine transform 
! this problem is equivalent to DST-I with N = n+1 

        k = 2*(ndim+1)
        z = f_out
        d(1) = cmplx(0., 0.)
        do i=1, ndim
            d(i+1)   = cmplx( z(i), 0)
            d(k-i+1) = cmplx(-z(i), 0)
        enddo
        d(k-ndim) = cmplx(0., 0.)

        call fft_eff(d)

        d(1: k-1) = d(2: k)
        do i=1, ndim
            z(i) = 0.5*aimag(d(i) - d(k-i))
        enddo
        f_out = z/2.
    endif

    return
    end function discrete_sine_transform

!---------------------------------------------------------------------
    subroutine get_zccurb_efff(rc_cur, zc_cur, z2c_cur, rgeoc, zgeoc, ahorc)

    use feqis_circuit, only: nrho, nteta, rpol, zpol
    use metric_coefficients_pbe, only: R_curr_0D, Z_curr_0D, dator

    real*8, intent(out) :: rc_cur, zc_cur, z2c_cur, rgeoc, zgeoc, ahorc

    integer :: i, j
    real*8 :: perimz, ahorc2, avgelem

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
    end subroutine get_zccurb_efff

!---------------------------------------------------------------------
    double precision function psib_ext_feqis

! gives back external flux on plasma boundary

    use feqis_circuit, only: nr, nz, nbnd, rbnd, zbnd, psibnd, psiextrz

    integer :: i
    double precision :: psiext_out, dlt, dllt, psi_ext_1, psiext1, psiext2

!cycle over boundary
    psiext_out = 0.
    dllt = 0.
    psi_ext_1 = interp2d_psi(rbnd(1), zbnd(1), psiextrz(1:nr, 1:nz))
    psiext1 = psi_ext_1
    do i=1, nbnd-1
        psiext2 = interp2d_psi(rbnd(i+1), zbnd(i+1), psiextrz(1:nr, 1:nz))
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
    use feqis_circuit, only: nr, nz, psirz

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
        u001 = interp2d_psi(x1, y1, psirz(1: nr, 1:nz))
        u002 = interp2d_psi(x2, y2, psirz(1: nr, 1:nz))
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

    use feqis_circuit, only: psibnd

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

    use feqis_circuit, only: nr2, nz2, nconduc, curconduc, psiextrz
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
    function least_square_biquad(r, z, u, n) result(derivs)

    integer, intent(in) :: n
    double precision, intent(in), dimension(n) :: r, z, u
    double precision, dimension(8) :: derivs
    
    integer :: k
    double precision :: det, det_r, det_z, rax, zax, uax
    double precision, dimension(6) :: B, cc
    double precision, dimension(21) :: sums
    double precision, dimension(6, 6) :: A, Ainv

    sums = 0.

    sums(1)  = sum(r**4)
    sums(2)  = sum(z**4)
    sums(3)  = sum(r**2 * z**2)
    sums(4)  = sum(r**2)
    sums(5)  = sum(z**2)
    sums(6)  = n + 0.
    sums(7)  = sum(r**3 * z)
    sums(8)  = sum(r**3)
    sums(9)  = sum(r**2 * z)
    sums(10) = sum(r * z**3)
    sums(11) = sum(r * z**2)
    sums(12) = sum(z**3)
    sums(13) = sum(r*z)
    sums(14) = sum(r)
    sums(15) = sum(z)
    sums(16) = sum(u * r**2)
    sums(17) = sum(u * z**2)
    sums(18) = sum(u*r*z)
    sums(19) = sum(u*r)
    sums(20) = sum(u*z)
    sums(21) = sum(u)

    B(1) = -2*sums(16)
    B(2) = -2*sums(17)
    B(3) = -2*sums(18)
    B(4) = -2*sums(19)
    B(5) = -2*sums(20)
    B(6) = -2*sums(21)

    A(1, 1) = 4*sums(1)
    A(1, 2) = 4*sums(3)
    A(1, 3) = 4*sums(7)
    A(1, 4) = 4*sums(8)
    A(1, 5) = 4*sums(9)
    A(1, 6) = 4*sums(4)

    A(2, 1) = A(1, 2)
    A(2, 2) = 4*sums(2)
    A(2, 3) = 4*sums(10)
    A(2, 4) = 4*sums(11)
    A(2, 5) = 4*sums(12)
    A(2, 6) = 4*sums(5)

    A(3, 1) = A(1, 3)
    A(3, 2) = A(2, 3)
    A(3, 3) = 4*sums(3)
    A(3, 4) = 4*sums(9)
    A(3, 5) = 4*sums(11)
    A(3, 6) = 4*sums(13)

    A(4, 1) = A(1, 4)
    A(4, 2) = A(2, 4)
    A(4, 3) = A(3, 4)
    A(4, 4) = 4*sums(4)
    A(4, 5) = 4*sums(13)
    A(4, 6) = 4*sums(14)

    A(5, 1) = A(1, 5)
    A(5, 2) = A(2, 5)
    A(5, 3) = A(3, 5)
    A(5, 4) = A(4, 5)
    A(5, 5) = 4*sums(5)
    A(5, 6) = 4*sums(15)

    A(6, 1) = A(1, 6)
    A(6, 2) = A(2, 6)
    A(6, 3) = A(3, 6)
    A(6, 4) = A(4, 6)
    A(6, 5) = A(5, 6)
    A(6, 6) = 4*sums(6)

!find coefficients
    Ainv  = inv_matrix(A, 6)

    do k=1, 6
        cc(k) = -2*sum(Ainv(k, 1: 6)*B(1: 6))
    enddo

!  magnetic axis

    det   =  4.d0*cc(1)*cc(2) - cc(3)**2
    det_r = -2.d0*cc(2)*cc(4) + cc(3)*cc(5)
    det_z = -2.d0*cc(1)*cc(5) + cc(3)*cc(4)

    rax = det_r/det
    zax = det_z/det

    uax = cc(1)*rax**2 + cc(2)*zax**2 + cc(3)*rax*zax + cc(4)*rax + cc(5)*zax + cc(6)

    derivs(1) = cc(4)
    derivs(2) = cc(5)
    derivs(3) = 2.*cc(1)
    derivs(4) = 2.*cc(2)
    derivs(5) = cc(3) 

    return
    end function least_square_biquad

!---------------------------------------------------------------------
    subroutine exact_biquad(rx_in, zx_in, u, ndim, rax, zax, uax, derivs, dr, dz)

    use errors_params, only: err_find_biquad

    integer, intent(in) :: ndim
    double precision, intent(in) :: dr, dz, rx_in, zx_in
    double precision, intent(in), dimension(ndim) :: u
    double precision, intent(out) :: rax, zax, uax
    double precision, intent(out), dimension(8) :: derivs

    integer :: k, niter, j_success
    double precision :: s_r, s_z, s_r2, s_z2, det
    double precision :: A(2, 2), B(2), c(9)

! Find coefficients
    do k=1, 9
        c(k) = sum(A_inv(k, 1:ndim) * u(1:ndim))
    enddo

    rax = 0.
    zax = 0.

!now find axis
    niter = 0
    s_r = 100000.
    s_z = 100000.

    do
        niter = niter + 1
        s_r2 = 2*C(1)*rax*zax**2 + 2*C(2)*rax*zax +   C(3)*zax**2  + C(4)*zax + 2*C(5)*rax + C(7)
        s_z2 = 2*C(1)*rax**2*zax +   C(2)*rax**2  + 2*C(3)*zax*rax + C(4)*rax + 2*C(6)*zax + C(8)
        A(1, 1) = 2*C(1)*zax**2  + 2*C(2)*zax + 2*C(5)
        A(1, 2) = 4*C(1)*rax*zax + 2*C(2)*rax + 2*C(3)*zax + C(4)
        A(2, 2) = 2*C(1)*rax**2  + 2*C(3)*rax + 2*c(6)
        A(2, 1) = 4*C(1)*rax*zax + 2*C(2)*rax + 2*C(3)*zax + C(4)
        B(1) = s_r2
        B(2) = s_z2
        det = (A(1, 1)*A(2, 2)) - (A(1, 2)*A(2, 1))
        s_r2 = 1/det*(A(2, 2)*B(1) - A(1, 2)*B(2))
        s_z2 = 1/det*(A(1, 1)*B(2) - A(2, 1)*B(1))
        rax = rax - s_r2
        zax = zax - s_z2

        s_r2 = s_r
        s_z2 = s_z
        s_r = 2*C(1)*rax*zax**2 + 2*C(2)*rax*zax +   C(3)*zax**2  + C(4)*zax + 2*C(5)*rax + C(7)
        s_z = 2*C(1)*rax**2*zax + C(2)*rax**2    + 2*C(3)*zax*rax + C(4)*rax + 2*C(6)*zax + C(8)

        if (abs(s_r) < err_find_biquad .and. abs(s_z) < err_find_biquad) then
            j_success = 1
            EXIT
        endif
        if (abs(rax) > 1. .or. abs(zax) > 1 .or. niter > 100000) then
            j_success = 0
            EXIT
        endif
    enddo

    if (j_success == 0) then
        rax    =  1.e6
        zax    =  1.e6
        uax    = -1.e6
        derivs =  1.e6
        return
    endif

! magnetic axis

    uax = c(1)*rax**2 * zax**2 + & 
          c(2)*rax**2 * zax    + &
          c(3)*rax    * zax**2 + &
          c(4)*rax    * zax    + &
          c(5)*rax**2 + &
          c(6)*zax**2 + &
          c(7)*rax    + &
          c(8)*zax    + &
          c(9)

    derivs(1) = 1/dr*(2*C(1)*rax*zax**2 + 2*C(2)*rax*zax +   C(3)*zax**2  + C(4)*zax + 2*C(5)*rax + C(7))
    derivs(2) = 1/dz*(2*C(1)*rax**2*zax +   C(2)*rax**2  + 2*C(3)*zax*rax + C(4)*rax + 2*C(6)*zax + C(8))
    derivs(3) = 1/dr**2*(2.*c(1)*zax**2 + 2*c(2)*zax + 2*c(5))
    derivs(4) = 1/dz**2*(2.*c(1)*rax**2 + 2*c(3)*rax + 2*c(6))
    derivs(5) = 1/dr/dz*(4*c(1)*rax*zax + 2*c(2)*rax + 2*c(3)*zax + c(4)) 

    rax = rax*dr + rx_in
    zax = zax*dz + zx_in

    return
    end subroutine exact_biquad

!---------------------------------------------------------------------
    subroutine exact_biquad_regress(rx_in, zx_in, u, ndim, rax, zax, uax, derivs, dr, dz)

    integer, intent(in) :: ndim
    double precision, intent(in) :: dr, dz, rx_in, zx_in
    double precision, intent(in), dimension(ndim) :: u
    double precision, intent(out) :: rax, zax, uax
    double precision, intent(out), dimension(ndim-1) :: derivs

    integer :: k
    double precision :: s_r2, s_z2, det
    double precision :: A(2, 2), B(2)
    double precision :: c(9)

! Transform
    rax = 0.
    zax = 0.

! Find coefficients
    do k=1, ndim
        c(k) = sum(A_inv(k, 1:ndim)*u(1:ndim))
    enddo

! Now find axis
    s_r2 = 2.*C(1)*rax*zax**2 + 2.*C(2)*rax*zax +    C(3)*zax**2  + C(4)*zax + 2.*C(5)*rax + C(7)
    s_z2 = 2.*C(1)*rax**2*zax +    C(2)*rax**2  + 2.*C(3)*zax*rax + C(4)*rax + 2.*C(6)*zax + C(8)
    A(1, 1) = 2.*C(1)*zax**2  + 2.*C(2)*zax + 2.*C(5)
    A(1, 2) = 4.*C(1)*rax*zax + 2.*C(2)*rax + 2.*C(3)*zax + C(4)
    A(2, 2) = 2.*C(1)*rax**2  + 2.*C(3)*rax + 2.*c(6)
    A(2, 1) = A(1, 2)
    B(1) = s_r2
    B(2) = s_z2
    det = (A(1, 1)*A(2, 2)) - (A(1, 2)*A(2, 1))
    s_r2 = 1./det*(A(2, 2)*B(1) - A(1, 2)*B(2))
    s_z2 = 1./det*(A(1, 1)*B(2) - A(2, 1)*B(1))
    rax = rax - s_r2
    zax = zax - s_z2

    uax = c(1)*rax**2 * zax**2 + & 
          c(2)*rax**2 * zax + &
          c(3)*rax * zax**2 + &
          c(4)*rax * zax + &
          c(5)*rax**2 + &
          c(6)*zax**2 + &
          c(7)*rax    + &
          c(8)*zax + &
          c(9)

    derivs(1) = 1./dr*(2.*C(1)*rax*zax**2 + 2.*C(2)*rax*zax +    C(3)*zax**2  + C(4)*zax + 2.*C(5)*rax + C(7))
    derivs(2) = 1./dz*(2.*C(1)*rax**2*zax +    C(2)*rax**2  + 2.*C(3)*zax*rax + C(4)*rax + 2.*C(6)*zax + C(8))
    derivs(3) = 1./dr**2*(2.*c(1)*zax**2  + 2.*c(2)*zax + 2.*c(5))
    derivs(4) = 1./dz**2*(2.*c(1)*rax**2  + 2.*c(3)*rax + 2.*c(6))
    derivs(5) = 1./dr/dz*(4.*c(1)*rax*zax + 2.*c(2)*rax + 2.*c(3)*zax + c(4)) 
     
    rax = rax*dr + rx_in
    zax = zax*dz + zx_in

    return
    end subroutine exact_biquad_regress

!---------------------------------------------------------------------
    subroutine nine_point_regression(r_in, z_in, pos_xpoint, ddpsi, f00)

    use feqis_circuit, only: nr1, nz1, r, z, dr, dz, psirz

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

    use feqis_circuit, only: psirz, nr1, nz1, r, z, dr, dz, psirz
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

    use feqis_circuit, only: nr, nz, dr, dz, psirz

    integer, parameter :: ndim=9
    double precision, intent(in) :: rx, zx
    double precision, intent(out) :: f00
    double precision, intent(out), dimension(2) :: pos_xpoint
    double precision, intent(out), dimension(ndim-1) :: ddpsi
    
    double precision, dimension(ndim) :: bub

    bub(1) = interp2d_psi(rx - dr, zx - dz, psirz(1:nr, 1:nz))
    bub(2) = interp2d_psi(rx     , zx - dz, psirz(1:nr, 1:nz))
    bub(3) = interp2d_psi(rx + dr, zx - dz, psirz(1:nr, 1:nz))
    bub(4) = interp2d_psi(rx - dr, zx     , psirz(1:nr, 1:nz))
    bub(5) = interp2d_psi(rx     , zx     , psirz(1:nr, 1:nz))
    bub(6) = interp2d_psi(rx + dr, zx     , psirz(1:nr, 1:nz))
    bub(7) = interp2d_psi(rx - dr, zx + dz, psirz(1:nr, 1:nz))
    bub(8) = interp2d_psi(rx     , zx + dz, psirz(1:nr, 1:nz))
    bub(9) = interp2d_psi(rx + dr, zx + dz, psirz(1:nr, 1:nz))

    call exact_biquad_regress(rx, zx, bub(1:ndim), ndim,  &
        pos_xpoint(1), pos_xpoint(2), f00, ddpsi, dr, dz)

    return
    end subroutine nine_point_regression_follow

!---------------------------------------------------------------------
    integer function closest_index(x_in, xmin, dx)

    double precision, intent(in) :: x_in, xmin, dx

    closest_index = nint((x_in - xmin)/dx + 1.) ! nint(1.8) = 2

    return
    end function closest_index

!---------------------------------------------------------------------
    integer function floor_index(x_in, xmin, dx)

    double precision, intent(in) :: x_in, xmin, dx

    floor_index = floor((x_in - xmin)/dx + 1.) ! floor(1.8) = 1

    return
    end function floor_index

!---------------------------------------------------------------------
    double precision function interp2d_psi(r_in, z_in, psi_in)
! ->psi at r_in, z_in
    use feqis_circuit, only: r, z, dr, dz

    double precision, intent(in) :: r_in, z_in
    double precision, intent(in), dimension(:, :) :: psi_in

    integer :: i, j
    integer, dimension(2) :: psi_shape
    double precision :: r1, r2, z1, z2, psi1, psi2, psi3, psi4

    psi_shape = SHAPE(psi_in)

    i = floor_index(r_in, r(1), dr)
    j = floor_index(z_in, z(1), dz)
    i = min(psi_shape(1) - 1, max(1, i))
    j = min(psi_shape(2) - 1, max(1, j))

    r1 = r(i)
    r2 = r(i+1)
    z1 = z(j)
    z2 = z(j+1)
    psi1 = psi_in(i, j)
    psi2 = psi_in(i+1, j)
    psi3 = psi_in(i, j+1)
    psi4 = psi_in(i+1, j+1)

    interp2d_psi = bilinear_interp(r1, r2, z1, z2, r_in, z_in, psi1, psi2, psi3, psi4)

    return
    end function interp2d_psi

!---------------------------------------------------------------------
    double precision function bilinear_interp(x1, x2, y1, y2, x, y, f11, f21, f12, f22)
! Bilinear interpolation

    double precision, intent(in) :: x1, x2, y1, y2, x, y, f11, f21, f12, f22

    bilinear_interp = &
         1./((x2 - x1)*(y2 - y1)) * ( &
         f11*(x2 - x )*(y2 - y ) + & 
         f21*(x  - x1)*(y2 - y ) + & 
         f12*(x2 - x )*(y  - y1) + &
         f22*(x  - x1)*(y  - y1) )

    return
    end function bilinear_interp

!---------------------------------------------------------------------
    double precision function ellE_green(X,  DL) ! gives back the first kind elliptic integral,  from K. Lackner,  T. Lunt,  IPP - Garching 2022

    double precision,  intent(in) :: X,  DL

    ellE_green = (((0.01736506451D0 *X + 0.04757383546D0)*X + 0.06260601220D0)*X + 0.44325141463D0)*X + 1.0D0 - &
                 (((0.00526449639D0 *X + 0.04069697526D0)*X + 0.09200180037D0)*X + 0.24998368310D0)*X*DL

    return
    end function ellE_green

!---------------------------------------------------------------------
    double precision function ellK_green(X,  DL) ! gives back the second kind elliptic integral,  from K. Lackner,  T. Lunt,  IPP - Garching 2022

    double precision,  intent(in) :: X,  DL

    ellK_green = ((( 0.01451196212D0*X + 0.03742563713D0)*X + 0.03590092383D0)*X + 0.09666344259D0)*X + 1.38629436112D0 - &
                (((( 0.00441787012D0*X + 0.03328355346D0)*X + 0.06880248576D0)*X + 0.12498593597D0)*X + 0.5D0)*DL

    return 
    end function ellK_green

!----------------------------------------------------------------------------------- 
    double precision function green_function(r1, z1, r2, z2)

    double precision, intent(in) :: r1, z1, r2, z2
    double precision :: TT, K, ELCK, ELCE, acl, alg

    K = sqrt(4.*r1*r2/ ((r2 + r1)**2 + (z2 - z1)**2))

    TT = 1. - K**2

    acl = tt
    alg = dlog(acl)

    ELCK = ellK_green(acl, alg) !              S21BBF(0.D0, TT, 1.D0, IFAILK)
    ELCE = ellE_green(acl, alg) ! ELCK-K**2/3.D0*S21BCF(0.D0, TT, 1.D0, IFAILE)

    green_function = ( (1.D0 - K**2/2.)*ELCK - ELCE )*( SQRT(r1*r2)/K )

    return
    end function green_function

!-----------------------------------------------------------------------------------
    double precision function green_function_includingsamepoint(r1, z1, r2, z2, dl, R0)

    double precision, intent(in) :: r1, r2, z1, z2, dl, r0
    double precision :: TT

    if (abs(r1 - r2) < 1.e-6 .and. abs(z1 - z2) < 1.e-6) then
        tt = dl/(4*R0)
        green_function_includingsamepoint = -tt*(log(tt**2) - 4*log(2.) + 2.)*r0**2/dl
    else
        green_function_includingsamepoint = green_function(r1, z1, r2, z2)
    endif

    return
    end function green_function_includingsamepoint

!-----------------------------------------------------------------------------------
    double precision function green_function_identity(r1, z1, dr, dz, ntype)

    integer, intent(in) :: ntype
    double precision, intent(in) :: r1, z1, dr, dz
    integer :: i, j
    double precision :: ELCK, ELCE, greenf

    if (ntype == 2) then
        green_function_identity = r1*(log(8.*r1/(0.2236*(dr + dz))) - 2.0) ! From A. Kavin,  used in SPIDER (A. A. Ivanov and S. Yu. Medvedev),  self inductance of a rectangular coil in toroidal direction
    else if (ntype == 1) then
        elce = SQRT(dr**2 + dz**2)
        greenf = 0.
        do i=-1, 1
            do j=-1, i
                if (i == j) then
                    greenf = greenf + (r1 + dr*i/3.d0) * (log( 8.*(r1 + dr*i/3.d0)/(elce/3.d0) ) - 0.5D0)
                else
                    elck = green_function(r1 + dr*i/3., dz*i/3., r1 + dr*j/3., dz*j/3.)
                    greenf = greenf + 4.*elck
                endif
            enddo
        enddo
        green_function_identity = greenf/9.d0
    endif

    return
    end function green_function_identity

!-----------------------------------------------------------------------------------
    double precision function green_function_non_identity(r1, z1, r2, z2, dr, dz, dr2, dz2, ntype1, ntype2)

    integer, intent(in) :: ntype1, ntype2
    double precision, intent(in) :: r1, z1, r2, z2, dr, dz, dr2, dz2

    integer :: i, j
    double precision :: ELCK, greenf

    greenf = 0.

    if (ntype1 == 1 .and. ntype2 == 1) then
        do i=-1, 1
            do j=-1, 1
                elck = green_function(r1+dr*i/3., z1+dz*i/3., r2+dr2*j/3., z2+dz2*j/3.)
                greenf = greenf + elck
            enddo
        enddo
        green_function_non_identity = greenf/9.d0
    else if (ntype1 == 1 .and. ntype2 == 2) then
        do i=-1, 1
            elck = green_function(r1+dr*i/3., z1+dz*i/3., r2, z2)
            greenf = greenf + elck
        enddo
        green_function_non_identity = greenf/3.d0
    else if (ntype1 == 2 .and. ntype2 == 1) then
        do j=-1, 1
            elck = green_function(r1, z1, r2+dr2*j/3., z2+dz2*j/3.)
            greenf = greenf + elck
        enddo
        green_function_non_identity = greenf/3.d0
    else if (ntype1 == 2 .and. ntype2 == 2) then
        green_function_non_identity = green_function(r1, z1, r2, z2)
    endif

    return
    end function green_function_non_identity

!-----------------------------------------------------------------------------------
    function boundary(green_in) result(green_out)

! new bc is integral_over_boundary of -Green * dg/dn * dl
    use pi_vars, only: GPI
    use feqis_circuit, only: nr1, nz1, nr2, nz2, i_dim2

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

! calculates  integral_over_boundary of -Green * dg/dn * dl
    use feqis_circuit, only: nr1, nr2, nz1, i_dim2, green_bnd_f, r, dr, dz

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

    use feqis_circuit, only: nr, nz, psirz

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
        z1 = interp2d_psi(t1, t2, psirz(1:nr, 1:nz))
        z2 = interp2d_psi(t3, t4, psirz(1:nr, 1:nz))
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

!---------------------------------------------------------------------
    subroutine find_closest_xpoints(rx, zx, ierr, n_add)

    use feqis_circuit,  only: rbnd, zbnd, nteta, r, z, nr1, nz1, & 
        lim_maxR, lim_minR, lim_minZ, lim_maxZ, dr, dz

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

!--------------------------------------------------------------------
    function solve_tridiag_fbe(A, B, C, R, Ngrid) result(f_out)

! Provides solution of the system:
!
!   Aj fj - 1  + Bj fj  + Cj fj + 1 = Rj
!
!   where j = 1..Ngrid
!   bcbound = 1  -> given f_NA1
!   eximp = 2: implicit

    integer, intent(in) :: Ngrid
    double precision, intent(in) , dimension(Ngrid) :: A, B, C, R
    double precision, dimension(Ngrid) :: f_out

    integer :: j, k
    double precision, dimension(Ngrid) :: alpha, beta

    alpha(1) = C(1)/B(1)
    beta (1) = R(1)/B(1)
    do j=2, Ngrid-1
        alpha(j) = C(j)/(B(j) - A(j)*alpha(j-1))
    enddo
    do j=2, Ngrid
        beta(j) = (R(j) - A(j)*beta(j-1))/(B(j) - A(j)*alpha(j-1))
    enddo

! Note that boundary value is assumed to be on the main last grid point, so 1 - dx/2. This has to be
! corrected later on... CEfable
!     f(j-1) = (f_bound - beta(j-1))/alpha(j-1)
! This one should be appropriate with extrapolation... but now go back to real b.c.
!     f(j-1) = (2./3.*f_bound - beta(j-1))/(alpha(j-1)-1./3.)
    f_out(Ngrid) = beta(Ngrid)
    do k=1, Ngrid-1
        j = Ngrid - k
        f_out(j) = beta(j) - alpha(j)*f_out(j + 1)
    enddo

    return
    end function solve_tridiag_fbe

!--------------------------------------------------------------------
    double precision function fill_in_current(r0, nx, ppp_2d, ffp_2d, un)

    integer, intent(in) :: nx
    double precision, intent(in) :: r0, un
    double precision, intent(in), dimension(nx) ::  ppp_2d, ffp_2d

    integer :: k
    double precision :: zeta

    if (un >  1.) then
        fill_in_current = 0.0
    elseif (un <= 0.) then
        fill_in_current = ppp_2d(1 )*r0 + ffp_2d(1 )/r0
    elseif (un == 1.) then
        fill_in_current = ppp_2d(nx)*r0 + ffp_2d(nx)/r0
    else
        zeta = un*(nx - 1.) + 1.
        k = min(floor(zeta), nx - 1)
        k = max(k, 1)
        fill_in_current = ( &
             (ppp_2d(k + 1)*(zeta - k) + ppp_2d(k)*(k + 1. - zeta))*r0 + &
             (ffp_2d(k + 1)*(zeta - k) + ffp_2d(k)*(k + 1. - zeta))/r0   )
    endif

    return
    end function fill_in_current

!--------------------------------------------------------------------
    function solve_circuit_equations(nc, im, rm, I0, V, dpc, tau, invertcommand) result(cur_conduc)

    integer, intent(in) :: nc, invertcommand
    double precision, intent(in) :: tau
    double precision, intent(in), dimension(nc) :: I0, V, dpc
    double precision, intent(in), dimension(nc, nc) :: im, rm
    double precision, dimension(nc) :: cur_conduc

    integer :: i
    double precision, dimension(nc) :: B
    double precision, dimension(200, 200) :: invmatrix
    double precision, dimension(nc, nc) :: matrix

    save invmatrix

! Equation is im*(i1 - i0)/tau + rm*i1 = v - dpc

    do i=1, nc
        b(i) = v(i) - dpc(i) + sum(im(i, 1:nc)*i0(1:nc))/tau
    enddo

    if (invertcommand == 1) then
        matrix(1:nc, 1:nc) = im(1:nc, 1:nc)/tau + rm(1:nc, 1:nc)
        invmatrix(1:nc, 1:nc) = inv_matrix(matrix, nc)
    endif

    do i=1, nc
        cur_conduc(i) = sum(invmatrix(i, 1:nc)*b(1:nc))
    enddo

    return
    end function solve_circuit_equations


end module feqis_tools
