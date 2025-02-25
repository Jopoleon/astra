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

!        call fft_eff(d)

        d(1: k-1) = d(2: k)
        do i=1, ndim
            z(i) = 0.5*aimag(d(i) - d(k-i))
        enddo
        f_out = z/2.
    endif

    return
    end function discrete_sine_transform

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
    integer function ceil_index(x_in, xmin, dx)

    double precision, intent(in) :: x_in, xmin, dx

    ceil_index = ceiling((x_in - xmin)/dx + 1.) ! ceil(1.8) = 2

    return
    end function ceil_index

!---------------------------------------------------------------------
    double precision function interp2d_psi(r_in, z_in, Rgrid, Zgrid, psi_in)
! ->psi at r_in, z_in

    double precision, intent(in) :: r_in, z_in
    double precision, intent(in), dimension(:) :: Rgrid, Zgrid
    double precision, intent(in), dimension(:, :) :: psi_in

    integer :: i, j
    integer, dimension(2) :: psi_shape
    double precision :: dr, dz, r1, r2, z1, z2, psi1, psi2, psi3, psi4

    psi_shape = SHAPE(psi_in)
    dr = Rgrid(2) - Rgrid(1)
    dz = Zgrid(2) - Zgrid(1)

    i = floor_index(r_in, Rgrid(1), dr)
    j = floor_index(z_in, Zgrid(1), dz)
    i = min(psi_shape(1) - 1, max(1, i))
    j = min(psi_shape(2) - 1, max(1, j))

    r1 = Rgrid(i)
    r2 = Rgrid(i+1)
    z1 = Zgrid(j)
    z2 = Zgrid(j+1)
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
    double precision, dimension(500, 500) :: invmatrix
    double precision, dimension(nc, nc) :: matrix

    save invmatrix

! Equation is im*(i1 - i0)/tau + rm*i1 = v - dpc

    do i=1, nc
        b(i) = tau*(v(i) - dpc(i)) + sum(im(i, 1:nc)*i0(1:nc))
    enddo

    if (invertcommand == 1) then
        matrix(1:nc, 1:nc) = im(1:nc, 1:nc) + tau*rm(1:nc, 1:nc)
        invmatrix(1:nc, 1:nc) = inv_matrix(matrix, nc)
    endif

    do i=1, nc
        cur_conduc(i) = sum(invmatrix(i, 1:nc)*b(1:nc))
    enddo

    return
    end function solve_circuit_equations


end module feqis_tools
