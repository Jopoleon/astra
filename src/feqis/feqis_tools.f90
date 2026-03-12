module feqis_tools

implicit none

double precision, dimension(:), allocatable :: costable
double precision, dimension(:, :), allocatable :: sintable

contains

!---------------------------------------------------------------------
    double precision function pol_angle(r0, z0, r, z)

    use pi_vars, only: GPI2

    double precision, intent(in) :: r0, z0, r, z

    pol_angle = ATAN2(z - z0, r - r0)
    if (pol_angle < 0) pol_angle = pol_angle + GPI2

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

    end function inv_matrix

!---------------------------------------------------------------------
    function discrete_sine_transform(ndim, f_in) result(f_out)

    integer, parameter :: dp=selected_real_kind(15, 300)

    integer, intent(in) :: ndim
    double precision,  intent(in), dimension(ndim) :: f_in
    double precision, dimension(ndim) :: f_out

    integer :: i, j, imethod1, icall, k
    double precision :: z(ndim)
    complex(kind=dp) :: d(2*(ndim+1))

    imethod1 = 1
    if (imethod1 == 1) then
        do i=1, ndim
            f_out(i) = SUM(f_in(1:ndim)*sintable(i, 1:ndim))
        enddo
    else if (imethod1 == 2) then
! fast sine transform
! this problem is equivalent to DST-I with N = n+1
        k = 2*(ndim+1)
        z = f_in
        d(1) = cmplx(0., 0.)
        do i=1, ndim
            d(i+1)   = cmplx( z(i), 0)
            d(k-i+1) = cmplx(-z(i), 0)
        enddo
        d(k-ndim) = cmplx(0., 0.)
        d(1: k-1) = d(2: k)
        do i=1, ndim
            f_out(i) = 0.25*aimag(d(i) - d(k-i))
        enddo
        f_out = f_out/2.
    endif

    end function discrete_sine_transform

!---------------------------------------------------------------------
    function least_square_biquad(r, z, u, n) result(derivs)

    integer, intent(in) :: n
    double precision, intent(in), dimension(n) :: r, z, u
    double precision, dimension(5) :: derivs

    integer :: k
    double precision :: det, det_r, det_z, rax, zax, uax
    double precision, dimension(6) :: B, coeff
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
    Ainv = inv_matrix(A, 6)

    do k=1, 6
        coeff(k) = -2*sum(Ainv(k, 1: 6)*B(1: 6))
    enddo

    derivs(1) = coeff(4)
    derivs(2) = coeff(5)
    derivs(3) = 2.*coeff(1)
    derivs(4) = 2.*coeff(2)
    derivs(5) = coeff(3)

    end function least_square_biquad

!---------------------------------------------------------------------
    subroutine expandCoeffs(coeff, r_in, z_in, u_out, deriv_out)

    double precision, intent(in) :: r_in, z_in
    double precision, intent(in), dimension(9) :: coeff
    double precision, intent(out) :: u_out
    double precision, intent(out), dimension(5) :: deriv_out

    u_out = coeff(1)*r_in**2 * z_in**2 + &
            coeff(2)*r_in**2 * z_in    + &
            coeff(3)*r_in    * z_in**2 + &
            coeff(4)*r_in    * z_in    + &
            coeff(5)*r_in**2 + &
            coeff(6)*z_in**2 + &
            coeff(7)*r_in    + &
            coeff(8)*z_in    + &
            coeff(9)

    deriv_out(1) = 2*coeff(1)*r_in*z_in**2 + 2*coeff(2)*r_in*z_in + coeff(3)*z_in**2   + coeff(4)*z_in + 2*coeff(5)*r_in + coeff(7)
    deriv_out(2) = 2*coeff(1)*r_in**2*z_in +   coeff(2)*r_in**2 + 2*coeff(3)*z_in*r_in + coeff(4)*r_in + 2*coeff(6)*z_in + coeff(8)
    deriv_out(3) = 2.*coeff(1)*z_in**2 + 2*coeff(2)*z_in + 2*coeff(5)
    deriv_out(4) = 2.*coeff(1)*r_in**2 + 2*coeff(3)*r_in + 2*coeff(6)
    deriv_out(5) = 4*coeff(1)*r_in*z_in + 2*coeff(2)*r_in + 2*coeff(3)*z_in + coeff(4)

    end subroutine expandCoeffs

!---------------------------------------------------------------------
    double precision function getHessian(derivs)

    double precision, intent(in), dimension(5) :: derivs
    getHessian = derivs(3)*derivs(4) - derivs(5)**2

    end function getHessian

!---------------------------------------------------------------------
    function getDerivs(derivs_in, dr, dz) result(derivs_out)

    double precision, intent(in) :: dr, dz
    double precision, dimension(5), intent(in) :: derivs_in

    double precision, dimension(5) :: derivs_out

    derivs_out(1) = derivs_in(1)/dr
    derivs_out(2) = derivs_in(2)/dz
    derivs_out(3) = derivs_in(3)/dr**2
    derivs_out(4) = derivs_in(4)/dz**2
    derivs_out(5) = derivs_in(5)/(dr*dz)
    
    end function getDerivs

!---------------------------------------------------------------------
    subroutine transformRZ(r_in, z_in, deriv_in, r_out, z_out)

    double precision, intent(in) :: r_in, z_in
    double precision, intent(in), dimension(5) :: deriv_in
    double precision, intent(out) :: r_out, z_out

    double precision :: hessian

    hessian = getHessian(deriv_in)
    r_out = r_in - (deriv_in(4)*deriv_in(1) - deriv_in(5)*deriv_in(2))/hessian
    z_out = z_in - (deriv_in(3)*deriv_in(2) - deriv_in(5)*deriv_in(1))/hessian

    end subroutine transformRZ
    
!---------------------------------------------------------------------
    subroutine exact_biquad(coeff, dr_out, dz_out, u_out, hessian, dr, dz)

    use errors_params, only: err_find_biquad

    integer, parameter :: n_iter=100000

    double precision, intent(in) :: dr, dz
    double precision, intent(in), dimension(9) :: coeff
    double precision, intent(out) :: dr_out, dz_out, u_out
    double precision, intent(out) :: hessian

    integer :: jiter
    double precision :: r_loc, z_loc
    double precision, dimension(5) :: d_dpsi, derivs

    r_loc = 0.
    z_loc = 0.
    do jiter=1, n_iter
        call expandCoeffs(coeff, r_loc, z_loc, u_out, d_dpsi)
        if (abs(d_dpsi(1)) < err_find_biquad .and. abs(d_dpsi(2)) < err_find_biquad) then ! Convergence
            dr_out = r_loc*dr
            dz_out = z_loc*dz
            derivs = getDerivs(d_dpsi, dr, dz)
            hessian = getHessian(derivs)
            return
        endif
        call transformRZ(r_loc, z_loc, d_dpsi, r_loc, z_loc)
        if (abs(r_loc) > 1. .or. abs(z_loc) > 1) then
            EXIT
        endif
    enddo
    dr_out  =  1.e6
    dz_out  =  1.e6
    u_out   = -1.e6
    hessian =  1.e6

    end subroutine exact_biquad

!---------------------------------------------------------------------
    subroutine exact_biquad_regress(coeff, dr_out, dz_out, u_out, derivs, dr, dz)

    double precision, intent(in) :: dr, dz
    double precision, intent(in), dimension(9) :: coeff
    double precision, intent(out) :: dr_out, dz_out, u_out
    double precision, intent(out), dimension(5) :: derivs

    double precision :: hessian, r_loc, z_loc
    double precision :: d_dpsi(5)

    r_loc = 0.
    z_loc = 0.
    call expandCoeffs(coeff, r_loc, z_loc, u_out, d_dpsi)
    call transformRZ(r_loc, z_loc, d_dpsi, r_loc, z_loc)
    call expandCoeffs(coeff, r_loc, z_loc, u_out, d_dpsi)
    derivs = getDerivs(d_dpsi, dr, dz)
    dr_out = r_loc*dr
    dz_out = z_loc*dz

    end subroutine exact_biquad_regress

!---------------------------------------------------------------------
    integer function closest_index(x_in, xmin, dx)

    double precision, intent(in) :: x_in, xmin, dx

    closest_index = nint((x_in - xmin)/dx + 1.) ! nint(1.8) = 2

    end function closest_index

!---------------------------------------------------------------------
    integer function floor_index(x_in, xmin, dx)

    double precision, intent(in) :: x_in, xmin, dx

    floor_index = floor((x_in - xmin)/dx + 1.) ! floor(1.8) = 1

    end function floor_index

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

    end function bilinear_interp

!---------------------------------------------------------------------
    double precision function ellE_green(X,  DL) ! gives back the first kind elliptic integral,  from K. Lackner,  T. Lunt,  IPP - Garching 2022

    double precision,  intent(in) :: X,  DL

    ellE_green = (((0.01736506451D0 *X + 0.04757383546D0)*X + 0.06260601220D0)*X + 0.44325141463D0)*X + 1.0D0 - &
                 (((0.00526449639D0 *X + 0.04069697526D0)*X + 0.09200180037D0)*X + 0.24998368310D0)*X*DL

    end function ellE_green

!---------------------------------------------------------------------
    double precision function ellK_green(X,  DL) ! gives back the second kind elliptic integral,  from K. Lackner,  T. Lunt,  IPP - Garching 2022

    double precision,  intent(in) :: X,  DL

    ellK_green = ((( 0.01451196212D0*X + 0.03742563713D0)*X + 0.03590092383D0)*X + 0.09666344259D0)*X + 1.38629436112D0 - &
                (((( 0.00441787012D0*X + 0.03328355346D0)*X + 0.06880248576D0)*X + 0.12498593597D0)*X + 0.5D0)*DL

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
    double precision :: den
    double precision, dimension(Ngrid) :: alpha, beta

    alpha(1) = C(1)/B(1)
    beta (1) = R(1)/B(1)
    do j=2, Ngrid
        den = B(j) - A(j)*alpha(j-1)
        alpha(j) = C(j)/den
        beta(j) = (R(j) - A(j)*beta(j-1))/den
    enddo

! Note that boundary value is assumed to be on the main last grid point, so 1 - dx/2. This has to be
! corrected later on... CEfable
!     f(j-1) = (f_bound - beta(j-1))/alpha(j-1)
! This one should be appropriate with extrapolation... but now go back to real b.c.
!     f(j-1) = (2./3.*f_bound - beta(j-1))/(alpha(j-1)-1./3.)
    f_out(Ngrid) = beta(Ngrid)
    do k=1, Ngrid-1
        j = Ngrid - k
        f_out(j) = beta(j) - alpha(j)*f_out(j+1)
    enddo

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

    end function solve_circuit_equations

end module feqis_tools
