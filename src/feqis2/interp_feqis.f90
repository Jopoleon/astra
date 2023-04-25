subroutine qinterp_feqis(x1, y1, Nx1, x2, y2, Nx2)

implicit none

integer, intent(in) :: Nx1, Nx2
double precision, intent(in) , dimension(Nx1) :: x1, y1
double precision, intent(in) , dimension(Nx2) :: x2
double precision, intent(out), dimension(Nx2) :: y2

integer :: i, j
double precision :: A, B, C, z1, z2, z3, t1, t2, t3, t4

do i=1, Nx2
    t4 = x2(i)

    do j=2, Nx1-1
        z1 = x1(j-1)
        z2 = x1(j)
        z3 = x1(j+1)

        if (t4 == z1) then
            y2(i) = y1(j-1)
            EXIT
        elseif (t4 == z2) then
            y2(i) = y1(j)
            EXIT
        elseif (t4 == z3) then
            y2(i) = y1(j+1)
            EXIT
        elseif (t4 > z1 .and. t4 < z3) then
            t1 = y1(j-1)
            t2 = y1(j)
            t3 = y1(j+1)
            A = (t3 - t2 - (z3 - z2)*(t1 - t2)/(z1 - z2)) / ((z3 - z2)*(z3 - z1))
            B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
            C = t2 - A*z2**2 - B*z2
            y2(i) = A*t4**2 + B*t4 + C
            EXIT
        elseif (t4 < z1 .and. j == 2) then
            z1 = x1(j-1)**2
            z2 = x1(j)**2
            t1 = y1(j-1)
            t2 = y1(j)
            A = 0.0
            B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
            C = t2 - A * z2**2 - B*z2
            y2(i) = A * t4**4 + B * t4**2 + C
            EXIT
        elseif (t4 > z3 .and. j == Nx1-1) then
            t1 = y1(j-1)
            t2 = y1(j)
            t3 = y1(j+1)
            A = (t3 - t2 - (z3 - z2)*(t1 - t2)/(z1 - z2)) / ((z3 - z2)*(z3 - z1))
            B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
            C = t2 - A * z2**2 - B*z2
            y2(i) = A * t4**2 + B*t4 + C
            EXIT
        endif
    enddo
enddo

return
end subroutine qinterp_feqis

!---------------------------------------------------------------------
subroutine linterp_feqis(x1, y1, Nx1, x2, y2, Nx2)

implicit none

integer, intent(in) :: Nx1, Nx2
double precision, intent(in) , dimension(Nx1) :: x1, y1
double precision, intent(in) , dimension(Nx2) :: x2
double precision, intent(out), dimension(Nx2) :: y2

integer :: i, j
double precision :: A, C, z1, z2, t1, t2, t4

do i=1, Nx2
    t4 = x2(i)
    do j=2, Nx1
        z1 = x1(j-1)
        z2 = x1(j)
        if (t4 == z1) then
            y2(i)=y1(j-1)
            EXIT
        elseif (t4 == z2) then
            y2(i)=y1(j)
            EXIT
        elseif (t4 > z1 .and. t4 < z2) then
            t1 = y1(j-1)
            t2 = y1(j)
            A = (t2 - t1)/(z2 - z1)
            C = t2 - A*z2
            y2(i) = A*t4 + C
            EXIT
        elseif (t4 < z1 .and.j == 2) then
            y2(i)=y1(j-1)
            EXIT
        elseif (t4 > z2 .and. j == nx1) then
            y2(i) = y1(j)
            EXIT
        endif
    enddo
enddo

return
end subroutine linterp_feqis

!---------------------------------------------------------------------
subroutine find_angle_ef(rt, zt, r, z, anglr)

use pi_grec_vars, only: GPI2

implicit none
double precision, intent(in) :: rt, zt, r, z
double precision, intent(out) :: anglr

anglr = ATAN2(z - zt, r - rt)
if (anglr < 0) anglr = anglr + GPI2

return
end subroutine find_angle_ef

!---------------------------------------------------------------------
subroutine polyfitcc_feqis(x, y, P)

implicit none

double precision, intent(in) , dimension(3) :: x, y
double precision, intent(out), dimension(3) :: P

double precision :: y21, y32, x21, x32, h21, h32

y32 = y(3) - y(2)
y21 = y(2) - y(1)
x32 = x(3) - x(2)
x21 = x(2) - x(1)
h32 = x(3) + x(2)
h21 = x(2) + x(1)

P(1) = (x21*y32 - x32*y21)/(x21*x32*(h32 - h21))
P(2) = y21/x21 - P(1)*h21
P(3) = y(3) - P(1)*x(3)**2 - P(2)*x(3)

return
end subroutine polyfitcc_feqis

!---------------------------------------------------------------------
subroutine inverse_matrix_equilef(mat, mat_inv, n)

! mat(n,n)     - array of coefficients for matrix A
! n            - dimension
! mat_inv(n,n) - inverse matrix of A
! Based on Doolittle LU factorization for Ax = B
! Alex G. December 2009. www2.odu.eud/~agodunov/computing/programs/book2/Ch06/Inverse.f90

implicit none 

integer, intent(in) :: n
double precision, intent(in) , dimension(n, n) :: mat
double precision, intent(out), dimension(n, n) :: mat_inv

integer :: i, j, k
double precision, dimension(n) :: b, d, x
double precision, dimension(n, n) :: a, L, U

L = 0.0
U = 0.0
b = 0.0
a = mat

! step 1: forward elimination
do k=1, n-1
    do i=k+1, n
        L(i, k) = a(i, k)/a(k, k)
        do j=k+1, n
            a(i, j) = a(i, j) - L(i, k)*a(k, j)
        enddo
    enddo
enddo

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1, n
    L(i, i) = 1.0
enddo
! U matrix is the upper triangular part of A
do j=1, n
    do i=1, j
        U(i, j) = a(i, j)
    enddo
enddo

! Step 3: compute columns of the inverse matrix C
do k=1, n
    b(k) = 1.0
    d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
    do i=2, n
        d(i) = b(i)
        do j=1, i-1
            d(i) = d(i) - L(i, j)*d(j)
        enddo
    enddo
! Step 3b: Solve Ux=d using the back substitution
    x(n) = d(n)/U(n, n)
    do i = n-1, 1, -1
        x(i) = d(i)
        do j=n, i+1, -1
            x(i) = x(i) - U(i, j)*x(j)
        enddo
        x(i) = x(i)/u(i, i)
    enddo
! Step 3c: fill the solutions x(n) into column k of C
    do i=1, n
        mat_inv(i, k) = x(i)
    enddo
    b(k) = 0.0
enddo

return
end subroutine inverse_matrix_equilef

!---------------------------------------------------------------------
subroutine interp_j_fromrhotorz

use ef_circuit, only: nrho, nteta, nr2, nz2, r, z, jrz, jrhoteta, rho, teta, raxp, zaxp

implicit none

integer i, j

! go from jrhoteta to jrz
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
subroutine curinterp_ef(r, z, jrho, rho, teta, rax, zax, nrho, nteta, j)

use pi_grec_vars, only: GPI2

implicit none

integer, intent(in) :: nrho, nteta
double precision, intent(in) :: r, z, rax, zax
double precision, intent(in), dimension(nteta) :: teta
double precision, intent(in), dimension(nrho, nteta) :: jrho, rho
double precision, intent(out) :: j

integer :: i, j1, j2, k1, k2, k3, k4
double precision :: anglr, rho0
double precision, dimension(1) :: jc
double precision, dimension(2) :: x1, y1, jt

call find_angle_ef(rax, zax, r, z, anglr)
rho0 = sqrt((r - rax)**2 + (z - zax)**2)
  
if (anglr < teta(1)) anglr = anglr + GPI2
do i=1, nteta-1
    if (anglr >= teta(i)) j1 = i
enddo
j2 = j1 + 1

if (rho0 > rho(nrho, j1) .or. rho0 > rho(nrho, j2)) then
    j = 0.
    return
endif

do i=1, nrho-1
    if (rho0 > rho(i, j1)) k1 = i
enddo
k2 = k1 + 1
do i=1, nrho-1
    if (rho0 > rho(i, j2)) k3 = i
enddo
k4 = k3 + 1

x1(1) = rho(k1, j1)
x1(2) = rho(k2, j1)
y1(1) = jrho(k1, j1)
y1(2) = jrho(k2, j1)

call linterp_feqis(x1, y1, 2, (/rho0/), jt(1:1), 1)
x1(1) = rho(k3, j2)
x1(2) = rho(k4, j2)
y1(1) = jrho(k3, j2)
y1(2) = jrho(k4, j2)

call linterp_feqis(x1, y1, 2, (/rho0/), jt(2:2), 1)
x1(1) = teta(j1)
x1(2) = teta(j2)
y1(1) = jt(1)
y1(2) = jt(2)
call linterp_feqis(x1, y1, 2, (/anglr/), jc, 1)

j = jc(1)

return
end subroutine curinterp_ef

!---------------------------------------------------------------------
subroutine discrete_sine_transform_ef(n, y)

use fft_mod_eff, only: sintable

implicit none

integer, parameter :: dp = kind(1.d0)

integer, intent(in) :: n
double precision, intent(inout), dimension(n) :: y
integer :: i, j, imethod1, k
double precision, dimension(n) :: z
complex(kind=dp), dimension(2*(n+1)) :: d

imethod1 = 1

if (imethod1 == 1) then
    z = 0.
    do i=1, n
        do j=1, n
            z(i) = z(i) + y(j)*sintable(i, j)
        enddo
    enddo
    y = z
else if (imethod1 == 2) then
! fast sine transform 
! this problem is equivalent to DST-I with N = n+1 
    k = 2*(n + 1)
    z = y
    d(1) = cmplx(0., 0.)
    do i=1, n
        d(  i+1) = cmplx( z(i), 0)
        d(k-i+1) = cmplx(-z(i), 0)
    enddo
    d(k-n) = cmplx(0., 0.)
    call fft_eff(d)
    d(1:k-1) = d(2:k)
    do i=1, n
        z(i) = 0.5*aimag(d(i) - d(k-i))
    enddo
    y = 0.5*z
endif

return
end subroutine discrete_sine_transform_ef

!---------------------------------------------------------------------
subroutine plasma_psi_to_coils_ef

use ef_circuit, only: nconduc, nr2, nz2, psiplasmatoconduc, jrz, area_eff
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
subroutine psiplex_calc_ef(dumz)

use pi_grec_vars, only: GPI2
use exchange_with_astra, only: psplex_from_fbe
use ef_circuit, only: nbnd, dr, rbnd, dz, zbnd

implicit none

double precision, intent(inout) :: dumz
integer :: i, j
double precision :: t1, t2, t3, t4, z1, z2, z3, z4, &
    x1, x2, x3, y1, y2, y3, arc1, arc2

write(*, *) 'spid par', psplex_from_fbe

if (psplex_from_fbe == 1) then
!psplexs = psplexs/(btor*rhos(iplas)**2)*qedge
    z4 = 0.
    z2 = 0.
    do j=1, nbnd-1
        x1 = rbnd(j)
        y1 = zbnd(j)
        x3 = rbnd(j+1)
        y3 = zbnd(j+1)
        x2 = 0.5*(rbnd(j+1) + rbnd(j))
        y2 = 0.5*(zbnd(j+1) + zbnd(j))
        arc2 = sqrt((x3 - x1)**2 + (y3 - y1)**2)
        do i=1, nbnd-1
            x1 = rbnd(i)
            y1 = zbnd(i)
            x3 = rbnd(i+1)
            y3 = zbnd(i+1)
            arc1 = sqrt((x3-x1)**2 + (y3-y1)**2)
            call green_function(x1, y1, x2, y2, z3)
            call find_fields_interp_ef_psionly(rbnd(i) + dr/2, zbnd(i), t1)
            call find_fields_interp_ef_psionly(rbnd(i), zbnd(i) + dz/2, t2)
            call find_fields_interp_ef_psionly(rbnd(i) - dr/2, zbnd(i), t3)
            call find_fields_interp_ef_psionly(rbnd(i), zbnd(i) - dz/2, t4)
            z1 = sqrt(((t3 - t1)/dr)**2 + ((t4 - t2)/dz)**2)
            z2 = z2 + z3/x1*z1*arc1*arc2
        enddo
        z4 = z4 + arc2
    enddo
    dumz = z2/z4*GPI2
else
    dumz = dumz
endif

return
end subroutine psiplex_calc_ef

!---------------------------------------------------------------------
subroutine psi_external_calc_ef

use ef_circuit, only: nr2, nz2, nconduc, psiextrz, curconduc
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
subroutine least_square_biquad_ef(r, z, u, n, c, rax, zax, uax, derivs)

implicit none

integer, intent(in) :: n
double precision, intent(in), dimension(n) :: r, z, u
double precision, intent(out) :: rax, zax, uax
double precision, intent(out), dimension(6) :: c
double precision, intent(out), dimension(8) :: derivs

integer :: k
double precision :: det, det_r, det_z
double precision, dimension(6) :: B, cc
double precision, dimension(21) :: sums
double precision, dimension(6, 6) :: A, Ainv

sums = 0.
sums( 1) = sum(r**4)
sums( 2) = sum(z**4)
sums( 3) = sum(r**2 * z**2)
sums( 4) = sum(r**2)
sums( 5) = sum(z**2)
sums( 6) = n + 0.
sums( 7) = sum(r**3 * z)
sums( 8) = sum(r**3)
sums( 9) = sum(r**2 * z)
sums(10) = sum(r * z**3)
sums(11) = sum(r * z**2)
sums(12) = sum(z**3)
sums(13) = sum(r*z)
sums(14) = sum(r)
sums(15) = sum(z)
sums(16) = sum(u * r**2)
sums(17) = sum(u * z**2)
sums(18) = sum(u * r*z)
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
call inverse_matrix_equilef(A, Ainv, 6)

do k=1, 6
    cc(k) = -2*sum(Ainv(k, 1:6)*B(1:6))
enddo
c(4) = cc(1)
c(5) = cc(2)
c(6) = cc(3)
c(2) = cc(4)
c(3) = cc(5)
c(1) = cc(6)

! magnetic axis

det = 4.d0*cc(1)*cc(2)-cc(3)**2
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
end subroutine least_square_biquad_ef

!---------------------------------------------------------------------
subroutine exact_biquad_ef(r, z, u, rax, zax, uax, derivs, dr, dz)

use errors_params, only: err_find_biquad

implicit none

integer, parameter :: ndim=9

double precision, intent(in) :: dr, dz
double precision, intent(in), dimension(ndim) :: r, z, u
double precision, intent(out) :: rax, zax, uax
double precision, intent(out), dimension(ndim-1) :: derivs

integer :: k, niter, j_success
double precision :: s_r, s_z, s_r2, s_z2, det, tolez
double precision, dimension(ndim) :: x, y, B, C
double precision, dimension(ndim, ndim) :: A, Ainv

tolez = err_find_biquad

!transformation
x = (r - r(5))/dr
y = (z - z(5))/dz

!find coefficients
call ainv_matrix_def(Ainv)

do k=1, ndim
    c(k) = sum(Ainv(k, 1:ndim)*u(1:ndim))
enddo

rax = x(5)
zax = y(5)

!now find axis
s_r = 100000.
s_z = 100000.
j_success = 0
do niter=1, 100000
    s_r2 = 2*C(1)*rax*zax**2 + 2*C(2)*rax*zax + C(3)*zax**2 + C(4)*zax + 2*C(5)*rax + C(7)
    s_z2 = 2*C(1)*rax**2*zax + C(2)*rax**2 + 2*C(3)*zax*rax + C(4)*rax + 2*C(6)*zax + C(8)
    A(1, 1) = 2*C(1)*zax**2 + 2*C(2)*zax + 2*C(5)
    A(1, 2) = 4*C(1)*rax*zax + 2*C(2)*rax + 2*C(3)*zax + C(4)
    A(2, 2) = 2*C(1)*rax**2 + 2*C(3)*rax + 2*c(6)
    A(2, 1) = 4*C(1)*rax*zax + 2*C(2)*rax + 2*C(3)*zax + C(4)
    B(1) = s_r2
    B(2) = s_z2
    det = (A(1, 1)*A(2, 2))-(A(1, 2)*A(2, 1))
    s_r2 = 1/det*(A(2, 2)*B(1)-A(1, 2)*B(2))
    s_z2 = 1/det*(A(1, 1)*B(2)-A(2, 1)*B(1))
    rax = rax-s_r2
    zax = zax-s_z2

    s_r2 = s_r
    s_z2 = s_z
    s_r = 2*C(1)*rax*zax**2 + 2*C(2)*rax*zax + C(3)*zax**2 + C(4)*zax + 2*C(5)*rax + C(7)
    s_z = 2*C(1)*rax**2*zax + C(2)*rax**2 + 2*C(3)*zax*rax + C(4)*rax + 2*C(6)*zax + C(8)

    if (abs(s_r) < tolez .and. abs(s_z) < tolez) then
        j_success = 1
        EXIT
    endif
    if (abs(rax) > 1. .or. abs(zax) > 1) then
        EXIT
    endif
enddo

if (j_success == 0) then
    rax    =  1.e6
    zax    =  1.e6
    uax    = -1.e6
    derivs =  1.e6
else
! magnetic axis
    uax = c(1) * rax**2 * zax**2 + c(2) * rax**2 * zax + c(3)*rax * zax**2 + c(4)*rax*zax + &
          c(5) * rax**2 + c(6) * zax**2 + c(7)*rax + c(8)*zax + c(9)

    derivs(1) = 1/dr*(2*C(1)*rax * zax**2 + 2*C(2)*rax*zax  + C(3) * zax**2 + C(4)*zax + 2*C(5)*rax + C(7))
    derivs(2) = 1/dz*(2*C(1) * rax**2 * zax + C(2) * rax**2 + 2*C(3)*zax*rax + C(4)*rax + 2*C(6)*zax + C(8))
    derivs(3) = 1/dr**2 * (2.*c(1) * zax**2 + 2*c(2)*zax + 2*c(5))
    derivs(4) = 1/dz**2 * (2.*c(1) * rax**2 + 2*c(3)*rax + 2*c(6))
    derivs(5) = 1/dr/dz * (4*c(1)*rax*zax + 2*c(2)*rax + 2*c(3)*zax + c(4)) 

    rax = rax*dr + r(5)
    zax = zax*dz + z(5)
endif

return
end subroutine exact_biquad_ef
         
!---------------------------------------------------------------------
subroutine exact_biquad_regress_ef(r, z, u, rax, zax, uax, derivs, dr, dz, rx, zx)

implicit none

integer, parameter :: ndim=9

double precision, intent(in) :: dr, dz, rx, zx
double precision, intent(in), dimension(ndim) :: r, z, u
double precision, intent(out) :: rax, zax, uax
double precision, intent(out), dimension(ndim-1) :: derivs

integer :: k
double precision :: s_r, s_z, s_r2, s_z2, xx, yy, det
double precision, dimension(ndim) :: x, y, B, C
double precision, dimension(ndim, ndim) :: A, Ainv

!transformation
x  = (r  - r(5))/dr
y  = (z  - z(5))/dz
xx = (rx - r(5))/dr
yy = (zx - z(5))/dz


!find coefficients
call ainv_matrix_def(Ainv)

do k=1, ndim
    c(k) = sum(Ainv(k, 1:ndim)*u(1:ndim))
enddo

rax = xx
zax = yy

!now find axis
s_r2 = 2*C(1)*rax * zax**2 + 2*C(2)*rax*zax + C(3) * zax**2 + C(4)*zax + 2*C(5)*rax + C(7)
s_z2 = 2*C(1) * rax**2 * zax + C(2) * rax**2 + 2*C(3)*zax*rax + C(4)*rax + 2*C(6)*zax + C(8)
A(1, 1) = 2*C(1) * zax**2 + 2*C(2)*zax + 2*C(5)
A(1, 2) = 4*C(1)*rax*zax  + 2*C(2)*rax + 2*C(3)*zax + C(4)
A(2, 2) = 2*C(1) * rax**2 + 2*C(3)*rax + 2*c(6)
A(2, 1) = 4*C(1)*rax*zax  + 2*C(2)*rax + 2*C(3)*zax + C(4)
B(1) = s_r2
B(2) = s_z2
det = (A(1, 1)*A(2, 2)) - (A(1, 2)*A(2, 1))
s_r2 = 1/det*(A(2, 2)*B(1) - A(1, 2)*B(2))
s_z2 = 1/det*(A(1, 1)*B(2) - A(2, 1)*B(1))
rax = rax - s_r2
zax = zax - s_z2

s_r = 2*C(1)*rax * zax**2 + 2*C(2)*rax*zax + C(3)*zax**2 + C(4)*zax + 2*C(5)*rax + C(7)
s_z = 2*C(1)*rax**2 * zax + C(2)*rax**2 + 2*C(3)*zax*rax + C(4)*rax + 2*C(6)*zax + C(8)

uax = c(1) * rax**2 * zax**2 + c(2) * rax**2 * zax + c(3)*rax * zax**2 + c(4)*rax*zax + &
      c(5) * rax**2 + c(6) * zax**2 + c(7)*rax +c(8)*zax + c(9)

derivs(1) = 1/dr*(2*C(1)*rax * zax**2 + 2*C(2)*rax*zax + C(3) * zax**2 + C(4)*zax + 2*C(5)*rax + C(7))
derivs(2) = 1/dz*(2*C(1) * rax**2 * zax + C(2) * rax**2 + 2*C(3)*zax*rax + C(4)*rax + 2*C(6)*zax + C(8))
derivs(3) = 1/dr**2*(2.*c(1) * zax**2 + 2*c(2)*zax + 2*c(5))
derivs(4) = 1/dz**2*(2.*c(1) * rax**2 + 2*c(3)*rax + 2*c(6))
derivs(5) = 1/dr/dz*(4*c(1)*rax*zax + 2*c(2)*rax + 2*c(3)*zax + c(4)) 
             
rax = rax*dr + r(5)
zax = zax*dz + z(5)

return
end subroutine exact_biquad_regress_ef

!---------------------------------------------------------------------             
real*8 function frlim_ef(dp, ylim, rx, zx, rm, zm)

implicit none

double precision, intent(in) :: ylim, rx, zx, rm, zm
double precision, intent(in), dimension(5) :: dp
double precision :: dxx, dxy, dyy, disc, cc, cdpls, cdmns, ang1, ang2, &
    c1, dl2x, c2, dl2y

Dxx = dp(3)
Dxy = dp(4)
Dyy = dp(5)
disc = (Dxy/Dyy)**2 - Dxx/Dyy

if (Disc < 0.) then
    cc = (zm-zx)/(rm-rx)
    cc = -1.d0/cc
else
    cdpls = -Dxy/Dyy + dsqrt(disc)
    cdmns = -Dxy/Dyy - dsqrt(disc)
    ang1  =  0.5d0*(datan(cdpls) + datan(cdmns))
    ang2  = -0.5d0*(datan(1.d0/cdpls) + datan(1.d0/cdmns))

! Calculate D2u/Dl2(direction ang1)
    c1 = dtan(ang1)
    Dl2x = Dyy*c1*c1 + 2.d0*Dxy*c1 + Dxx

! Calculate D2u/Dl2(direction ang2)
    c2 = dtan(ang2)
    Dl2y = Dyy*c2*c2 + 2.*Dxy*c2 + Dxx

    if (Dl2x < 0.) then
        cc = c1
    else if (Dl2y < 0.) then
        cc = c2
    else
        cc = (zm-zx)/(rm-rx)
        cc = -1.d0/cc
    endif
endif

frlim_ef = rx + (ylim-zx)/cc

return
end function frlim_ef

!---------------------------------------------------------------------
subroutine nine_point_regression(r0, z0, pos_xpoint, ddpsi, f00)

use ef_circuit, only: nr1, nz1, r, dr, z, dz, psirz

implicit none

integer, parameter :: ndim=9

double precision, intent(in) :: r0, z0
double precision, intent(out) :: f00
double precision, intent(out), dimension(2) :: pos_xpoint
double precision, intent(out), dimension(ndim-1) :: ddpsi

integer :: iax, jax, i, j, k, d, i1, i2, i3, i4
double precision, dimension(ndim) :: xub, bub, yub

call find_actual_index_ef(r0, z0, iax, jax)
!find true axis

i3 = -1
i1 = -1
i4 =  1
i2 =  1
if (iax ==   2) i3 = -1
if (jax ==   2) i1 = -1
if (iax == nr1) i4 =  1
if (jax == nz1) i2 =  1

d = (i4 - i3 + 1)*(i2 - i1 + 1)

k = 0
do j=i3, i4
    do i=i1, i2
        k = k + 1
        xub(k) = r(iax+i)
        yub(k) = z(jax+j)
        bub(k) = psirz(iax+i, jax+j)
    enddo
enddo

call exact_biquad_ef(xub(1:d), yub(1:d), bub(1:d), pos_xpoint(1), pos_xpoint(2), f00, ddpsi, dr, dz)

return
end subroutine nine_point_regression

!---------------------------------------------------------------------
subroutine nine_point_regression_follow(rx, zx, pos_xpoint, ddpsi, f00)

use ef_circuit, only: dr, dz

implicit none

integer, parameter :: ndim=9

double precision, intent(in) :: rx, zx
double precision, intent(out) :: f00
double precision, intent(out), dimension(2) :: pos_xpoint
double precision, intent(out), dimension(ndim-1) :: ddpsi

double precision, dimension(ndim) :: xub, bub, yub

call find_fields_interp_ef_psionly(rx - dr, zx - dz, bub(1))
call find_fields_interp_ef_psionly(rx, zx - dz, bub(2))
call find_fields_interp_ef_psionly(rx + dr, zx - dz, bub(3))
call find_fields_interp_ef_psionly(rx - dr, zx, bub(4))
call find_fields_interp_ef_psionly(rx, zx, bub(5))
call find_fields_interp_ef_psionly(rx + dr, zx, bub(6))
call find_fields_interp_ef_psionly(rx - dr, zx + dz, bub(7))
call find_fields_interp_ef_psionly(rx, zx + dz, bub(8))
call find_fields_interp_ef_psionly(rx + dr, zx + dz, bub(9))

xub(1) = rx - dr
yub(1) = zx - dz
xub(2) = rx
yub(2) = zx - dz
xub(3) = rx + dr
yub(3) = zx - dz
xub(4) = rx - dr
yub(4) = zx
xub(5) = rx
yub(5) = zx
xub(6) = rx + dr
yub(6) = zx
xub(7) = rx - dr
yub(7) = zx + dz
xub(8) = rx
yub(8) = zx + dz
xub(9) = rx + dr
yub(9) = zx + dz

call exact_biquad_regress_ef(xub, yub, bub, pos_xpoint(1), pos_xpoint(2), f00, ddpsi, dr, dz, rx, zx)

return
end subroutine nine_point_regression_follow

!---------------------------------------------------------------------
subroutine find_actual_index_ef(r0, z0, i, j)

use ef_circuit, only: rmin, dr, zmin, dz

implicit none

double precision, intent(in) :: r0, z0
integer, intent(out) :: i, j

i = nint((r0 - rmin)/dr + 1.)
j = nint((z0 - zmin)/dz + 1.)

return
end subroutine find_actual_index_ef

!---------------------------------------------------------------------
subroutine find_fields_interp_ef_psionly(r0, z0, psi0) !give back psi, br, bz at r0, z0

use ef_circuit, only: rmin, r, dr, zmin, z, dz, psirz

implicit none
double precision, intent(in)  :: r0, z0
double precision, intent(out) :: psi0
integer i, j
double precision x1, x2, y1, y2
double precision z1, z2, z3, z4

i = floor((r0 -rmin)/dr + 1.)
j = floor((z0 -zmin)/dz + 1.)
x1 = r(i)
x2 = r(i+1)
y1 = z(j)
y2 = z(j+1)
z1 = psirz(i, j)
z2 = psirz(i+1, j)
z3 = psirz(i, j+1)
z4 = psirz(i+1, j+1)

!bilinear interpolation
call bilinear_average_ef(x1, x2, y1, y2, r0, z0, z1, z2, z3, z4, psi0)

return
end subroutine find_fields_interp_ef_psionly

!---------------------------------------------------------------------
subroutine find_fields_interp_ef_green(r0, z0, psi0, iconduc)

use ef_circuit, only: rmin, r, dr, zmin, z, dz
use green_matrix, only: greeni

implicit none

integer, intent(in) :: iconduc
double precision, intent(in) :: r0, z0
double precision, intent(out) :: psi0

integer :: i, j
double precision :: x1, x2, y1, y2, z1, z2, z3, z4

i = floor((r0 - rmin)/dr + 1.)
j = floor((z0 - zmin)/dz + 1.)
x1 = r(i)
x2 = r(i+1)
y1 = z(j)
y2 = z(j+1)
z1 = greeni(i, j, iconduc)
z2 = greeni(i + 1, j, iconduc)
z3 = greeni(i, j + 1, iconduc)
z4 = greeni(i + 1, j + 1, iconduc)

!bilinear interpolation
call bilinear_average_ef(x1, x2, y1, y2, r0, z0, z1, z2, z3, z4, psi0)

return
end subroutine find_fields_interp_ef_green

!---------------------------------------------------------------------
subroutine bilinear_average_ef(x1, x2, y1, y2, x, y, f11, f21, f12, f22, f_out)

implicit none

double precision, intent(in)  :: x1, x2, y1, y2, x, y, f11, f21, f12, f22
double precision, intent(out) :: f_out

f_out = 1./((x2 - x1)*(y2 - y1)) * &
    ( f11*(x2 - x)*(y2 - y) + f21*(x - x1)*(y2 - y) + & 
      f12*(x2 - x)*(y - y1) + f22*(x - x1)*(y - y1) )

return
end subroutine bilinear_average_ef

!---------------------------------------------------------------------
function ellE_green(X, DL)
! First kind elliptic integral, from K. Lackner, T. Lunt, IPP - Garching 2022

implicit none

double precision, intent(in) :: X, DL
double precision :: ellE_green

ellE_green = (((0.01736506451D0 *X + 0.04757383546D0)*X + 0.06260601220D0)*X + 0.44325141463D0)*X + 1.0D0 - &
             (((0.00526449639D0 *X + 0.04069697526D0)*X + 0.09200180037D0)*X + 0.24998368310D0)*X*DL

return
end function ellE_green

!---------------------------------------------------------------------
function ellK_green(X, DL)
! Second kind elliptic integral, from K. Lackner, T. Lunt, IPP - Garching 2022

implicit none

double precision, intent(in) :: X, DL
double precision :: ellK_green

ellK_green = ((( 0.01451196212D0*X + 0.03742563713D0)*X + 0.03590092383D0)*X + 0.09666344259D0)*X + 1.38629436112D0 - &
             ((((0.00441787012D0*X + 0.03328355346D0)*X + 0.06880248576D0)*X + 0.12498593597D0)*X + 0.5D0)*DL

return
end function ellK_green

!---------------------------------------------------------------------
subroutine green_function(r1, z1, r2, z2, greenf)

implicit none

double precision, intent(in) :: r1, z1, r2, z2
double precision, intent(out) ::greenf
double precision :: TT, K, ELCK, ELCE, ellk_green, elle_green, acl, alg

K = sqrt(4.*r1*r2/((r2 + r1)**2 + (z2 - z1)**2) )

TT = 1. - K**2

acl = tt
alg = dlog(acl)

ELCK = ellK_green(acl, alg)
ELCE = ellE_green(acl, alg)

greenf = ( (1.D0 - K**2/2.)*ELCK - ELCE )*(SQRT(r1*r2)/K )

return
end subroutine green_function

!---------------------------------------------------------------------
subroutine green_function_identity(r1, z1, greenf, dr, dz)

implicit none

double precision, intent(in)  :: r1, z1, dr, dz
double precision, intent(out) :: greenf

greenf = r1*(log(8.*r1/(0.2236*(dr + dz))) - 2.0) ! From A. Kavin, used in spider, self inductance of a rectangular coil in toroidal direction

return
end subroutine green_function_identity

!---------------------------------------------------------------------
subroutine boundary_ef(g)

! new bc is integral_over_boundary of -Green * dg/dn * dl
use dimensions_ef_parameters, only: i_dim2
use pi_grec_vars, only: GPI
use ef_circuit, only: nr1, nz1, nr2, nz2

implicit none

double precision, intent(inout), dimension(i_dim2, i_dim2) :: g

integer :: i, jcounty
double precision, dimension(i_dim2, 4) :: integr

integr = 0.
jcounty = 0
! lower side
do i=2, nr1
    call bgint_ef(integr(i, 1), g, jcounty)
enddo
! right side
do i=2, nz1
    call bgint_ef(integr(i, 2), g, jcounty)
enddo
! upper side
do i=2, nr1
    call bgint_ef(integr(i, 3), g, jcounty)
enddo
! left side
do i=2, nz1
    call bgint_ef(integr(i, 4), g, jcounty)
enddo
g(2:nr1,   1) = integr(2:nr1, 1)/GPI
g(nr2, 2:nz1) = integr(2:nz1, 2)/GPI
g(2:nr1, nz2) = integr(2:nr1, 3)/GPI
g(1,   2:nz1) = integr(2:nz1, 4)/GPI

return
end subroutine boundary_ef

!---------------------------------------------------------------------
subroutine bgint_ef(bgintsol, g, jcounty)

! calculates  integral_over_boundary of -Green * dg/dn * dl for point r0, z0
use dimensions_ef_parameters, only: i_dim2
use ef_circuit, only: nr1, nz1, nr2, nz2, dr, r, dz, green_bnd_f

implicit none

double precision, intent(in), dimension(i_dim2, i_dim2) :: g(i_dim2, i_dim2)
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
end subroutine bgint_ef
