subroutine find_angle_ef(rt, zt, r, z, angle)

use pi_grec_vars, only: GPI2
use debugger_ef, only: markloc_ef, debug

double precision, intent(in) :: rt, zt, r, z
double precision, intent(out) :: angle

angle = ATAN2(z-zt, r-rt)
if (angle < 0) angle = angle + GPI2

return
end subroutine find_angle_ef

!------------------------------------------------------------
subroutine qinterp_ef_feqis(x1, y1, Nx1, x2, y2, Nx2)

use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: Nx1, Nx2
double precision, intent(in) , dimension(Nx1) :: x1, y1
double precision, intent(in) , dimension(Nx2) :: x2
double precision, intent(out), dimension(Nx2) :: y2

integer :: i, j
double precision :: A, B, C, z1, z2, z3, t1, t2, t3, t4

call markloc_ef('qinterp_ef_feqis', debug_lev=debug-1)

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
            C = t2 - A*z2**2 - B*z2
            y2(i) = A*t4**4 + B*t4**2 + C
            EXIT
        elseif (t4 > z3 .and. j == Nx1-1) then
            t1 = y1(j-1)
            t2 = y1(j)
            t3 = y1(j+1)
            A = (t3 - t2 - (z3 - z2)*(t1 - t2)/(z1 - z2)) / ((z3 - z2)*(z3 - z1))
            B = (t1 - t2)/(z1 - z2) - A*(z1 + z2)
            C = t2 - A*z2**2 - B*z2
            y2(i) = A*t4**2 + B*t4 + C
            EXIT
        endif
    enddo
enddo

return
end subroutine qinterp_ef_feqis

!------------------------------------------------------------
subroutine linterp_ef_feqis(x1, y1, Nx1, x2, y2, Nx2)

use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: Nx1, Nx2
double precision, intent(in), dimension(Nx1) :: x1, y1
double precision, intent(in), dimension(Nx2) :: x2
double precision, intent(out), dimension(Nx2) :: y2

integer :: i, j
double precision :: A, C, z1, z2, t1, t2, t4

call markloc_ef('linterp_ef_feqis', debug_lev=debug-1)

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
end subroutine linterp_ef_feqis

!------------------------------------------------------------
subroutine polyfitcc_feqis(x, y, P)

use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in) , dimension(3) :: x, y
double precision, intent(out), dimension(3) :: P

double precision :: y21, y32, x21, x32, h21, h32

call markloc_ef('polyfitcc_feqis', debug_lev=debug)

y32=y(3)-y(2)                 
y21=y(2)-y(1)                 
x32=x(3)-x(2)                 
x21=x(2)-x(1)                 
h32=x(3)+x(2)                 
h21=x(2)+x(1)                 

P(1)=(x21*y32-x32*y21)/(x21*x32*(h32-h21))
P(2)=y21/x21-P(1)*h21
P(3)=y(3)-P(1)*x(3)**2.-P(2)*x(3)

return
end subroutine polyfitcc_feqis

!------------------------------------------------------------
subroutine inverse_matrix_feqis(a, c, n)
! a(n, n) - array of coefficients for matrix A
! n      - dimension
! c(n, n) - inverse matrix of A

use debugger_ef, only: markloc_ef, debug

implicit none 

integer, intent(in) :: n
double precision, intent(out)  , dimension(n, n) :: c
double precision, intent(inout), dimension(n, n) :: a

integer :: i, j, k
double precision :: coeff
double precision, dimension(n) :: b, d, x
double precision, dimension(n, n) :: L, U

call markloc_ef('inverse_matrix_feqis', debug_lev=debug)

L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1,  n-1
    do i=k+1, n
        coeff=a(i, k)/a(k, k)
        L(i, k) = coeff
        do j=k+1, n
            a(i, j) = a(i, j)-coeff*a(k, j)
        enddo
    enddo
enddo

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1, n
    L(i, i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1, n
    do i=1, j
        U(i, j) = a(i, j)
    enddo
enddo

! Step 3: compute columns of the inverse matrix C
do k=1, n
    b(k)=1.0
    d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
    do i=2, n
        d(i)=b(i)
        do j=1, i-1
            d(i) = d(i) - L(i, j)*d(j)
        enddo
    enddo
! Step 3b: Solve Ux=d using the back substitution
    x(n)=d(n)/U(n, n)
    do i = n-1, 1, -1
        x(i) = d(i)
        do j=n, i+1, -1
            x(i)=x(i)-U(i, j)*x(j)
        enddo
        x(i) = x(i)/u(i, i)
    enddo
! Step 3c: fill the solutions x(n) into column k of C
    do i=1, n
        c(i, k) = x(i)
    enddo
    b(k)=0.0
enddo

return
end subroutine inverse_matrix_feqis

!------------------------------------------------------------
subroutine interp_j_fromrhotorz

use ef_circuit, only: nr2, nz2, r, z, jrhoteta, nrho, nteta, &
    rho, teta, raxp, zaxp, jrz
use debugger_ef, only: markloc_ef, debug

implicit none

integer :: i, j

call markloc_ef('interp_j_fromrhotorz', debug_lev=debug)

! go from jrhoteta to jrz
jrz=0.
do j=1, nz2
    do i=1, nr2
        call curinterp_ef(r(i), z(j), jrhoteta(1: nrho, 1: nteta+1),  & 
            rho(1: nrho, 1: nteta+1), teta(1: nteta+1), raxp, zaxp, nrho, nteta+1, jrz(i, j))
    enddo
enddo

return
end subroutine interp_j_fromrhotorz

!------------------------------------------------------------
subroutine curinterp_ef(r, z, jrho, rho, teta, rax, zax, nrho, nteta, j)

use pi_grec_vars, only: GPI2
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: nrho, nteta
double precision, intent(in) :: r, z, rax, zax
double precision, intent(in), dimension(nteta) :: teta
double precision, intent(in), dimension(nrho, nteta) :: jrho, rho
double precision, intent(out) :: j

integer :: i, j1, j2, k1, k2, k3, k4
double precision :: anglr, rho0
double precision, dimension(1) :: jarr
double precision, dimension(2) :: x1, y1, jt

call markloc_ef('curinterp_ef', debug_lev=debug)

call find_angle_ef(rax, zax, r, z, anglr)
rho0=sqrt((r-rax)**2.0+(z-zax)**2.0)
  
if (anglr < teta(1)) anglr=anglr+GPI2
do i=1, nteta-1
    if (anglr.ge.teta(i)) j1=i
enddo
j2=j1+1

if (rho0 > rho(nrho, j1) .or. rho0 > rho(nrho, j2)) then
    j=0.
    return
endif

do i=1, nrho-1
    if (rho0 > rho(i, j1)) k1=i
enddo    
    k2=k1+1
do i=1, nrho-1
    if (rho0 > rho(i, j2)) k3=i
enddo    
k4=k3+1

x1(1)=rho(k1, j1)
x1(2)=rho(k2, j1)
y1(1)=jrho(k1, j1)
y1(2)=jrho(k2, j1)
call linterp_ef_feqis(x1, y1, 2, (/rho0/), jt(1), 1)

x1(1)=rho(k3, j2)
x1(2)=rho(k4, j2)
y1(1)=jrho(k3, j2)
y1(2)=jrho(k4, j2)
call linterp_ef_feqis(x1, y1, 2, (/rho0/), jt(2), 1)

x1(1)=teta(j1)
x1(2)=teta(j2)
y1(1)=jt(1)
y1(2)=jt(2)
call linterp_ef_feqis(x1, y1, 2, (/anglr/), jarr, 1)
j = jarr(1)

return
end subroutine curinterp_ef

!------------------------------------------------------------
subroutine discrete_sine_transform_ef(n, y_in, y_out)

use pi_grec_vars, only: GPI
use debugger_ef, only: markloc_ef, debug

implicit none

integer :: n
double precision,  intent(in) , dimension(n) :: y_in
double precision,  intent(out), dimension(n) :: y_out

integer :: i, j

call markloc_ef('discrete_sine_transform', debug_lev=debug)

y_out=0.
do i=1, n
    do j=1, n
        y_out(i)=y_out(i)+y_in(j)*sin(i*j*GPI/(n+1))
    enddo
enddo

return
end subroutine discrete_sine_transform_ef

!------------------------------------------------------------
subroutine plasma_psi_to_coils_ef

use ef_circuit, only: nr2, nz2, nconduc, jrz, area_eff, dr, dz, &
    psiplasmatoconduc
use green_matrix, only: greeni
use debugger_ef, only: markloc_ef, debug

implicit none

integer :: i

call markloc_ef('plasma_psi_to_coils_ef', debug_lev=debug)

do i=1, nconduc
    psiplasmatoconduc(i)=sum(jrz(1: nr2, 1: nz2)*area_eff(1: nr2, 1: nz2)* &
        greeni(1: nr2, 1: nz2, i))*dr*dz
enddo

return
end subroutine plasma_psi_to_coils_ef

!------------------------------------------------------------
subroutine psiplex_calc_ef

use ef_circuit, only: psplex
use debugger_ef, only: markloc_ef, debug

implicit none

call markloc_ef('psiplex_calc_ef', debug_lev=0)

psplex=2.

return
end subroutine psiplex_calc_ef

!------------------------------------------------------------
subroutine psi_external_calc_ef

use ef_circuit, only: nr2, nz2, nconduc, curconduc, psiextrz
use green_matrix, only: greeni
use debugger_ef, only: markloc_ef, debug

implicit none

integer :: i, j

call markloc_ef('psi_external_calc_ef', debug_lev=debug)

do j=1, nz2    
    do i=1, nr2    
        psiextrz(i, j)=sum(curconduc(1: nconduc)*greeni(i, j, 1: nconduc))
    enddo
enddo

return
end subroutine psi_external_calc_ef

!------------------------------------------------------------
subroutine generate_zlim_potential(r, z, nt, rl, zl, zlimp)

use exchange_with_astra, only: use_zlim_pot
use pi_grec_vars, only: GPI2
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: nt
double precision, intent(in) :: r, z
double precision, intent(in), dimension(nt) :: rl, zl
double precision, intent(out) :: zlimp

integer :: j
double precision :: discrim, rt, zt
double precision, dimension(nt) :: anglr, anglf, anglt, diffa

call markloc_ef('generate_zlim_potential', debug_lev=debug)

zlimp=-1.

rt=r
zt=z

!check if sol current filament is inside plasma
do j=1, nt
    call find_angle_ef(rt, zt, rl(j), zl(j), anglr(j))
enddo

do j=1, nt
    if (j > 1) then
        diffa(j)=abs(anglr(j)-anglr(j-1))
        if (diffa(j) > 1.) diffa(j)=abs(diffa(j)-GPI2)
    endif
    if (j == 1) then
        diffa(j)=abs(anglr(j)-anglr(nt))
        if (diffa(j) > 1.) diffa(j)=abs(diffa(j)-GPI2)
    endif
enddo

anglt(1: nt)=sin(anglr(1: nt))*diffa(1: nt)
anglf(1: nt)=cos(anglr(1: nt))*diffa(1: nt)    

discrim=abs(sum(anglt(1: nt)))/sum(diffa(1: nt))+ &
        abs(sum(anglf(1: nt)))/sum(diffa(1: nt))

if (discrim < 0.6) zlimp=1.
if (use_zlim_pot == 0) zlimp=1.

return
end subroutine generate_zlim_potential

!------------------------------------------------------------
subroutine least_square_biquad_ef(r, z, u, n, c, rax, zax, uax, derivs)

use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: n
double precision, intent(in), dimension(n) :: r, z, u
double precision, intent(out) :: rax, zax, uax, derivs(8)

integer :: k
double precision :: s_r, s_z, s_r2, s_z2, s_rz, s_r3, s_rz2, s_r2z, s_z3, &
    s_r4, s_r2z2, s_r3z, s_z4, s_rz3, s_u, s_ur, s_uz, s_ur2, s_uz2, s_urz, &  
    det, det_r, det_z
double precision, dimension(6) :: B, c
double precision, dimension(6, 6) :: A, Ainv

call markloc_ef('least_square_biquad_ef', debug_lev=debug)

s_r= 0.d0
s_z= 0.d0
s_r2=0.d0
s_z2=0.d0
s_rz=0.d0
s_r3= 0.d0
s_rz2=0.d0
s_r2z=0.d0
s_z3=0.d0
s_r4=  0.d0
s_r2z2=0.d0
s_r3z= 0.d0
s_z4= 0.d0
s_rz3=0.d0
s_u= 0.d0
s_ur=0.d0
s_uz=0.d0
s_ur2= 0.d0
s_uz2=0.d0
s_urz=0.d0

do k=1, n
    s_r=  s_r + (r(k)-r(1))
    s_z=  s_z + (z(k)-z(1))
    s_r2= s_r2 + (r(k)-r(1))**2
    s_z2= s_z2 + (z(k)-z(1))**2
    s_rz= s_rz + (r(k)-r(1))*(z(k)-z(1))
    s_u=  s_u +u(k)
    s_r3=  s_r3 + (r(k)-r(1))**3
    s_rz2= s_rz2 + (r(k)-r(1))*(z(k)-z(1))**2
    s_r2z= s_r2z + (r(k)-r(1))**2*(z(k)-z(1))
    s_ur=  s_ur + u(k)*(r(k)-r(1))
    s_z3= s_z3 + (z(k)-z(1))**3
    s_uz= s_uz + u(k)*(z(k)-z(1))
    s_r4=   s_r4 + (r(k)-r(1))**4
    s_r2z2= s_r2z2 + (r(k)-r(1))**2*(z(k)-z(1))**2
    s_r3z=  s_r3z +  (r(k)-r(1))**3*(z(k)-z(1))
    s_ur2=  s_ur2 +  u(k)*(r(k)-r(1))**2
    s_z4=  s_z4 + (z(k)-z(1))**4
    s_rz3= s_rz3 + (r(k)-r(1))*(z(k)-z(1))**3
    s_uz2= s_uz2 + u(k)*(z(k)-z(1))**2 
    s_urz= s_urz + u(k)*(r(k)-r(1))*(z(k)-z(1))
enddo

! Creating matrix A and right-hand B:

A(1, 1) = n+0.
A(1, 2) = s_r
A(1, 3) = s_z
A(1, 4) = s_r2
A(1, 5) = s_z2
A(1, 6) = s_rz
B(1) = s_u

A(2, 1) = s_r
A(2, 2) = s_r2
A(2, 3) = s_rz
A(2, 4) = s_r3
A(2, 5) = s_rz2
A(2, 6) = s_r2z
B(2) = s_ur

A(3, 1) = s_z
A(3, 2) = s_rz
A(3, 3) = s_z2
A(3, 4) = s_r2z
A(3, 5) = s_z3
A(3, 6) = s_rz2
B(3) = s_uz

A(4, 1) = s_r2
A(4, 2) = s_r3
A(4, 3) = s_r2z
A(4, 4) = s_r4
A(4, 5) = s_r2z2
A(4, 6) = s_r3z
B(4) = s_ur2

A(5, 1) = s_z2
A(5, 2) = s_rz2
A(5, 3) = s_z3
A(5, 4) = s_r2z2
A(5, 5) = s_z4
A(5, 6) = s_rz3
B(5) = s_uz2

A(6, 1) = s_rz
A(6, 2) = s_r2z
A(6, 3) = s_rz2
A(6, 4) = s_r3z
A(6, 5) = s_rz3
A(6, 6) = s_r2z2
B(6) = s_urz

! Find coefficients
call inverse_matrix_feqis(A, Ainv, 6)
do k=1, 6
    c(k)=sum(Ainv(k, 1: 6)*B(1: 6))
enddo

! Magnetic axis

det=4.d0*c(4)*c(5)-c(6)**2

det_r=-2.d0*c(2)*c(5)+c(6)*c(3)
det_z=-2.d0*c(3)*c(4)+c(6)*c(2)

rax=det_r/det
zax=det_z/det

uax=c(1)+c(2)*rax+c(3)*zax+c(4)*rax**2+c(5)*zax**2+c(6)*rax*zax

rax=rax+r(1)
zax=zax+z(1)

derivs(1)=c(2)
derivs(2)=c(3)
derivs(3)=2.*c(4)
derivs(4)=2.*c(5)
derivs(5)=c(6) 

return
end subroutine least_square_biquad_ef
             
!------------------------------------------------------------
double precision function frlim_ef(dp, ylim, rx, zx, rm, zm)

use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in) :: dp(5), ylim, rx, zx, rm, zm
double precision :: dxx, dxy, dyy, disc, cc, cdpls, cdmns, &
    ang1, ang2, c1, dl2x, c2, dl2y

call markloc_ef('frlim_ef', debug_lev=debug)

Dxx=dp(3)
Dxy=dp(4)
Dyy=dp(5)

disc=(Dxy/Dyy)**2 - Dxx/Dyy

if(Disc < 0.) then
    cc=(zm-zx)/(rm-rx)
    cc=-1.d0/cc
else
    cdpls = -Dxy/Dyy + dsqrt(disc)
    cdmns = -Dxy/Dyy - dsqrt(disc)

    ang1 = 0.5d0*(datan(cdpls)+datan(cdmns))
    ang2 =-0.5d0*(datan(1.d0/cdpls)+datan(1.d0/cdmns))

! D2u/Dl2(direction ang1)

    c1=dtan(ang1)
    Dl2x=Dyy*c1*c1+2.d0*Dxy*c1+Dxx

! D2u/Dl2(direction ang2)

    c2=dtan(ang2)
    Dl2y=Dyy*c2*c2+2.*Dxy*c2+Dxx

    if(Dl2x < 0.) then
        cc=c1
    elseif(Dl2y < 0.) then
        cc=c2
    else
        cc=(zm-zx)/(rm-rx)
        cc=-1.d0/cc
    endif
endif

frlim_ef=rx+(ylim-zx)/cc

return
end function frlim_ef

!------------------------------------------------------------
subroutine nine_point_regression(r0, z0, pos_xpoint, ddpsi, f00)

use ef_circuit, only: r, z, psirz
use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in) :: r0, z0
double precision, intent(out) :: ddpsi(8), f00, pos_xpoint(2)
integer :: iax, jax
double precision :: xub(9), bub(9), yub(9), c(6)

call markloc_ef('nine_point_regression', debug_lev=debug)

call find_actual_index_ef(r0, z0, iax, jax)
xub(1)=r(iax-1)
xub(2)=r(iax)
xub(3)=r(iax+1)
xub(4)=r(iax)
xub(5)=r(iax)
xub(6)=r(iax-1)
xub(7)=r(iax-1)
xub(8)=r(iax+1)
xub(9)=r(iax+1)
yub(1)=z(jax)
yub(2)=z(jax)
yub(3)=z(jax)
yub(4)=z(jax-1)
yub(5)=z(jax+1)
yub(6)=z(jax-1)
yub(7)=z(jax+1)
yub(8)=z(jax-1)
yub(9)=z(jax+1)
bub(1)=psirz(iax-1, jax)
bub(2)=psirz(iax, jax)
bub(3)=psirz(iax+1, jax)
bub(4)=psirz(iax, jax-1)
bub(5)=psirz(iax, jax+1)
bub(6)=psirz(iax-1, jax-1)
bub(7)=psirz(iax-1, jax+1)
bub(8)=psirz(iax+1, jax-1)
bub(9)=psirz(iax+1, jax+1)

call least_square_biquad_ef(xub, yub, bub, 9, c, pos_xpoint(1), pos_xpoint(2), f00, ddpsi)

return
end subroutine nine_point_regression

!------------------------------------------------------------
subroutine find_actual_index_ef(r0, z0, i, j)

use ef_circuit, only: rmin, zmin, dr, dz
use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in) :: r0, z0
integer, intent(out) :: i, j

call markloc_ef('find_actual_index_ef', debug_lev=debug)

i=nint((r0-rmin)/dr+1.)    
j=nint((z0-zmin)/dz+1.)    

return
end subroutine find_actual_index_ef

!------------------------------------------------------------
subroutine find_floor_index_ef(r0, z0, i, j)

use ef_circuit, only: rmin, zmin, dr, dz
use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in) :: r0, z0
integer, intent(out) :: i, j

call markloc_ef('find_floor_index_ef', debug_lev=debug)

i=floor((r0-rmin)/dr+1.)
j=floor((z0-zmin)/dz+1.)

return
end subroutine find_floor_index_ef

!------------------------------------------------------------
subroutine find_fields_interp_ef(r0, z0, psi0, br0, bz0, brr, brz, bzr, bzz)

use ef_circuit, only: dr, dz, r, z, psirz
use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in) :: r0, z0
double precision, intent(out) :: psi0, br0, bz0, brr, brz, bzr, bzz
integer :: i, j
double precision :: x1, x2, y1, y2, z1, z2, z3, z4, t1, t2, t3, t4

call markloc_ef('find_fields_interp_ef', debug_lev=debug)

call find_floor_index_ef(r0,z0,i,j)

x1=r(i)
x2=r(i+1)
y1=z(j)
y2=z(j+1)
z1=psirz(i, j)
z2=psirz(i+1, j)
z3=psirz(i, j+1)
z4=psirz(i+1, j+1)

!bilinear interpolation
call bilinear_average_ef(x1, x2, y1, y2, r0, z0, z1, z2, z3, z4, psi0)

!br
t1=-1./x1*(psirz(i, j+1)-psirz(i, j-1))/dz/2.
t2=-1./x2*(psirz(i+1, j+1)-psirz(i+1, j-1))/dz/2.
t3=-1./x1*(psirz(i, j+2)-psirz(i, j))/dz/2.
t4=-1./x2*(psirz(i+1, j+2)-psirz(i+1, j))/dz/2.
call bilinear_average_ef(x1, x2, y1, y2, r0, z0, t1, t2, t3, t4, br0)
br0=-br0

!bz
t1=1./x1*(psirz(i+1, j)-psirz(i-1, j))/dr/2.
t2=1./x2*(psirz(i+2, j)-psirz(i, j))/dr/2.
t3=1./x1*(psirz(i+1, j+1)-psirz(i-1, j+1))/dr/2.
t4=1./x2*(psirz(i+2, j+1)-psirz(i, j+1))/dr/2.
call bilinear_average_ef(x1, x2, y1, y2, r0, z0, t1, t2, t3, t4, bz0)
bz0=-bz0

!brr
t1=(-1./r(i+1)*(psirz(i+1, j+1)-psirz(i+1, j-1))/dz/2.+1./r(i-1)*(psirz(i-1, j+1)-psirz(i-1, j-1))/dz/2.)/dr/2.
t2=(-1./r(i+2)*(psirz(i+2, j+1)-psirz(i+2, j-1))/dz/2.+1./r(i)*(psirz(i, j+1)-psirz(i, j-1))/dz/2.)/dr/2.
t3=(-1./r(i+1)*(psirz(i+1, j+2)-psirz(i+1, j))/dz/2.+1./r(i-1)*(psirz(i-1, j+2)-psirz(i-1, j))/dz/2.)/dr/2.
t4=(-1./r(i+2)*(psirz(i+2, j+2)-psirz(i+2, j))/dz/2.+1./r(i)*(psirz(i, j+2)-psirz(i, j))/dz/2.)/dr/2.
call bilinear_average_ef(x1, x2, y1, y2, r0, z0, t1, t2, t3, t4, brr)
brr=-brr

!brz
t1=(-1./r(i)*(psirz(i, j+1)-psirz(i, j))/dz+1./r(i)*(psirz(i, j)-psirz(i, j-1))/dz)/dz
t2=(-1./r(i+1)*(psirz(i+1, j+1)-psirz(i+1, j))/dz+1./r(i+1)*(psirz(i+1, j)-psirz(i+1, j-1))/dz)/dz
t3=(-1./r(i)*(psirz(i, j+2)-psirz(i, j+1))/dz+1./r(i)*(psirz(i, j+1)-psirz(i, j))/dz)/dz
t4=(-1./r(i+1)*(psirz(i+1, j+2)-psirz(i+1, j+1))/dz+1./r(i+1)*(psirz(i+1, j+1)-psirz(i+1, j))/dz)/dz
call bilinear_average_ef(x1, x2, y1, y2, r0, z0, t1, t2, t3, t4, brz)
brz=-brz

!bzr
t1=(+2./(r(i)+r(i+1))*(psirz(i+1, j)-psirz(i, j))/dr-2./(r(i)+r(i-1))*(psirz(i, j)-psirz(i-1, j))/dr)/dr
t2=(+2./(r(i+2)+r(i+1))*(psirz(i+2, j)-psirz(i+1, j))/dr-2./(r(i)+r(i+1))*(psirz(i+1, j)-psirz(i, j))/dr)/dr
t3=(+2./(r(i)+r(i+1))*(psirz(i+1, j+1)-psirz(i, j+1))/dr-2./(r(i)+r(i-1))*(psirz(i, j+1)-psirz(i, j+1))/dr)/dr
t4=(+2./(r(i+2)+r(i+1))*(psirz(i+2, j+1)-psirz(i+1, j+1))/dr-2./(r(i)+r(i+1))*(psirz(i+1, j+1)-psirz(i, j+1))/dr)/dr
call bilinear_average_ef(x1, x2, y1, y2, r0, z0, t1, t2, t3, t4, bzr)
bzr=-bzr

!bzz
t1=(+1./(r(i))*(psirz(i+1, j+1)-psirz(i-1, j+1))/dr/2.-1./(r(i))*(psirz(i+1, j-1)-psirz(i-1, j-1))/dr/2.)/dz/2.
t2=(+1./(r(i+1))*(psirz(i+2, j+1)-psirz(i, j+1))/dr/2.-1./(r(i+1))*(psirz(i+2, j-1)-psirz(i, j-1))/dr/2.)/dz/2.
t3=(+1./(r(i))*(psirz(i+1, j+2)-psirz(i-1, j+2))/dr/2.-1./(r(i))*(psirz(i+1, j)-psirz(i-1, j))/dr/2.)/dz/2.
t4=(+1./(r(i+1))*(psirz(i+2, j+2)-psirz(i, j+2))/dr/2.-1./(r(i+1))*(psirz(i+2, j)-psirz(i, j))/dr/2.)/dz/2.
call bilinear_average_ef(x1, x2, y1, y2, r0, z0, t1, t2, t3, t4, bzz)
bzz=-bzz

return
end subroutine find_fields_interp_ef

!------------------------------------------------------------
subroutine find_fields_interp_ef_psionly(r0, z0, psi0) !give back psi, br, bz at r0, z0

use ef_circuit, only: r, z, psirz
use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in)  :: r0, z0
double precision, intent(out) :: psi0

integer :: i, j
double precision :: x1, x2, y1, y2, z1, z2, z3, z4

call markloc_ef('find_fields_interp_ef_psionly', debug_lev=debug)

call find_floor_index_ef(r0,z0,i,j)

x1=r(i)
x2=r(i+1)
y1=z(j)
y2=z(j+1)
z1=psirz(i, j)
z2=psirz(i+1, j)
z3=psirz(i, j+1)
z4=psirz(i+1, j+1)

!bilinear interpolation
call bilinear_average_ef(x1, x2, y1, y2, r0, z0, z1, z2, z3, z4, psi0)

return
end subroutine find_fields_interp_ef_psionly

!------------------------------------------------------------
subroutine find_fields_interp_ef_psionly_neg(r0, z0, psi0) !give back psi, br, bz at r0, z0

use ef_circuit, only: r, z, psirz
use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in)  :: r0, z0
double precision, intent(out) :: psi0

integer :: i, j
double precision :: x1, x2, y1, y2, z1, z2, z3, z4

call markloc_ef('find_fields_interp_ef_psionly_neg', debug_lev=debug)

call find_floor_index_ef(r0,z0,i,j)

x1=r(i)
x2=r(i+1)
y1=z(j)
y2=z(j+1)
z1=-psirz(i, j)
z2=-psirz(i+1, j)
z3=-psirz(i, j+1)
z4=-psirz(i+1, j+1)

!bilinear interpolation
call bilinear_average_ef(x1, x2, y1, y2, r0, z0, z1, z2, z3, z4, psi0)

return
end subroutine find_fields_interp_ef_psionly_neg

!------------------------------------------------------------
subroutine bilinear_average_ef(x1, x2, y1, y2, x, y, f11, f21, f12, f22, f0)

use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in)  :: x1, x2, y1, y2, x, y, f11, f21, f12, f22
double precision, intent(out) :: f0

call markloc_ef('bilinear_average_ef', debug_lev=debug)

f0 = 1./((x2 - x1)*(y2 - y1)) * ( &
     f11*(x2 - x )*(y2 - y) + f21*(x - x1)*(y2 - y) + & 
     f12*(x2 - x )*(y - y1) + f22*(x - x1)*(y - y1) )

return
end subroutine bilinear_average_ef

!------------------------------------------------------------
subroutine green_function(r1, z1, r2, z2, greenf)

use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in)  :: r1, z1, r2, z2
double precision, intent(out) :: greenf

integer :: ifailk, ifaile
double precision :: TT, K, ELCK, ELCE, s21bbf, s21bcf

call markloc_ef('green_function', debug_lev=debug)

K = sqrt(4.*r1*r2/((r2 + r1)**2 + (z2 - z1)**2))
TT = 1. - K**2.
ELCK = S21BBF(0.D0, TT, 1.D0, IFAILK)
ELCE = ELCK - K**2./3.D0*S21BCF(0.D0, TT, 1.D0, IFAILE)
greenf = ( (1.D0 - 0.5*K**2)*ELCK - ELCE )*( SQRT(r1*r2)/K )

return
end subroutine green_function

!------------------------------------------------------------
subroutine green_function_identity(r1, greenf, dr, dz)
! From A. Kavin, used in spider: self inductance of a rectangular
! coil in toroidal direction

use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in)  :: r1, dr, dz
double precision, intent(out) :: greenf

call markloc_ef('green_function_identity', debug_lev=debug)

greenf = r1*(log(8.*r1/(0.2236*(dr + dz))) - 2.0) 

return
end subroutine green_function_identity

!------------------------------------------------------------
subroutine boundary_ef(g)

! new bc is integral_over_boundary of -Green * dg/dn * dl
use pi_grec_vars,  only: GPI
use ef_circuit, only: i_dim2, nr2, nz2
use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(out), dimension(i_dim2, i_dim2) :: g
integer :: i, jcounty
double precision :: integr(i_dim2, 4)

call markloc_ef('boundary_ef', debug_lev=debug)

integr=0.
jcounty=0
! lower side
do i=1, nr2
    call bgint_ef(integr(i, 1), g, jcounty)
enddo
! right side
do i=1, nz2
    call bgint_ef(integr(i, 2), g, jcounty)
enddo
! upper side
do i=1, nr2
    call bgint_ef(integr(i, 3), g, jcounty)
enddo
! left side
do i=1, nz2
    call bgint_ef(integr(i, 4), g, jcounty)
enddo
g(1: nr2, 1)   = integr(1: nr2, 1)/GPI
g(nr2, 1: nz2) = integr(1: nz2, 2)/GPI
g(1: nr2, nz2) = integr(1: nr2, 3)/GPI
g(1, 1: nz2)   = integr(1: nz2, 4)/GPI

return
end subroutine boundary_ef

!------------------------------------------------------------
subroutine bgint_ef(bgintsol, g, jcounty)
! calculates  integral_over_boundary of -Green * dg/dn * dl for point r0, z0

use ef_circuit, only: i_dim2, nr1, nz1, nr2, nz2, &
    r, dr, dz, green_bnd_f
use debugger_ef, only: markloc_ef, debug

implicit none

double precision, intent(in) :: g(i_dim2, i_dim2)
double precision, intent(out) :: bgintsol
integer, intent(inout) :: jcounty

integer :: j
double precision :: dgdn(i_dim2), greenf

call markloc_ef('bgint_ef', debug_lev=debug)

bgintsol=0.
! lower side
do j=1, nr1
    jcounty=jcounty+1
    greenf=green_bnd_f(jcounty)
    dgdn(j)=-(g(j, 2)-g(j, 1)+g(j+1, 2)-g(j+1, 1))/2./dz*greenf*dr/(r(j)+dr/2.)
enddo
bgintsol=bgintsol-sum(dgdn(1: nr1))
!right side    
do j=1, nz1
    jcounty=jcounty+1
    greenf=green_bnd_f(jcounty)
    dgdn(j)=(g(nr2, j)-g(nr1, j)+g(nr2, j+1)-g(nr1, j+1))/2./dr*greenf*dz/(r(nr2)+r(nr1))*2.
enddo
bgintsol=bgintsol-sum(dgdn(1: nz1))
! upper side
do j=1, nr1
    jcounty=jcounty+1
    greenf=green_bnd_f(jcounty)
    dgdn(j)=(g(j, nz2)-g(j, nz1)+g(j+1, nz2)-g(j+1, nz1))/2./dz*greenf*dr/(r(j)+dr/2.)
enddo
bgintsol=bgintsol-sum(dgdn(1: nr1))
!left side    
do j=1, nz1
    jcounty=jcounty+1
    greenf=green_bnd_f(jcounty)
    dgdn(j)=-(g(2, j)-g(1, j)+g(2, j+1)-g(1, j+1))/2./dr*greenf*dz/(r(1)+r(2))*2.
enddo
bgintsol=bgintsol-sum(dgdn(1: nz1))

return
end subroutine bgint_ef
