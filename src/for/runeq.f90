!---------------------------------------------------------------------
subroutine RUNEQ(G_new, H_new, G_old, H_old, y_old, W_in, V_in, unit_coeff, G11, A_in, B_in, R_in, Src_y_in, Src_in, &
    rbdot, bbdot, Ngridb, Ngrid, dx, dt, x_in, imethod, bctype, bc_values, y_out, Q_out, adcmp_term, mphit)
!---------------------------------------------------------------------
! WARNING: at the moment Qb is explicit, no option for QNNB, QETB, QITB is given at the moment!
!
! also, schemes are not appropriate for A = 0, also the Grigori's scheme with additional D and V for
! stiff transport is not yet implemented
!
!  Inputs: G, V, Src are on main grid
!          A, B, R, G11 are on shifted grid
!      dx, dt, x (main grid, 1:Ngrid), should be RHO
! R_in is non-zero only in the momentum equation (residual stress)
!
!  bctype=2,3 if YB isn't set and QB is set
!   then C(2)=QB, otherwise, if bctype=1,4, use yb = y(Ngrid)
!   In the case bctype=2,3, solves up to Ngridb
!   In the case bctype=1,4, solves up to Ngridb-1 and uses Ngridb as b.c.
!   If bctype = 3, then do mixed b.c., where C(5)*y(b)+C(6)*y(b-1) = C(7)
!
!  Outputs: y(1:Ngrid), Q(1:Ngrid).
!
!  *N are at t+dt
!  *O are at t
!
! implicit adiabatic compression 
!
! 1/G d/dt(G*H*y)+1/V d/dx(unit_coeff*(Q-bdot*(V/W*N-H/M)*x*y-rdot*x*H/M*y)) = 
!        (S-bdot*M*N*x*d/dx(V/W)-(rdot-bdot)/V*G*H*d/dx(x*V/G))*y+P
!
!  rbdot = phibdot /(2*phib), bdot = Bdot/(2*B)
!
! S = S - bdot*M*N*x*d/dx(V/W) - (rdot-bdot)/V*G*H*d/dx (x*V/G)
! B = B + 1/G11*(bdot*V/W*N+(rdot-bdot)*H/M)*x
!
! y is the quantity and Q = -G11*(A*dy/dx + B*y) + G11*R
!
! imethod can be:
!      21 - implicit, centered differencies
!      22 - implicit, power law scheme
!      23 - exponential scheme
!      31 - Crank-Nicholson, centered differencies
!      32 - Crank-Nicholson, power law scheme
!      33 - Crank-Nicholson, exponential scheme
!
!  bctype = 1 -> f_bound
!  bctype = 2 -> Qbound
!  bctype = 3 -> mixed
!  bctype = 4 -> Q/y BC

use numerical_tools, only: extrap, deriv, gradient

implicit none

integer, intent(in) :: Ngrid, imethod, Ngridb, bctype
double precision, intent(in) :: dx, dt, rbdot, bbdot, unit_coeff, bc_values(5)
double precision, intent(in), dimension(Ngrid) :: G_new, H_new, G_old, H_old, &
   y_old, V_in, G11, A_in, B_in, R_in, Src_y_in, Src_in, mphit, x_in, W_in
double precision, intent(out), dimension(Ngrid) :: y_out, Q_out, adcmp_term

integer :: j
double precision :: theta, bc_val(3)
double precision, dimension(Ngrid) :: x_b, Gmid, Hmid, NN, NO, &
    Rsource, Pdot_1, Pdot_2, Src_new, ytmp
double precision, dimension(Ngridb) :: Gmix, AA, BB, CC, RR, g_v, vta, gsdt, gsydt, xi, fxi, gxi

Gmid = 0.5*(G_old + G_new)
Hmid = 0.5*(H_old + H_new)
x_b = x_in + 0.5*dx

! Compute source: Rsource = - 1/V d/dx (unit_coeff*G11*R), Rsource is on main grid
! Note that R_in (residual term) is always 0 except for the momentum equation
call DERIV(x_b, x_in(Ngrid), 2, unit_coeff*G11*R_in, Rsource, 1, Ngrid, 0)

if (sum(mphit) == 0.) then
    ytmp = y_old
else
    ytmp = mphit
endif

call GRADIENT(x_in, H_new*W_in*x_in*ytmp, Pdot_1, Ngrid) ! H_new*W_in = N_in*unit_coeff, as in pyparse YWNB = YWHN*YWWB/unit_coeff
call GRADIENT(x_in, Gmid*Hmid*ytmp, Pdot_2, Ngrid)
Pdot_1 = Pdot_1/W_in
Pdot_2 = Pdot_2*x_in/Gmid

adcmp_term = bbdot*Pdot_1 + (rbdot - bbdot)*Pdot_2  ! output, used only for FP ans Src_new
Src_new = Src_in - Rsource/V_in + adcmp_term

! Default is imethod=22 (INUME1-4 in const.f90)
SELECT CASE(imethod)
CASE(21: 23)
    theta = 1.
CASE(31: 33)
    theta = 0.5
END SELECT

! Define cd, or power law
do j=1, Ngridb
    xi(j) = dx*B_in(j)/max(A_in(j), 1.E-16)
enddo
call calc_fxi(imethod, Ngridb, xi, fxi)
gxi = fxi - xi

! 1/G d/dt (N*y) + 1/V d/dx (unit_coeff*G11*(-A/dx*(fxi, gxi, ntilde))) = S*y + P

NN = G_new*H_new
NO = G_old*H_old
Gmix = G_new(1: Ngridb)*theta + G_old(1: Ngridb)*(1. - theta)
g_v = dt/dx**2 * Gmix/V_in(1: Ngridb)
vta = unit_coeff*G11(1: Ngridb)*A_in(1: Ngridb)
gsdt  = dt*Gmix*Src_new(1: Ngridb)
gsydt = dt*Gmix*Src_y_in(1: Ngridb)

AA = -theta*g_v*vta*fxi
CC(1) = 0.
CC(2: Ngridb) = -g_v(2: Ngridb)*vta(1: Ngridb-1)*gxi(1: Ngridb-1)

BB(1) = NN(1) + theta*(-gsydt(1) + g_v(1)*vta(1)*gxi(1))
RR(1) = NO(1)*y_old(1) + gsdt(1) + (1. - theta) * ( gsydt(1)*y_old(1) + g_v(1)*vta(1)*(fxi(1)*y_old(2) - gxi(1)*y_old(1)) )
do j=2, Ngridb-1
    BB(j) = NN(j) + theta * ( -gsydt(j) + g_v(j)*(vta(j)*gxi(j) + vta(j-1)*fxi(j-1)) )
    RR(j) = NO(j)*y_old(j) + gsdt(j) + (1. - theta) * (  gsydt(j)*y_old(j) + g_v(j) * &
        ( vta(j) * (fxi(j)*y_old(j+1) - gxi(j)*y_old(j)) - vta(j-1) * (fxi(j-1)*y_old(j) - gxi(j-1)*y_old(j-1)) )  )
enddo

SELECT CASE(bctype)
CASE(1)
    bc_val(1) = y_old(Ngridb)
CASE(3)
    bc_val = bc_values(3: 5)
CASE(2, 4)
    bc_val(1) = NN(Ngridb) + theta * ( -gsydt(Ngridb) + g_v(Ngridb)*vta(Ngridb-1)*fxi(Ngridb-1) )
    if (bctype == 4) bc_val(1) = bc_val(1) + theta*g_v(Ngridb)*bc_values(2)*unit_coeff*dx
    bc_val(2) = CC(Ngridb)
    bc_val(3) = NO(Ngridb)*y_old(Ngridb) + gsydt(Ngridb) + (1. - theta) * &
        (  gsydt(Ngridb)*y_old(Ngridb) - g_v(Ngridb) * &
        (vta(Ngridb-1) * (fxi(Ngridb-1)*y_old(Ngridb) - gxi(Ngridb-1)*y_old(Ngridb-1)) )  )
    if (bctype == 2) bc_val(3) = bc_val(3) + g_v(Ngridb)*bc_values(1)*unit_coeff*dx
END SELECT

! Tridiagonal solver
call TRIDIAG(Ngridb, bctype, bc_val, AA, BB, CC, RR, y_out(1: Ngridb))

! Compute flux from solution (flux is on shifted grid)
do j=1, Ngridb-1
    Q_out(j) = G11(j)*(-A_in(j)/dx * (fxi(j)*y_out(j+1) - gxi(j)*y_out(j)) + R_in(j))
enddo

SELECT CASE(bctype)
CASE(1, 3)
    Q_out(Ngridb) = EXTRAP(x_in(1: Ngridb-1), Q_out(1: Ngridb-1), x_in(Ngridb), Ngridb-1, 2, .false.)
CASE(2)
    Q_out(Ngridb) = bc_values(1)
CASE(4)
    Q_out(Ngridb) = bc_values(2)*y_out(Ngridb)
END SELECT

return
end subroutine RUNEQ

!---------------------------------------------------------------------
subroutine calc_fxi(imethod, NgridB, xi, fxi)

implicit none

integer, intent(in) :: imethod, Ngridb
double precision, intent(in), dimension(Ngridb) :: xi
double precision, intent(out), dimension(Ngridb) :: fxi

integer :: j

fxi = 0.

SELECT CASE(imethod)
CASE(21, 31)
    fxi = 1. + 0.5*xi
CASE(22, 32)
    do j=1, Ngridb
        if (xi(j) >= -10. .and. xi(j) < 0) then
            fxi(j) = (1. + 0.1*xi(j))**5
        endif
        if (xi(j) >= 0. .and. xi(j) <= 10) then
            fxi(j) = (1. - 0.1*xi(j))**5 + xi(j)
        endif
        if (xi(j) > 10.) then
            fxi(j) = xi(j)
        endif
    enddo
CASE(23, 33)
    do j=1, Ngridb
        if (xi(j) == 0) then
            fxi(j) = 1.
        else
            fxi(j) = xi(j)/(1. - exp(-xi(j)))
        endif
    enddo
 END SELECT

return
end subroutine calc_fxi

!---------------------------------------------------------------------
subroutine TRIDIAG(Ngridb, bctype, bc_val, AA_in, BB_in, CC_in, RR_in, y_out)
!---------------------------------------------------------------------
! Tridiagonal solver. It solves the system:
!
!   Aj fj+1  + Bj fj  + Cj fj-1 = Rj
!   for j=1, Ngridb

implicit none

integer, intent(in) :: Ngridb, bctype
double precision, intent(in), dimension(3) :: bc_val
double precision, intent(in), dimension(Ngridb) :: AA_in, BB_in, CC_in, RR_in
double precision, intent(out), dimension(Ngridb) :: y_out

integer :: j
double precision, dimension(Ngridb-1) :: alpha, beta

alpha(1) = -BB_in(1)/AA_in(1)
beta(1)  =  RR_in(1)/AA_in(1)
do j=2, Ngridb-1
    alpha(j) = -CC_in(j)/(AA_in(j)*alpha(j-1)) - BB_in(j)/AA_in(j)
    beta(j) = RR_in(j)/AA_in(j) + CC_in(j)/AA_in(j)*beta(j-1)/alpha(j-1)
enddo

! Note that boundary value is assumed to be on the main last grid point, so 1-dx/2.
! This has to be corrected later on... 

if (bctype == 1) then
    y_out(Ngridb) = bc_val(1)
else
    y_out(Ngridb) = (bc_val(3) + bc_val(2)*beta(Ngridb-1)/alpha(Ngridb-1))/(bc_val(1) + bc_val(2)/alpha(Ngridb-1))
endif
do j=Ngridb-1, 1, -1
    y_out(j) = (y_out(j+1) - beta(j))/alpha(j)
enddo

return
end subroutine TRIDIAG

!---------------------------------------------------------------------
double precision function GETPEI(j)
!------------------
! Implicit scheme !
!------------------

use status_inc, only: PEIQI, NE, NI, TE, TI, ZMAIN, AMAIN, NMAIN, & 
    NIZ1, NIZ2, NIZ3, ZIM1, ZIM2, ZIM3, PET, PIT
use const_inc, only: IPROT, AIM1, AIM2, AIM3

implicit none

integer, intent(in) :: j
double precision :: COULG, SUZPEI, t1, t2

GETPEI = PEIQI(j)
t1 = PEIQI(j)
t2 = TE(j) - TI(j)

if (t1 == 0. .and. t2 /= 0.) then
    COULG = 15.9 - 0.5*LOG(NE(J)) + log(TE(J))
    if (nint(abs(IPROT)) == 2 .or. nint(ABS(IPROT)) == 4) then
        SUZPEI = NMAIN(J)/AMAIN(J)*ZMAIN(J)**2 + &
            NIZ1(J)/AIM1*ZIM1(J)**2 + &
            NIZ2(J)/AIM2*ZIM2(J)**2 + &
            NIZ3(J)/AIM3*ZIM3(J)**2
        GETPEI = 0.00246*COULG*NE(J)*SUZPEI/(TE(J)*SQRT(TE(J)))
    else
        GETPEI = 0.00246*COULG*NE(J)*NI(J)*ZMAIN(J)**2/(AMAIN(J)*TE(J)*SQRT(TE(J)))
    endif
endif

! This will reflect equipartition in PETOT and PITOT after equations are solved, for plotting and post-processing.
PET(J) = PET(J) - GETPEI*t2/TE(J)
PIT(J) = PIT(J) + GETPEI*t2/TI(J)

return
end function GETPEI

!---------------------------------------------------------------------
subroutine RUNEQ_TETI(G_new, H1_new, H2_new, G_old, H1_old, H2_old, & 
    Y1_old, Y2_old, W1, W2, V, unit_coeff, G11, A1_in, A2_in, B1_in, B2_in, &
    S1_in, S2_in, P1, P2, rbdot, bbdot, Ngridb, Ngrid, dx, dt, &
    x_in, imethod, y1, y2, Q1_out, Q2_out)

use numerical_tools, only: deriv, extrap, gradient

implicit none

integer, intent(in) :: Ngrid, imethod, Ngridb
double precision, intent(in) :: rbdot, bbdot, dx, dt, unit_coeff
double precision, intent(in), dimension(Ngrid) :: G_new, H1_new, H2_new, G_old, &
    H1_old, H2_old, Y1_old, Y2_old, W1, W2, V, G11, &
    A1_in, A2_in, B1_in, B2_in, S1_in, S2_in, P1, P2, x_in
double precision, intent(out)  , dimension(Ngrid) :: Q1_out, Q2_out
double precision, intent(inout), dimension(Ngrid) :: y1, y2

integer :: j
double precision :: theta, f_bound1, f_bound2
double precision, dimension(Ngrid) :: x_b, &
    Vtilde, N1_new, N1_old, N2_new, N2_old
double precision, dimension(Ngridb) :: A1, A2, T12, &
    S1_new, S2_new, xi1, fxi1, gxi1, xi2, fxi2, gxi2, &
    Gmix, g_v, AA1, BB1, CC1, RR1, TT1, AA2, BB2, CC2, RR2, TT2
double precision, external :: GETPEI

x_b = x_in + 0.5*dx
do j=1, Ngridb
    A1(j) = max(A1_in(j),  1.E-16)
    A2(j) = max(A2_in(j),  1.E-16)
enddo

SELECT CASE(imethod)
CASE(21: 23)
    theta = 1.
CASE(31: 33)
    theta = 0.5
END SELECT

Vtilde = unit_coeff*G11

N1_new = G_new*H1_new
N2_new = G_new*H2_new
N1_old = G_old*H1_old
N2_old = G_old*H2_old

do j=1, Ngridb
    T12(j) = 625.*GETPEI(j)
enddo

S1_new = S1_in(1: Ngridb) - T12
S2_new = S2_in(1: Ngridb) - T12

f_bound1 = Y1_old(Ngridb)
f_bound2 = Y2_old(Ngridb)

! Define cd, or power law
xi1 = dx*B1_in(1: Ngridb)/A1
xi2 = dx*B2_in(1: Ngridb)/A2
call calc_fxi(imethod, Ngridb, xi1, fxi1)
call calc_fxi(imethod, Ngridb, xi2, fxi2)
gxi1 = fxi1 - xi1
gxi2 = fxi2 - xi2

Gmix = G_new(1: Ngridb)*theta + G_old(1: Ngridb) * (1. - theta)
g_v = dt/dx**2 * Gmix/V(1: Ngridb)
AA1 = -g_v*Vtilde(1: Ngridb)*A1(1: Ngridb)*fxi1*theta
AA2 = -g_v*Vtilde(1: Ngridb)*A2(1: Ngridb)*fxi2*theta
TT1 = -theta*dt*Gmix*T12(1: Ngridb)
CC1(1) = 0.
CC2(1) = 0.
CC1(2: Ngridb) = -g_v(2: Ngridb)*Vtilde(1: Ngridb-1)*A1(1: Ngridb-1)*gxi1(1: Ngridb-1)*theta
CC2(2: Ngridb) = -g_v(2: Ngridb)*Vtilde(1: Ngridb-1)*A2(1: Ngridb-1)*gxi2(1: Ngridb-1)*theta

BB1(1) = N1_new(1) + theta * ( -dt*Gmix(1)*S1_new(1) + g_v(1)*Vtilde(1)*A1(1)*gxi1(1) )
BB2(1) = N2_new(1) + theta * ( -dt*Gmix(1)*S2_new(1) + g_v(1)*Vtilde(1)*A2(1)*gxi2(1) )

RR1(1) = N1_old(1)*Y1_old(1) + Gmix(1)*P1(1)*dt + (1 - theta) * &
    ( dt*Gmix(1)*S1_new(1)*Y1_old(1) + g_v(1)*Vtilde(1)*A1(1) * (fxi1(1)*Y1_old(2) - gxi1(1)*Y1_old(1)) )
RR2(1) = N2_old(1)*Y2_old(1) + Gmix(1)*P2(1)*dt + (1 - theta) * &
    ( dt*Gmix(1)*S2_new(1)*Y2_old(1) + g_v(1)*Vtilde(1)*A2(1) * (fxi2(1)*Y2_old(2) - gxi2(1)*Y2_old(1)) )

do j=2, Ngridb-1
    BB1(j) = N1_new(j) + theta * (  -dt*Gmix(j)*S1_new(j) + g_v(j) * &
        ( Vtilde(j)*A1(j)*gxi1(j) + Vtilde(j-1)*A1(j-1)*fxi1(j-1) )  )
    BB2(j) = N2_new(j) + theta * ( -dt*Gmix(j)*S2_new(j) + g_v(j) * &
        (Vtilde(j)*A2(j)*gxi2(j) + Vtilde(j-1)*A2(j-1)*fxi2(j-1)) )
    RR1(j) = N1_old(j)*Y1_old(j) + Gmix(j)*P1(j)*dt + (1 - theta) * &
        (  dt*Gmix(j)*S1_new(j)*Y1_old(j) + g_v(j) * &
        (Vtilde(j)*A1(j) * (fxi1(j)*Y1_old(j+1) - gxi1(j)*Y1_old(j)) - &
        Vtilde(j-1)*A1(j-1) * (fxi1(j-1)*Y1_old(j) - gxi1(j-1)*Y1_old(j-1)))  )
    RR2(j) = N2_old(j)*Y2_old(j) + Gmix(j)*P2(j)*dt + (1 - theta) * &
        (  dt*Gmix(j)*S2_new(j)*Y2_old(j) + g_v(j) * &
        (Vtilde(j)*A2(j) * (fxi2(j)*Y2_old(j+1) - gxi2(j)  *Y2_old(j)) - &
         Vtilde(j-1)*A2(j-1) * (fxi2(j-1)*Y2_old(j) - gxi2(j-1)*Y2_old(j-1)))  )
enddo

! Main call to tridiagonal solver
call TRIDIAG_TETI(AA1, BB1, CC1, RR1, y1, TT1, Ngridb, &
    f_bound1, AA2, BB2, CC2, RR2, TT1, y2, f_bound2)

! Compute flux from solution, flux is on shifted grid
do j=1, Ngridb-1
    Q1_out(j) = G11(j)*(-A1(j)/dx*(fxi1(j)*y1(j+1) - gxi1(j)*y1(j)))
    Q2_out(j) = G11(j)*(-A2(j)/dx*(fxi2(j)*y2(j+1) - gxi2(j)*y2(j)))
enddo
Q1_out(Ngridb) = EXTRAP(x_in(1: Ngridb-1), Q1_out(1: Ngridb-1), x_in(Ngridb), Ngridb-1, 1, .false.)
Q2_out(Ngridb) = EXTRAP(x_in(1: Ngridb-1), Q2_out(1: Ngridb-1), x_in(Ngridb), Ngridb-1, 1, .false.)

return
end subroutine RUNEQ_TETI

!---------------------------------------------------------------------
subroutine TRIDIAG_TETI(A1, B1, C1, R1, f1, T1, Ngridb, &
    f_bound1, A2, B2, C2, R2, T2, f2, f_bound2)
!---------------------------------------------------------------------
! Inputs: A(1:Ngridb), B(1:Ngridb), C(1:Ngridb), R(1:Ngridb)
! Outputs: f(1:Ngridb)
! 
! Provides solution of the system:
!
!   A1j f1j+1  + B1j f1j  + C1j f1j-1 + T1j f2j = R1j
!   A2j f2j+1  + B2j f2j  + C2j f2j-1 + T2j f1j = R2j
!
!   where j = 1..Ngridb
!
!   given f_NA1
!---------------------------------------------------------------------

implicit none

integer, intent(in) :: Ngridb
double precision, intent(in) :: f_bound1, f_bound2
double precision, intent(in), dimension(Ngridb) :: &
    A1, B1, C1, R1, T1, &
    A2, B2, C2, R2, T2
double precision, intent(out), dimension(Ngridb) :: f1, f2

integer :: j
double precision :: detjm1, k1, k2, k3, k4, k31, k32
double precision, dimension(Ngridb-1) :: alpha, beta, gamma, &
    delta, epsilon, theta

alpha(1)   = -B1(1)/A1(1)
beta(1)    =  R1(1)/A1(1)
epsilon(1) = -T1(1)/A1(1)      
gamma(1)   = -B2(1)/A2(1)
delta(1)   =  R2(1)/A2(1)
theta(1)   = -T2(1)/A2(1)

do j=2, Ngridb-1
    detjm1 = (alpha(j-1)*gamma(j-1) - theta(j-1)*epsilon(j-1))
    k1 = gamma(j-1)/detjm1
    k3 = alpha(j-1)/detjm1
    k2 = -epsilon(j-1)/detjm1
    k4 = -theta(j-1)/detjm1
    k31 = gamma(j-1)/detjm1*(-beta(j-1)) + epsilon(j-1)/detjm1*(delta(j-1))
    k32 = alpha(j-1)/detjm1*(-delta(j-1)) + theta(j-1)/detjm1*(beta(j-1))

    alpha(j)   = -C1(j)*k1/A1(j) - B1(j)/A1(j)
    beta(j)    =  R1(j)/A1(j)    - C1(j)/A1(j)*k31
    epsilon(j) = -T1(j)/A1(j)    - C1(j)*k2/A1(j)
    gamma(j)   = -C2(j)*k3/A2(j) - B2(j)/A2(j)
    delta(j)   =  R2(j)/A2(j) - C2(j)/A2(j)*k32
    theta(j)   = -T2(j)/A2(j) - C2(j)*k4/A2(j)
enddo

! Note that boundary value is assumed to be on the main last grid point, so 1-dx/2.
! This has to be corrected later on

detjm1 = (alpha(Ngridb-1)*gamma(Ngridb-1) - theta(Ngridb-1)*epsilon(Ngridb-1))
f1(Ngridb-1) =  gamma(Ngridb-1)/detjm1*(f_bound1 - beta(Ngridb-1)) - &
              epsilon(Ngridb-1)/detjm1*(f_bound2 - delta(Ngridb-1))
f2(Ngridb-1) = -theta(Ngridb-1)/detjm1*(f_bound1 - beta(Ngridb-1)) + &
                alpha(Ngridb-1)/detjm1*(f_bound2 - delta(Ngridb-1))
do j=Ngridb-2, 1, -1
     detjm1 = (alpha(j)*gamma(j) - theta(j)*epsilon(j))
     f1(j) =  gamma(j)/detjm1*(f1(j+1) - beta(j)) - &
            epsilon(j)/detjm1*(f2(j+1) - delta(j))
     f2(j) = -theta(j)/detjm1*(f1(j+1) - beta(j)) + &
              alpha(j)/detjm1*(f2(j+1) - delta(j))
enddo
f1(Ngridb) = f_bound1
f2(Ngridb) = f_bound2

return
end subroutine TRIDIAG_TETI
