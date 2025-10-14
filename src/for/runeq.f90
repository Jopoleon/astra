!---------------------------------------------------------------------
subroutine RUNEQ(GN, HN, GO, HO, YO, N, W, V, M, G11, A_in, B_in, R_in, Src_y_in, Src_in, &
    rbdot, bbdot, Ngridb, Ngrid, dx, dt, &
    x, imethod, bctype, bc_values, y_out, Q_out, adcmp_term, mphit)
!---------------------------------------------------------------------
!
! WARNING: at the moment Qb is explicit, no option for QNNB, QETB, QITB is given at the moment!
!
! also, schemes are not appropriate for A = 0, also the Grigori's scheme with additional D and V for
! stiff transport is not yet implemented
!
!  Inputs: G, V, A, B, R, S, G11, M, P in terms of Ngrid
!  Note that G, V, M, S, P are on main grid
!   while A, B, R, G11 are on shifted grid
!      dx, dt, x (main grid, 1:Ngrid), should be RHO, imethod
!      C: boundary conditions on y or on Q  
!  bctype=2,3 if (yb isn't set) .and. (QB is set)
!   then C(2)=QB, otherwise, if bctype=1,4, use yb = y(Ngrid)
!   In the case bctype=2,3, solves up to Ngrid
!   In the case bctype=1,4, solves up to Ngrid-1 and uses Ngrid as b.c.
!   If bctype = 3, then do mixed b.c., where C(5)*y(b)+C(6)*y(b-1) = C(7)
!
!  Outputs: y(1:Ngrid), Q(1:Ngrid).
!
!  *N are at t+dt
!  *O are at t
!
!  Equation is:
!
! explicit adiabatic compression
!
! 1/G d/dt(G*H*y)+1/V d/dx(V*M*Q) = S*y+P + bdot/W*d/dx (V*M*N*x*y)+(rdot-bdot)*x/G*d/dx (G*H*y) 
!
! implicit adiabatic compression 
!
! 1/G d/dt(G*H*y)+1/V d/dx(V*M*(Q-bdot*(V/W*N-H/M)*x*y-rdot*x*H/M*y)) = 
!        (S-bdot*M*N*x*d/dx(V/W)-(rdot-bdot)/V*G*H*d/dx(x*V/G))*y+P
!
!  rbdot = phibdot /(2*phib), bdot = Bdot/(2*B)
!
! So S = S - bdot*M*N*x*d/dx(V/W) - (rdot-bdot)/V*G*H*d/dx (x*V/G)      for implicit
!
! So B = B + 1/G11*(bdot*V/W*N+(rdot-bdot)*H/M)*x                       for implicit
!
! and P = P + bdot/W*d/dx (V*M*N*x*y)+(rdot-bdot)*x/G*d/dx (G*H*y)      for explicit
!
!  y is the quantity and Q = -G11*(A*dy/dx + B*y) + G11*R
!
!  imethod can be:
!      11 - explicit, centered differencies
!      12 - explicit, power law scheme
!      13 - explicit, exponential scheme
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

use numerical_tools, only: extrap, deriv, grid2grid

implicit none

integer, intent(in) :: Ngrid, imethod, Ngridb, bctype
double precision, intent(in) :: rbdot, bbdot, bc_values(5)
double precision, intent(in), dimension(Ngrid) :: GN, HN, GO, HO, &
   YO, V, G11, A_in, B_in, R_in, Src_y_in, Src_in, mphit, x, M, N, W
double precision, intent(out)  , dimension(Ngrid) :: y_out, Q_out, adcmp_term

integer :: eximp, j
double precision :: dx, dt, theta, f_bound, Qbound, Mbound(3)
double precision, dimension(Ngrid) :: x_b, G, H, A, NN, NO, dum1, &
    Rsource, Rsource2, Pdot_1, Pdot_2, Vtilde, Src_new, ydummy, ytmp, xi, fxi, gxi
double precision, dimension(Ngridb) :: dum1b

do j=1, Ngrid
    G(j) = 0.5*(GO(j) + GN(j))
    H(j) = 0.5*(HO(j) + HN(j))
    x_b(j) = x(j) + 0.5*dx
    A(j) = max(A_in(j), 1.E-16)
enddo

! Define explicit or implicit; default is imethod=22 (INUME1-4 in const.f90)
SELECT CASE(imethod)
CASE(11: 13)
    theta = 0.
    eximp = 1
CASE(21: 23)
    theta = 1.
    eximp = 2
CASE(31: 33)
    theta = 0.5
    eximp = 2
END SELECT
   
! Map V*M on tildes

dum1 = V*M
call GRID2GRID(1, x, dum1, Vtilde, Ngrid, 1)

! Compute source: Rsource = - 1/V d/dx (V*M*G11*R), Rsource is on main grid
Rsource = Vtilde*G11*R_in
call DERIV(x_b, x, 2, Rsource, Rsource2, 1, Ngrid, 0)
Src_new = Src_in - Rsource2/V

Pdot_1 = 0.
Pdot_2 = 0.
if (sum(mphit) == 0.) then
    ytmp = YO
else
    ytmp = mphit
endif
dum1 = V*M*N*x*ytmp
call DERIV(x, x_b, 1, dum1, ydummy, 1, Ngrid, 1)
call GRID2GRID(2, x_b, ydummy, Pdot_1, Ngrid, 0)
Pdot_1 = Pdot_1/W
dum1 = G*H*ytmp
call DERIV(x, x_b, 1, dum1, ydummy, 1, Ngrid, 1)
call GRID2GRID(2, x_b, ydummy, Pdot_2, Ngrid, 0)
Pdot_2 = Pdot_2*x/G

NN = GN*HN
NO = GO*HO
Vtilde = Vtilde*G11
! Src = Src + bdot/W*d/dx (V*M*N*x*y)+(rdot-bdot)*x/G*d/dx (G*H*y)    for explicit
adcmp_term = bbdot*Pdot_1 + (rbdot - bbdot)*Pdot_2  ! used only for FP
Src_new = Src_new + adcmp_term

! Check boundary condition, note that Qbound = Qbound/G11(b) since G11 is absorbed in Vtilde
Mbound  = 0.
f_bound = 0.
Qbound  = 0.
SELECT CASE(bctype)
CASE(1)
    f_bound = y_out(Ngridb)
CASE(2)
    Qbound = bc_values(1)/G11(Ngridb)
CASE(3)
    Mbound(1: 3) = bc_values(3: 5)
CASE(4)
    Qbound = bc_values(2)/G11(Ngridb)
END SELECT
 
! Define cd, or power law
xi = dx*B_in/A
SELECT CASE(imethod)
CASE(11, 21, 31)
    fxi = 1. + 0.5*xi
CASE(12, 22, 32)
    do j=1, Ngrid
        if (xi(j) < -10.) then
            fxi(j) = 0.
        endif
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
CASE(13, 23, 33)
    do j=1, Ngrid
        if (xi(j) == 0) then
            fxi(j) = 1.
        else
            fxi(j) = xi(j)/(1. - exp(-xi(j)))
        endif
    enddo
END SELECT
gxi = fxi - xi

! Call solver with this equation, note that A and B are on shifted grids
!
! Vtilde = (V*M)_shifted*G11
! 1/G d/dt (N*y) + 1/V d/dx (Vtilde*(-A/dx*(fxi, gxi, ntilde))) = S*y + P 

dum1b = GN(1: Ngridb)*theta + GO(1: Ngridb)*(1. - theta)
call SOLVER(dx, dt, dum1b, NN(1: Ngridb), NO(1: Ngridb), &
    V(1: Ngridb), Vtilde(1: Ngridb), A(1: Ngridb), &
    fxi(1: Ngridb), gxi(1: Ngridb), &
    Src_y_in(1: Ngridb), Src_new(1: Ngridb), &
    Ngridb, theta, eximp, YO(1: Ngridb), &
    f_bound, Qbound, bctype, Mbound, y_out(1: Ngridb))

! Compute flux from solution (flux is on shifted grid)
do j=1, Ngridb-1
    Q_out(j) = G11(j)*(-A(j)/dx * (fxi(j)*y_out(j+1) - gxi(j)*y_out(j)) + R_in(j))
enddo

SELECT CASE(bctype)
CASE(1, 3)
    call EXTRAP(x(1: Ngridb-1), Q_out(1: Ngridb-1), x(Ngridb), Ngridb-1, Q_out(Ngridb), 2, Ngridb-1)
CASE(2)
    Q_out(Ngridb) = Qbound*G11(Ngridb)
CASE(4)
    Q_out(Ngridb) = Qbound*G11(Ngridb)*y_out(Ngridb)
END SELECT

return
end subroutine RUNEQ

!---------------------------------------------------------------------
subroutine SOLVER(dx, dt, G, NN, NO, V, Vtilde, A, fxi, gxi, Src_y_in, Src_in, & 
    Ngridb, theta, eximp, YO, f_bound, Qbound, bctype, Mbound, y_out)
!---------------------------------------------------------------------
! Build up matrices to be passed ot TRIDIAGS

implicit none

integer, intent(in) :: Ngridb, bctype, eximp
double precision, intent(in) :: dx, dt, theta, f_bound, Qbound, Mbound(3)

double precision, intent(in), dimension(Ngridb) :: YO, V, A, &
    Src_y_in, Src_in, G, Vtilde, NN, NO, fxi, gxi
double precision, intent(out), dimension(Ngridb) :: y_out

integer :: j
double precision :: dt_dx2
double precision, dimension(Ngridb) :: AA, BB, CC, RR

dt_dx2 = dt/dx**2

if (eximp == 1) then ! Explicit scheme
    RR(1) = NO(1)*YO(1) + G(1)*Src_in(1)*dt + dt*G(1)*Src_y_in(1)*YO(1) + &
        dt_dx2 * G(1)/V(1) * Vtilde(1) * ( A(1)*(fxi(1)*YO(2) - gxi(1)*YO(1)) )
    y_out(1) = RR(1)/NN(1)
    do j=2, Ngridb-1
        RR(j) = NO(j)*YO(j) + G(j)*Src_in(j)*dt + dt*G(j)*Src_y_in(j)*YO(j) + dt_dx2*G(j)/V(j) * ( &
            Vtilde(j)  *A(j)   * (fxi(j)   * YO(j+1) - gxi(j)  * YO(j)) - &
            Vtilde(j-1)*A(j-1) * (fxi(j-1) * YO(j)   - gxi(j-1)* YO(j-1)) )
        y_out(j) = RR(j)/NN(j)
    enddo

    SELECT CASE(bctype)
    CASE(1)
        y_out(Ngridb) = f_bound
    CASE(2)
        RR(Ngridb) = NO(Ngridb)*YO(Ngridb) + G(Ngridb)*Src_in(Ngridb)*dt + dt*G(Ngridb)*Src_y_in(Ngridb)*YO(Ngridb) - &
            dt_dx2*G(Ngridb)/V(Ngridb) * ( Vtilde(Ngridb)*Qbound*dx + &
            Vtilde(Ngridb-1)*A(Ngridb-1) * (fxi(Ngridb-1)*YO(Ngridb) - gxi(Ngridb-1)*YO(Ngridb-1)) )
        y_out(Ngridb) = RR(Ngridb)/NN(Ngridb)
    CASE(3)
        y_out(Ngridb) = Mbound(3)/Mbound(1)
    CASE(4)
        BB(Ngridb) = NN(Ngridb) + dt_dx2 * G(Ngridb)/V(Ngridb)*Vtilde(Ngridb)*Qbound*dx
        RR(Ngridb) = NO(Ngridb)*YO(Ngridb) + G(Ngridb)*Src_in(Ngridb)*dt + dt*G(Ngridb)*Src_y_in(Ngridb)*YO(Ngridb) - &
            dt_dx2*G(Ngridb)/V(Ngridb)*Vtilde(Ngridb-1)*A(Ngridb-1) * (fxi(Ngridb-1)*YO(Ngridb) - gxi(Ngridb-1)*YO(Ngridb-1))
        y_out(Ngridb) = RR(Ngridb)/BB(Ngridb)
    END SELECT

else ! Implicit, ASTRA default

    AA(Ngridb) = 0.
    AA(1: Ngridb-1) = -dt_dx2*theta *G(1: Ngridb-1)/V(1: Ngridb-1)*Vtilde(1: Ngridb-1)*A(1: Ngridb-1)*fxi(1: Ngridb-1)
    BB(1) = NN(1) + theta * (-dt*G(1)*Src_y_in(1) + dt_dx2*G(1)/V(1)*Vtilde(1)*A(1)*gxi(1) )
    CC(1) = 0.0
    RR(1) = NO(1)*YO(1) + G(1)*Src_in(1)*dt + (1. - theta) * &
        ( dt*G(1)*Src_y_in(1)*YO(1) + dt_dx2*G(1)/V(1)*Vtilde(1)*A(1)*(fxi(1)*YO(2) - gxi(1)*YO(1)) )

    do j=2, Ngridb-1
        BB(j) = NN(j) + theta * (  -dt*G(j)*Src_y_in(j) + dt_dx2*G(j)/V(j) * &
            (Vtilde(j)*A(j)*gxi(j) + Vtilde(j-1)*A(j-1)*fxi(j-1) )  )
        CC(j) = -dt_dx2*theta*G(j)/V(j)*Vtilde(j-1)*A(j-1)*gxi(j-1)
        RR(j) = NO(j)*YO(j) + G(j)*Src_in(j)*dt + (1. - theta) * &
            (  dt*G(j)*Src_y_in(j)*YO(j) + dt_dx2*G(j)/V(j) * &
            ( Vtilde(j)*A(j) * (fxi(j)*YO(j+1) - gxi(j)*YO(j)) - Vtilde(j-1)*A(j-1) * (fxi(j-1)*YO(j) - gxi(j-1)*YO(j-1)) )  )
    enddo

    SELECT CASE(bctype)
    CASE(1)
        BB(Ngridb) = 0.
        CC(Ngridb) = 0.
        RR(Ngridb) = 0.
    CASE(2)
        BB(Ngridb) = NN(Ngridb) + theta * ( -dt*G(Ngridb)*Src_y_in(Ngridb) + &
            dt_dx2*G(Ngridb)/V(Ngridb)*Vtilde(Ngridb-1)*A(Ngridb-1)*fxi(Ngridb-1) )
        CC(Ngridb) = -dt_dx2*theta*G(Ngridb)/V(Ngridb)*Vtilde(Ngridb-1)*A(Ngridb-1)*gxi(Ngridb-1)
        RR(Ngridb) = NO(Ngridb)*YO(Ngridb) + G(Ngridb)*Src_in(Ngridb)*dt + (1. - theta) * &
            (dt*G(Ngridb)*Src_y_in(Ngridb)*YO(Ngridb) - dt_dx2*G(Ngridb)/V(Ngridb) * &
            (Vtilde(Ngridb-1)*A(Ngridb-1) * (fxi(Ngridb-1)*YO(Ngridb) - gxi(Ngridb-1)*YO(Ngridb-1)) ) ) + &
            dt_dx2*G(Ngridb)/V(Ngridb)*Vtilde(Ngridb)*Qbound*dx
    CASE(3)
        BB(Ngridb) = Mbound(1)
        CC(Ngridb) = Mbound(2)   
        RR(Ngridb) = Mbound(3)
    CASE(4)
        BB(Ngridb) = NN(Ngridb) + theta * ( -dt*G(Ngridb)*Src_y_in(Ngridb) + &
            dt_dx2*G(Ngridb)/V(Ngridb)*Vtilde(Ngridb-1)*A(Ngridb-1)*fxi(Ngridb-1) ) + &
            theta*dt_dx2*G(Ngridb)/V(Ngridb)*Vtilde(Ngridb)*Qbound*dx
        CC(Ngridb) = -dt_dx2*theta*G(Ngridb)/V(Ngridb)*Vtilde(Ngridb-1)*A(Ngridb-1)*gxi(Ngridb-1)
        RR(Ngridb) = NO(Ngridb)*YO(Ngridb) + G(Ngridb)*Src_in(Ngridb)*dt + (1. - theta) * &
            (  dt*G(Ngridb)*Src_y_in(Ngridb)*YO(Ngridb) - dt_dx2*G(Ngridb)/V(Ngridb) * &
            ( Vtilde(Ngridb-1)*A(Ngridb-1) * (fxi(Ngridb-1)*YO(Ngridb) - gxi(Ngridb-1)*YO(Ngridb-1)) )  )
    END SELECT
    call TRIDIAG(Ngridb, bctype, f_bound, AA, BB, CC, RR, y_out)

endif

return
end subroutine SOLVER

!---------------------------------------------------------------------
subroutine TRIDIAG(Ngridb, bctype, f_bound, A_in, B_in, C_in, R_in, f_out)
!---------------------------------------------------------------------
! Tridiagonal solver. It solves the system:
!
!   Aj fj+1  + Bj fj  + Cj fj-1 = Rj
!   for j=1, Ngridb
!
!   bctype = 1  -> given f_NA1
!   bctype = 2  -> given Gamma_NA1
!   bctype = 3  -> mixed b.c.

implicit none

integer, intent(in) :: Ngridb, bctype
double precision, intent(in) :: f_bound
double precision, intent(in) , dimension(Ngridb) :: A_in, B_in, C_in, R_in
double precision, intent(out), dimension(Ngridb) :: f_out

integer :: j, k
double precision, dimension(Ngridb-1) :: alpha, beta

alpha(1) = -B_in(1)/A_in(1)
beta(1)  =  R_in(1)/A_in(1)
do j=2, Ngridb - 1
    alpha(j) = -C_in(j)/(A_in(j)*alpha(j-1)) - B_in(j)/A_in(j)
    beta(j) = R_in(j)/A_in(j) + C_in(j)/A_in(j)*beta(j-1)/alpha(j-1)
enddo

! Note that boundary value is assumed to be on the main last grid point, so 1-dx/2.
! This has to be corrected later on... 

if (bctype == 1) then
    f_out(Ngridb-1) = (f_bound - beta(Ngridb-1))/alpha(Ngridb-1)
    do k=1, Ngridb-2
        j = Ngridb - 2 - k + 1
        f_out(j) = (f_out(j+1) - beta(j))/alpha(j)
    enddo
    f_out(Ngridb) = f_bound
else
    f_out(Ngridb) = (R_in(Ngridb) + C_in(Ngridb)*beta(Ngridb-1)/alpha(Ngridb-1))/(B_in(Ngridb) + C_in(Ngridb)/alpha(Ngridb-1))
    do k=1, Ngridb-1
        j = Ngridb - 1 - k + 1
        f_out(j) = (f_out(j+1) - beta(j))/alpha(j)
    enddo
endif

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
subroutine RUNEQTIMP(GN, H1N, H2N, GO, H1O, H2O, & 
    Y1O, Y2O, N1, N2, W1, W2, V, M, G11, A1_in, A2_in, B1_in, B2_in, R1, R2, &
    S1_in, S2_in, P1, P2, T12, T21, rbdot, bbdot, Ngridb, Ngrid, dx, dt, &
    x, imethod, bctype, bcvalue, y1, y2, Q1_out, Q2_out)
!---------------------------------------------------------------------
! WARNING: at the moment Qb is explicit, no option for QNNB, QETB, QITB is given at the moment!
!
!  Inputs: G, V, A, B, R, S, G11, M, P in terms of Ngrid
!  Note that G, V, M, S, P are on main grid
!   while A, B, R, G11 are on shifted grid
!      dx, dt, x (main grid, 1:Ngrid), should be RHO, imethod
!      C: boundary conditions on y or on Q  
!  C(1) = HRO
!  bc_flag < 0 if (yb isn't set) .and. (QB is set)
!   then C(2)=QB, otherwise, if bc_flag>0, use yb = y(Ngrid)
!   In the case bc_flag < 0, solves up to Ngrid
!   In the case bc_flag > 0, solves up to Ngrid-1 and uses Ngrid as b.c.
!
!  Outputs: y(1:Ngrid), Q(1:Ngrid).  
!
!  GN, GO are G at t+dt and at t
!  HN, HO are H at t+dt and at t
!  YO  is y at t
!
!  Equations are:
!
!    1/G d/dt (G*H1*y1) + 1/V d/dx (V*M*Q1) = S1*y1+P1+T12*y2 + rbdot*xhat*1/G d/dx (Gtilde*Htilde1*y1)
!    1/G d/dt (G*H2*y2) + 1/V d/dx (V*M*Q2) = S2*y2+P2+T21*y1 + rbdot*xhat*1/G d/dx (Gtilde*Htilde2*y2)
!
!explicit adiabatic compression
!
! 1/G d/dt(G*H1*y1)+1/V d/dx(V*M*Q1) = S1*y1+P1+T12*y2 + bdot/W1*d/dx (V*M*N1*x*y1)+(rdot-bdot)*x/G*d/dx (G*H1*y1) 
!
!implicit adiabatic compression 
!
! 1/G d/dt(G*H*y)+1/V d/dx(V*M*(Q-bdot*(V/W*N-H/M)*x*y-rdot*x*H/M*y)) = 
!                                             (S-bdot*M*N*x*d/dx(V/W)-(rdot-bdot)/V*G*H*d/dx(x*V/G))*y+P
!
!  rbdot = phibdot /(2*phib), bdot = Bdot/(2*B)
!
! So S = S - bdot*M*N*x*d/dx(V/W) - (rdot-bdot)/V*G*H1*d/dx (x*V/G)  for implicit
!
! So B = B + 1/G11*(bdot*V/W*N+(rdot-bdot)*H1/M)*x               for implicit
!
! and P1 = P1 + bdot/W1*d/dx (V*M*N1*x*y1)+(rdot-bdot)*x/G*d/dx (G*H1*y1)         for explicit
!
! basically we do as in RUNEQ but y = [y1 y2] and so on (bigger matrix)
!
!  y is the quantity and Qj = -G11*(Aj*dyj/dx + Bj*yj) + G11*Rj
!
! Note that Bj -> Bj + rbdot/roc * sx * Htildej/(G11*M)
!
! and rbdot*xhat*1/G d/dx (Gtilde*Htildej*yj)  --> rbdot/roc*G*Hj/V*d/dx(sx*Vtilde/Gtilde)*yj
!
!
!  imethod can be: 21 - implicit, centered differencies
!          22 - implicit, power law scheme   
!                  31 - Crank-Nicholson, centered differencies
!          32 - Crank-Nicholson, power law scheme   
!
!  bctype = 1 -> f_bound
!  bctype = 2 -> Qbound
!---------------------------------------------------------------------

use numerical_tools, only: deriv, extrap, grid2grid

implicit none

integer, intent(in) :: Ngrid, imethod, Ngridb, bctype(2)
double precision, intent(in) :: rbdot, bbdot, dx, dt
double precision, intent(in), dimension(Ngrid) :: GN, H1N, H2N, GO, &
    H1O, H2O, Y1O, Y2O, N1, N2, W1, W2, V, M, G11, &
    A1_in, A2_in, B1_in, B2_in, &
    R1, R2, S1_in, S2_in, P1, P2, x, bcvalue(2)
double precision, intent(out)  , dimension(Ngrid) :: T12, T21, Q1_out, Q2_out
double precision, intent(inout), dimension(Ngrid) :: y1, y2

integer :: j, NgridS(2)
double precision :: theta, f_bound1, Qbound1, f_bound2, Qbound2
double precision, dimension(Ngrid) :: x_b, S1_new, P1_new, S2_new, P2_new, &
    Vtilde, Rsource1, Rsource2, Rsource3, Rsource4, N1N, N1O, N2N, N2O, ydummy, &
    xi1, fxi1, gxi1, xi2, fxi2, gxi2, A1, A2, Pdot_11, Pdot_12, Pdot_21, Pdot_22
double precision, external :: GETPEI

do j=1, Ngrid
    A1(j) = max(A1_in(j),  1.E-16)
    A2(j) = max(A2_in(j),  1.E-16)
    x_b(j) = x(j) + 0.5*dx
enddo

! Define explicit or implicit
SELECT CASE(imethod)
CASE(21: 23)
    theta = 1.
CASE(31: 33)
    theta = 0.5
END SELECT

! Map V*M, G, H on tildes
! So B = B + 1/G11*(bdot*V/W*N*x+(rdot-bdot)*H/M*x)                           for implicit
! So S = S - bdot*M*N*x*d/dx(V/W) - (rdot-bdot)/V*G*H*d/dx (x*V/G)      for implicit
! and P = P + bdot/W*d/dx (V*M*N*x*y)+(rdot-bdot)*x/G*d/dx (G*H*y)         for explicit

call GRID2GRID(1, x, V*M, Vtilde, Ngrid, 1)
 
! Compute source: Rsource = - 1/V d/dx (V*M*G11*R), Rsource is on main grid
do j=1, Ngrid
    Rsource1(j) = Vtilde(j)*G11(j)*R1(j)
    Rsource3(j) = 0.
enddo
call DERIV(x_b, x, 2, Rsource1, Rsource3, 1, Ngrid, 0)
do j=1, Ngrid
    Rsource3(j) = -Rsource3(j)/V(j)
    P1_new(j) = P1(j) + Rsource3(j)
    T12(j) = 625*GETPEI(j)
    Rsource2(j) = Vtilde(j)*G11(j)*R2(j)
    Rsource4(j) = 0.
enddo

call DERIV(x_b, x, 2, Rsource2, Rsource4, 1, Ngrid, 0)
do j=1, Ngrid
    Rsource4(j) = -Rsource4(j)/V(j)
    P2_new(j) = P2(j) + Rsource4(j)
    T21(j) = T12(j)
    S1_new(j) = S1_in(j) - T12(j)
    S2_new(j) = S2_in(j) - T21(j)
enddo

! Compute additionals
Pdot_11 = 0.
Pdot_21 = 0.
Pdot_12 = 0.
Pdot_22 = 0.

call DERIV(x, x_b, 1, V*M*N1*x*Y1O, ydummy, 1, Ngrid, 1)
call GRID2GRID(2, x_b, ydummy, Pdot_11, Ngrid, 0)
Pdot_11 = Pdot_11/W1

call DERIV(x, x_b, 1, V*M*N2*x*Y2O, ydummy, 1, Ngrid, 1)
call GRID2GRID(2, x_b, ydummy, Pdot_12, Ngrid, 0)
Pdot_12 = Pdot_12/W2

call DERIV(x, x_b, 1, GN*H1N*Y1O, ydummy, 1, Ngrid, 1)
call GRID2GRID(2, x_b, ydummy, Pdot_21, Ngrid, 0)
Pdot_21 = Pdot_21*x/GN

call DERIV(x, x_b, 1, GN*H2N*Y2O, ydummy, 1, Ngrid, 1)
call GRID2GRID(2, x_b, ydummy, Pdot_22, Ngrid, 0)
Pdot_22 = Pdot_22*x/GN

do j=1, Ngrid
    N1N(j) = GN(j)*H1N(j)
    N2N(j) = GN(j)*H2N(j)
    N1O(j) = GO(j)*H1O(j)
    N2O(j) = GO(j)*H2O(j)
    Vtilde(j) = Vtilde(j)*G11(j)
! P = P + bdot/W*d/dx (V*M*N*x*y)+(rdot-bdot)*x/G*d/dx (G*H*y)   for explicit
    P1_new(j) = P1_new(j) + bbdot*Pdot_11(j) + (rbdot - bbdot)*Pdot_21(j)
! P = P + bdot/W*d/dx (V*M*N*x*y)+(rdot-bdot)*x/G*d/dx (G*H*y)   for explicit
    P2_new(j) = P2_new(j) + bbdot*Pdot_12(j) + (rbdot - bbdot)*Pdot_22(j)
enddo

! Check boundary condition, note that Qbound = Qbound/G11(b) since G11 is absorbed in Vtilde
if (bctype(1) == 1) then
    f_bound1 = y1(Ngridb)
    NgridS(1) = Ngridb - 1
    Qbound1 = 0.0
else if (bctype(1) == 2) then
    f_bound1 = 0.0
    NgridS(1) = Ngridb
    Qbound1 = bcvalue(1)/G11(Ngridb)
endif

if (bctype(2) == 1) then
    f_bound2 = y2(Ngridb)
    NgridS(2) = Ngridb - 1
    Qbound2 = 0.0
else if (bctype(2) == 2) then
    f_bound2 = 0.0
    NgridS(2) = Ngridb
    Qbound2 = bcvalue(2)/G11(Ngridb)
endif 
      
! Define cd, or power law
SELECT CASE(imethod)

CASE(11, 21, 31)
    do j=1, Ngrid
        xi1(j) = dx*B1_in(j)/A1(j)
        fxi1(j) = 1. + 0.5*xi1(j)
        gxi1(j) = 1. - 0.5*xi1(j)
        xi2(j)  = dx*B2_in(j)/A2(j)
        fxi2(j) = 1. + 0.5*xi2(j)
        gxi2(j) = 1. - 0.5*xi2(j)
    enddo

CASE(12, 22, 32)
    do j=1, Ngrid
        xi1(j) = dx*B1_in(j)/A1(j)
        if (xi1(j) < -10.) then
            fxi1(j) = 0
        else if (xi1(j) < 0) then
            fxi1(j) = (1. + 0.1*xi1(j))**5.0
        else if (xi1(j) <= 10) then
            fxi1(j) = (1. - 0.1*xi1(j))**5.0 + xi1(j)
        else
            fxi1(j) = xi1(j)
        endif
        gxi1(j) = fxi1(j)-xi1(j)

        xi2(j) = dx*B2_in(j)/A2(j)
        if (xi2(j) < -10.) then
            fxi2(j) = 0
        else if (xi2(j) < 0) then
            fxi2(j) = (1. + 0.1*xi2(j))**5.0
        else if (xi2(j) <= 10) then
            fxi2(j) = (1. - 0.1*xi2(j))**5.0 + xi2(j)
        else
            fxi2(j) = xi2(j)
        endif
        gxi2(j) = fxi2(j) - xi2(j)
    enddo

CASE(13, 23, 33)
    do j=1, Ngrid
        xi1(j) = dx*B1_in(j)/A1(j)
        xi2(j) = dx*B2_in(j)/A2(j)
        if (xi1(j) == 0.) then
            fxi1(j) = 1.
            gxi1(j) = 1.
        else
            fxi1(j) = xi1(j)/(1. - exp(-xi1(j)))
            gxi1(j) = fxi1(j) - xi1(j)
        endif
        if (xi2(j) == 0.) then
            fxi2(j) = 1.
            gxi2(j) = 1.
        else
            fxi2(j) = xi2(j)/(1. - exp(-xi2(j)))
            gxi2(j) = fxi2(j) - xi2(j)
        endif
    enddo

END SELECT

call SOLVERIMP(dx, dt, &
    GN(1: Ngridb)*theta + GO(1: Ngridb)*(1 - theta), & 
    N1N(1: Ngridb), N1O(1: Ngridb), N2N(1: Ngridb), & 
    N2O(1: Ngridb), V(1: Ngridb), Vtilde(1: Ngridb), & 
    A1(1: Ngridb), A2(1: Ngridb), & 
    fxi1(1: Ngridb), fxi2(1: Ngridb), gxi1(1: Ngridb), & 
    gxi2(1: Ngridb), S1_new(1: Ngridb), & 
    S2_new(1: Ngridb), & 
    P1_new(1: Ngridb), P2_new(1: Ngridb), & 
    T12(1: Ngridb), T21(1: Ngridb), & 
    Ngridb, NgridS, theta, Y1O(1: Ngridb), & 
    Y2O(1: Ngridb), & 
    f_bound1, Qbound1, f_bound2, Qbound2, & 
    bctype, y1(1: Ngridb), y2(1: Ngridb))

! Compute flux from solution, flux is on shifted grid
do j=1, Ngridb-1
    Q1_out(j) = G11(j)*(-A1(j)/dx*(fxi1(j)*y1(j+1) - gxi1(j)*y1(j)) + R1(j))
    Q2_out(j) = G11(j)*(-A2(j)/dx*(fxi2(j)*y2(j+1) - gxi2(j)*y2(j)) + R2(j))
enddo
if (bctype(1) == 1) then
    call EXTRAP(x(1: NgridS(1)), Q1_out(1: NgridS(1)), x(Ngridb), & 
        NgridS(1), Q1_out(Ngridb), 1, NgridS(1))
else if (bctype(1) == 2) then
! Restore G11 in Qbound
    Q1_out(Ngridb) = Qbound1*G11(Ngridb)
endif

if (bctype(2) == 1) then
    call EXTRAP(x(1: NgridS(2)), Q2_out(1: NgridS(2)), x(Ngridb), & 
        NgridS(2), Q2_out(Ngridb), 1, NgridS(2))
else if (bctype(2) == 2) then
! Restore G11 in Qbound
    Q2_out(Ngridb) = Qbound2*G11(Ngridb)
endif

return
end subroutine RUNEQTIMP

!---------------------------------------------------------------------
subroutine SOLVERIMP(dx, dt, G, N1N, N1O, N2N, N2O, V, Vtilde, & 
    A1, A2, fxi1, fxi2, gxi1, gxi2, S1, S2, P1, P2, T12, T21, &
    Ngrid, NgridS, theta, Y1O, Y2O, f_bound1, Qbound1, &
    f_bound2, Qbound2, bctype, y1, y2)
!---------------------------------------------------------------------
! Build up matrices to be passed ot TRIDIAG2

implicit none

integer, intent(in) :: Ngrid, NgridS(2), bctype(2)
double precision, intent(in) :: dx, dt, theta, f_bound1, f_bound2, &
    Qbound1, Qbound2
double precision, intent(in), dimension(Ngrid) :: G, &
    N1N, N1O, N2N, N2O, V, Vtilde, A1, A2, fxi1, fxi2, &
    gxi1, gxi2, S1, S2, P1, P2, T12, T21, Y1O, Y2O
double precision, intent(out), dimension(Ngrid) :: y1, y2

integer :: j
double precision :: dt_dx2
double precision, dimension(Ngrid) :: &
    AA1, BB1, CC1, RR1, TT1, &
    AA2, BB2, CC2, RR2, TT2

!    1/G d/dt (N1*y1) + 1/V d/dx (Vtilde*(-A1/dx*(fxi1, gxi1, ntilde1))) = S1*y1+P1+T1*y2
!    1/G d/dt (N2*y2) + 1/V d/dx (Vtilde*(-A2/dx*(fxi2, gxi2, ntilde2))) = S2*y2+P2+T2*y1

! implicit
j = 1
dt_dx2 = dt/dx**2

AA1(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j) * A1(j)*fxi1(j) * theta
BB1(j) = N1N(j) + theta * ( -dt*G(j)*S1(j) - dt_dx2 * G(j)/V(j) * &
    Vtilde(j) * (A1(j)*(-gxi1(j))) )
TT1(j) = theta*(-dt*G(j)*T12(j))
CC1(j) = 0.0      
RR1(j) = N1O(j)*Y1O(j) + G(j)*P1(j)*dt + (1 - theta) * &
    ( dt*G(j)*S1(j)*Y1O(j) + dt_dx2 * G(j)/V(j) * &
    Vtilde(j) * (A1(j)*(fxi1(j)*Y1O(j+1) - gxi1(j)*Y1O(j))) )

AA2(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j) * A2(j)*fxi2(j) * theta
BB2(j) = N2N(j) + theta * ( -dt*G(j)*S2(j) - dt_dx2 * G(j)/V(j) * &
    Vtilde(j) * (A2(j)*(-gxi2(j))) )
TT2(j) = theta*(-dt*G(j)*T21(j))
CC2(j) = 0.0      
RR2(j) = N2O(j)*Y2O(j) + G(j)*P2(j)*dt + (1 - theta)* &
    ( dt*G(j)*S2(j)*Y2O(j) + dt_dx2 * G(j)/V(j) * &
    Vtilde(j) *(A2(j)*(fxi2(j)*Y2O(j+1) - gxi2(j)*Y2O(j))) )

do j=2, NgridS(1)-1
    AA1(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j) * A1(j)*fxi1(j) * theta
    BB1(j) = N1N(j) + theta * (  -dt*G(j)*S1(j) - dt_dx2 * G(j)/V(j) * &
        ( Vtilde(j)  *(A1(j)  *(-gxi1(j))) &
         -Vtilde(j-1)*(A1(j-1)*(fxi1(j-1))) )  )
    TT1(j) = theta*(-dt*G(j)*T12(j))
    CC1(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j-1)*A1(j-1)*gxi1(j-1) * theta
    RR1(j) = N1O(j)*Y1O(j) + G(j)*P1(j)*dt + (1 - theta) * &
        (  dt*G(j)*S1(j)*Y1O(j) + dt_dx2 * G(j)/V(j) * &
        (Vtilde(j)  *(A1(j)  *(fxi1(j)*Y1O(j+1) - gxi1(j)  *Y1O(j))) - &
        Vtilde(j-1)*(A1(j-1)*(fxi1(j-1)*Y1O(j) - gxi1(j-1)*Y1O(j-1))))  )
enddo

do j=2, NgridS(2)-1
    AA2(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j) *A2(j)*fxi2(j) * theta
    BB2(j) = N2N(j) + theta * ( -dt*G(j)*S2(j) - dt_dx2 * G(j)/V(j) * &
        (Vtilde(j)  *(A2(j)  *(-gxi2(j))) &
        -Vtilde(j-1)*(A2(j-1)*(fxi2(j-1)))) )
    TT2(j) = theta*(-dt*G(j)*T21(j))
    CC2(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j-1) *A2(j-1)*gxi2(j-1) * theta
    RR2(j) = N2O(j)*Y2O(j) + G(j)*P2(j)*dt + (1 - theta) * &
        (  dt*G(j)*S2(j)*Y2O(j) + dt_dx2 * G(j)/V(j) * &
        (Vtilde(j)  *(A2(j)  *(fxi2(j)*Y2O(j+1) - gxi2(j)  *Y2O(j))) - &
         Vtilde(j-1)*(A2(j-1)*(fxi2(j-1)*Y2O(j) - gxi2(j-1)*Y2O(j-1))))  )
enddo

! bc type 1 (yb)
if (bctype(1) == 1) then
    j = Ngrids(1)
    AA1(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j) *A1(j)*fxi1(j) * theta
    BB1(j) = N1N(j) + theta * ( -dt*G(j)*S1(j) - dt_dx2 * G(j)/V(j) * &
        (Vtilde(j)  *(A1(j)  *(-gxi1(j))) - &
         Vtilde(j-1)*(A1(j-1)*(fxi1(j-1)))) )
    TT1(j) = theta*(-dt*G(j)*T12(j))
    CC1(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j-1) * A1(j-1)*gxi1(j-1) * theta
    RR1(j) = N1O(j)*Y1O(j) + G(j)*P1(j)*dt + (1 - theta)* &
        (  dt*G(j)*S1(j)*Y1O(j) + dt_dx2 * G(j)/V(j) * &
        (Vtilde(j)  *(A1(j)  *(fxi1(j)*Y1O(j+1) - gxi1(j)  *Y1O(j))) - &
         Vtilde(j-1)*(A1(j-1)*(fxi1(j-1)*Y1O(j) - gxi1(j-1)*Y1O(j-1))))  )
endif

! bc type 2 (Qb)
if (bctype(1) == 2) then
    j = Ngrids(1)
    AA1(j) = 0.0      
    BB1(j) = N1N(j) + theta * ( -dt*G(j)*S1(j) - dt_dx2 * G(j)/V(j) * &
        (Vtilde(j)  *(A1(j)  *(-0.0*gxi1(j))) - &
         Vtilde(j-1)*(A1(j-1)*(fxi1(j-1)))) )
    TT1(j) = theta*(-dt*G(j)*T12(j))
    CC1(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j-1) * A1(j-1)*gxi1(j-1) * theta
    RR1(j) = N1O(j)*Y1O(j) + G(j)*P1(j)*dt + (1 - theta)* &
        (  dt*G(j)*S1(j)*Y1O(j) + dt_dx2 * G(j)/V(j) * &
        (-Vtilde(j-1)*(A1(j-1)*(fxi1(j-1)*Y1O(j) - gxi1(j-1)*Y1O(j-1))))  ) + &
        dt_dx2 * G(j)/V(j) * (-Vtilde(j)*Qbound1*dx)
endif

! bc type 1 (yb)
if (bctype(2) == 1) then
    j = Ngrids(2)
    AA2(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j)*A2(j)*fxi2(j) * theta
    BB2(j) = N2N(j) + theta * ( -dt*G(j)*S2(j) - dt_dx2 * G(j)/V(j) * &
        (Vtilde(j)  *(A2(j)  *(-gxi2(j))) - &
         Vtilde(j-1)*(A2(j-1)*(fxi2(j-1)))) )
    TT2(j) = theta*(-dt*G(j)*T21(j))
    CC2(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j-1)*A2(j-1)*gxi2(j-1) * theta
    RR2(j) = N2O(j)*Y2O(j) + G(j)*P2(j)*dt + (1 - theta) * &
        (  dt*G(j)*S2(j)*Y2O(j) + dt_dx2 * G(j)/V(j) * &
        (Vtilde(j)  *(A2(j)  *(fxi2(j)  *Y2O(j+1) - gxi2(j)*Y2O(j))) - &
         Vtilde(j-1)*(A2(j-1)*(fxi2(j-1)*Y2O(j) - gxi2(j-1)*Y2O(j-1))))  )
endif
  
! bc type 2 (Qb)
if (bctype(2) == 2) then
    j = Ngrids(2)
    AA2(j) = 0.0
    BB2(j) = N2N(j) + theta * ( -dt*G(j)*S2(j) - dt_dx2 * G(j)/V(j) * &
        (Vtilde(j)  *(A2(j)  *(-0.0*gxi2(j))) - &
         Vtilde(j-1)*(A2(j-1)*(fxi2(j-1))))  )
    TT2(j) = theta*(-dt*G(j)*T21(j))
    CC2(j) = -dt_dx2 * G(j)/V(j)*Vtilde(j-1)*A2(j-1)*gxi2(j-1) * theta
    RR2(j) = N2O(j)*Y2O(j) + G(j)*P2(j)*dt + (1 - theta)* &
        (  dt*G(j)*S2(j)*Y2O(j) + dt_dx2 * G(j)/V(j) * &
        (-Vtilde(j-1)*(A2(j-1)*(fxi2(j-1)*Y2O(j) - gxi2(j-1)*Y2O(j-1))))  ) + &
        dt_dx2 * G(j)/V(j)*(-Vtilde(j)*Qbound2*dx)
endif

! Main call to tridiagonal solver
call TRIDIAG2(AA1, BB1, CC1, RR1, y1, TT1, Ngrid, bctype, &
    f_bound1, AA2, BB2, CC2, RR2, TT2, y2, f_bound2)

return
end subroutine SOLVERIMP

!---------------------------------------------------------------------
subroutine TRIDIAG2(A1, B1, C1, R1, f1, T1, Ngrid, bctype, &
    f_bound1, A2, B2, C2, R2, T2, f2, f_bound2)
!---------------------------------------------------------------------
! Inputs: A(1:Ngrid), B(1:Ngrid), C(1:Ngrid), R(1:Ngrid)
! Outputs: f(1:Ngrid)
! 
! Provides solution of the system:
!
!   A1j f1j+1  + B1j f1j  + C1j f1j-1 + T1j f2j = R1j
!   A2j f2j+1  + B2j f2j  + C2j f2j-1 + T2j f1j = R2j
!
!   where j = 1..Ngrid
!
!   bctype = 1  -> given f_NA1
!   bctype = 2  -> given Gamma_NA1
!
!  implicit
!---------------------------------------------------------------------

implicit none

integer, intent(in) :: Ngrid, bctype(2)
double precision, intent(in) :: f_bound1, f_bound2
double precision, intent(in), dimension(Ngrid) :: &
    A1, B1, C1, R1, T1, &
    A2, B2, C2, R2, T2
double precision, intent(out), dimension(Ngrid) :: f1, f2

integer :: j, k
double precision :: detjm1, k1, k2, k3, k4, k31, k32
double precision, dimension(Ngrid-1) :: alpha, beta, gamma, &
    delta, epsilon, theta

! Implicit

j = 1
alpha(j)   = -B1(j)/A1(j)
beta(j)    =  R1(j)/A1(j)
epsilon(j) = -T1(j)/A1(j)      
gamma(j)   = -B2(j)/A2(j)
delta(j)   =  R2(j)/A2(j)
theta(j)   = -T2(j)/A2(j)
   
do j=2, Ngrid-1
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
j = Ngrid
if (bctype(1) == 1 .and. bctype(2) == 1) then
    detjm1 = (alpha(j-1)*gamma(j-1) - theta(j-1)*epsilon(j-1))
    f1(j-1) =  gamma(j-1)/detjm1*(f_bound1 - beta(j-1)) - &
             epsilon(j-1)/detjm1*(f_bound2 - delta(j-1))
    f2(j-1) = -theta(j-1)/detjm1*(f_bound1 - beta(j-1)) + &
               alpha(j-1)/detjm1*(f_bound2 - delta(j-1))
    do k=1, Ngrid-2
        j = Ngrid - 2 - k + 1
        detjm1 = (alpha(j)*gamma(j) - theta(j)*epsilon(j))
        f1(j) =  gamma(j)/detjm1*(f1(j+1) - beta(j)) - &
               epsilon(j)/detjm1*(f2(j+1) - delta(j))
        f2(j) = -theta(j)/detjm1*(f1(j+1) - beta(j)) + &
                 alpha(j)/detjm1*(f2(j+1) - delta(j))
    enddo
    f1(Ngrid) = f_bound1
    f2(Ngrid) = f_bound2
endif

return
end subroutine TRIDIAG2
