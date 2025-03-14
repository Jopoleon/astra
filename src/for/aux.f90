integer function IFSTEP(IFCONV, updwno)

! IFCONV dummy parameter
! Input
!       LEQ(1)  LEQ(2)  LEQ(3)  LEQ(4)  LEQ(5)  LEQ(6-9)
!  NE  TE  TI  CU equil    free
!       LEQ(10) LEQ(11) LEQ(12), etc
!        F0  F1  F2
!       LEQ(j) - usage of j-th equation in the model
!        =-1 no equation (default)
!        = 0 type: AS
!        = 1 type: EQ
!        = 2 type: FU
!        = 3 type: heat conductivity flux + heat convection
! Returned value
!  1 - time step made, normal exit
!  0 - time step will be repeated
!---------------------------------------------------------------------

use const_inc, only: TAUINC, DELVAR, TAU, TAUPRP, TAUMIN, TAUMAX, &
    DTOUT, DPOUT, NA, NB1, LEQ, NSTEPS, ROC, ROCO, &
    FTN, FTO
use status_inc, only: NEO, NIO, TEO, TIO, FJO, FPO, VRO, UPARO, &
    NE, NI, TE, TI, FJ, FP, VR, UPAR
use debugger, only: markloc

implicit none

integer, intent(in) :: IFCONV
double precision, intent(in) :: updwno

integer :: j, jj
double precision :: CTAU, TAUO, YY, TAUN

call markloc('IFSTEP')

CTAU = 1./TAUINC

do j = 1, NA
    if (LEQ(1) > 0) CTAU = MAX(CTAU, ABS(NEO(j)/NE(j) - 1.)/DELVAR)
    if (LEQ(2) > 0) CTAU = MAX(CTAU, ABS(TEO(j)/TE(j) - 1.)/DELVAR)
    if (LEQ(3) > 0) CTAU = MAX(CTAU, ABS(TIO(j)/TI(j) - 1.)/DELVAR)
    do jj=0, 9
        if (LEQ(jj+10) > 0) then
            YY = 0.5*(abs(FJO(j, jj)) + abs(FJ(j, jj)))
            if (YY < 1.E-6) then  ! Allow zero FJ
                YY = abs(FJO(j, jj) - FJ(j, jj))
            else
                YY = abs(FJO(j, jj) - FJ(j, jj))/YY
            endif
            CTAU = MAX(CTAU, YY/DELVAR)
        endif
    enddo
enddo

TAUO   = TAU
TAUPRP = TAUO
TAUN   = TAU
TAUN   = MIN(TAUMAX, TAUN/CTAU, DTOUT, DPOUT)

TAU = MAX(TAUMIN, TAUN)   ! due to DELVAR & TAUINC

if (TAU >= TAUO) then
    IFSTEP = 1
    NSTEPS = NSTEPS+1
    return  ! -> proceed to the next time step
endif

! Reset all functions F(t) to old values

do j=1, NB1
    NE(j) = NEO(j)
    NI(j) = NIO(j)
    TE(j) = TEO(j)
    TI(j) = TIO(j)
    FP(j) = FPO(j)
    VR(j) = VRO(j)
    UPAR(j) = UPARO(j)
    do jj=0, 9
        FJ(j, jj) = FJO(j, jj)
    enddo
enddo

! Reset ROC
ROC = ROCO
FTO = FTN

call CUOFP  ! Restore CU & MU, to be checked after VR is done EFable
IFSTEP = 0

return
end function IFSTEP

!---------------------------------------------------------------------
integer function IFTREQ(YACC)
!---------------------------------------------------------------------
!    Control parameter:
! NITREQ == 1 iterations switched off (timing VR(t) downshifted)
!       IFTREQ == 2 is returned
! NITREQ >= 2 iterations switched on  (VR(t) properly included)
!
!    Input:
! YACC - requested accuracy
! 
! ITREQ = 0 is set in BLOCKDATA, reset in stepup, incremented in METRIC
! ITREQ == 1 save metric parameters (G11, G22, VR, ROC) and exit
!            IFTREQ == 2 is returned
!
! ITREQ >= 1 & ITREQ <= NITREQ - Mark time evolution section
!
!---------------------------------------------------------------------
!    Returned value:
! IFTREQ == 0 - No convergence
! IFTREQ == 1 - Convergence condition is achieved
! IFTREQ == 2 - No convergence and maximum iteration number is achieved
!---------------------------------------------------------------------

use parameter_inc, only: NRD
use status_inc, only: G11, G22, VR, FP
use const_inc, only: NITREQ, IPART, ITREQ, ROC, PSIEXT, PSPLEX, NA1, NB1
use debugger, only: markloc

implicit none

integer, parameter :: ITREQMIN=1, ITREQMAX=200
double precision, intent(in) :: YACC

integer :: j, NTREQ

double precision :: Y1, Y2, YV, YI, YR, YER(ITREQMAX), YPSE, YPSP
double precision, dimension(NRD) :: YWA, YWB, YWC, YWD

save YER, YR, YWA, YWB, YWC, YWD
save YPSE, YPSP

call markloc('IFTREQ')

NTREQ = nint(NITREQ)
! Disable NITREQ setting for the initial phase
if (IPART == 1) NTREQ = ITREQMAX
IFTREQ = 2
if (NTREQ == 1) return ! Iterations are off, no check

if (NTREQ > ITREQMAX) then  ! 
     write(*,'(2A,I5)') " >>> Warning >>> Max number of transport/equilibrium", &
         " iterations is reduced to ", ITREQMAX
     NTREQ = ITREQMAX
endif
IFTREQ = 0

if (ITREQ /= 0) then ! From 2nd iteration

    if (ITREQ >= ITREQMAX) then
        write(*,'(2A,2I5)') " >>> Warning >>> Max number of initialization ", &
             " iterations is overwhelmed, stopping ", ITREQ, ITREQMAX
        stop
    endif

! Convergence check:  enabled if NTREQ > 1
    YR   = abs(YR/ROC      - 1.d0)
    YPSE = abs(YPSE/PSIEXT - 1.d0)
    YPSP = abs(YPSP/PSPLEX - 1.d0)
    Y1 = 0.d0
    Y2 = 0.d0
    YV = 0.d0
    YI = 0.d0
    do j=1, NA1  ! Analize NA1new vs NA1old
        Y1 = max(Y1, abs(YWA(j)/G11(j) - 1.d0))
        Y2 = max(Y2, abs(YWB(j)/G22(j) - 1.d0))
        YV = max(YV, abs(YWC(j)/VR(j) - 1.d0))
        YI = max(YI, abs(YWD(j)/FP(j) - 1.d0))
    enddo
    if (Y1 > max(Y2, YV, YR, YI)) then
        YER(ITREQ) = Y1
        j = 1
    endif
    if (Y2 > max(Y1, YV, YR, YI)) then
        YER(ITREQ) = Y2
        j = 2
    endif
    if (YV > max(Y1, Y2, YR, YI)) then
        YER(ITREQ) = YV
        j = 3
    endif
    if (YR > max(Y1, Y2, YV, YI)) then
        YER(ITREQ) = YR
        j = 4
    endif
    if (YI > max(Y1, Y2, YV, YR)) then
        YER(ITREQ) = YI
        j = 5
    endif

! Check convergence:
! Note! the tolerance, YACC(alias ATREQ), should be >> than ACEQLB (1.d-6)
    IFTREQ = 0
    if (ITREQ == NTREQ    ) IFTREQ = 2 ! Go to next t-step anyway
    if (YACC >= YER(ITREQ)) IFTREQ = 1 ! Converged -> next t-step
    if (ITREQ < ITREQMIN  ) IFTREQ = 0 ! force another iteration if below minimum

    if (IFTREQ > 0) return ! converged

endif

ITREQ = ITREQ + 1 

YR   = ROC
YPSE = PSIEXT
YPSP = PSPLEX
do j=1, NB1   ! Store some metric data
    YWA(j) = G11(j)
    YWB(j) = G22(j)
    YWC(j) = VR(j)
    YWD(j) = FP(j)
enddo

return
end function IFTREQ

!---------------------------------------------------------------------
subroutine OLDNEW

use status_inc, only: NE, TE, TI, FP, VR, FJO, UPAR, &
    UPS0, UPS1, UPS2, NEO, TEO, TIO, FPO, VRO, FJ, &
    UPARO, UPS0O, UPS1O, UPS2O
use const_inc, only: NB1, BTN, FTO, FTN, ROCO, ROC, BTOR, BTN
use debugger, only: markloc

implicit none

integer :: j, jj

call markloc('OLDNEW')

do j=1, NB1
    NEO(j) = NE(j)
    TEO(j) = TE(j)
    TIO(j) = TI(j)
    FPO(j) = FP(j)
    VRO(j) = VR(j)
    UPARO(j) = UPAR(j)
    UPS0O(j) = UPS0(j)
    UPS1O(j) = UPS1(j)
    UPS2O(j) = UPS2(j)
    do jj=0, 9
        FJO(j, jj) = FJ(j, jj)
    enddo
enddo

BTN = BTOR
FTN = FTO
ROCO = ROC

return
end subroutine OLDNEW

!---------------------------------------------------------------------
subroutine SMOOTH(ALFA, n_in, y_in, x_in, n_out, y_out, x_out)
!---------------------------------------------------------------------
!  Subroutine minimizes the value of functional
!  INTEGRAL(alfa*P(x)*(dU/dx)**2+(U-F)**2)*dx, 
!  where y_in(n_in) is a function, given on the grid x_in(n_in)
! P(x) is equal to unit now
! ALFA=alfa<<0.01*x_in(n_in)**2 is regularizator
! n_in - number of old grid points
! N=<NRD - number of new grid points
! 0<=x_in(n_in) - old grid |     both grids are arbitrary
! 0<=x_out(N)  - new grid |     but x_in(n_in)=x_out(N)
! y_in(n_in) - origin function, given on the grid x_in(n_in)
! y_out(N) - smoothed function on grid x_out(N)
!  The result is function y_out(x_out), given on the new grid
! with additional conditions:
! dy_out/dx(x=0)=0 - cylindrical case and
! y_out(x_out(N))=y_in(x_in(n_in))
!---------------------------------------------------------------------

use parameter_inc, only: NRD
use debugger, only: astra_stop

implicit none

integer, intent(in) :: n_in, n_out
double precision, intent(in) :: ALFA, x_in(*), y_in(*), x_out(*)
double precision, intent(out) :: y_out(*)

integer :: I, J
double precision :: YF, YX, YP, YQ, YD, FJ, P(NRD), dx
character(len=132) :: err_msg

if (n_out > NRD .or. n_in <= 0) then
    write(err_msg, '(A, 1X, i, 1X, i)') ' >>> SMOOTH: array is out of limits', NRD, n_out
    call err_catch_a
    call astra_stop(err_msg)
endif
if (n_in == 1) then
    do j=1, n_out
        y_out(j) = y_in(1)
    enddo
    return
endif
if (n_in == 2) then
    do j=1, n_out
        y_out(j) = (y_in(2)*(x_out(j) - x_in(1)) - y_in(1)*(x_out(j) - x_in(2)))/(x_in(2) - x_in(1))
    enddo
    return
endif
if (n_out < 2) then
    call err_catch_a
    call astra_stop(' >>> SMOOTH: no output grid is provided')
endif
if (abs(x_in(n_in) - x_out(n_out)) > x_out(n_out)/n_out) then
    write(*, *) '>>> SMOOTH: grids are not aligned'
    write(*, '(1A23, I4, F8.4)') '     Old grid size/edge', n_in , x_in(n_in)
    write(*, '(1A23, I4, F8.4)') '     New grid size/edge', n_out, x_out(n_out)
    call err_catch_a
    call astra_stop
endif
do j=2, n_out
    dx = x_out(j) - x_out(j-1)
    if (dx <= 0.) then
        write(*, *)'>>> SMOOTH: new grid is not increasing monotonically'
        write(*, '(A, I4, A, F8.4)')'Node ', j-1, '   Value', x_out(j-1)
        write(*, '(A, I4, A, F8.4)')'Node ', j  , '   Value', x_out(j)
        call err_catch_a
        call astra_stop
    endif
    P(j) = ALFA/dx/x_in(n_in)**2
enddo
P(1)  = 0.
y_out(1) = 0.
i = 1
YF = (y_in(2) - y_in(1))/(x_in(2) - x_in(1))
YX = 2./(x_out(2) + x_out(1))
YP = 0.
YQ = 0.
do j=1, n_out-1
    if (x_in(i) <= x_out(j)) then
        do
            i = i + 1
            i = min(i, n_in)
            if (i == n_in .or. x_in(i) >= x_out(j)) EXIT
        enddo
        YF = (y_in(i) - y_in(i-1))/(x_in(i) - x_in(i-1))
    endif
    FJ = y_in(i) + YF*(x_out(j) - x_in(i))
    YD = 1. + YX*(YP + P(j+1))
    P(j) = YX*P(j+1)/YD
    y_out(j) = (FJ + YX*YQ)/YD
    if (j /= n_out-1) then
        YX = 2./(x_out(j+2) - x_out(j))
        YP = (1. - P(j))*P(j+1)
        YQ = y_out(j)*P(j+1)
    endif
enddo

y_out(n_out) = y_in(n_in)
do j=n_out-1, 1, -1
    y_out(j) = P(j)*y_out(j+1) + y_out(j)
enddo

return
end subroutine SMOOTH

!---------------------------------------------------------------------
subroutine SMAP(ALFA, n_in, x_in, n_out, x_out, y_out)
! Similar to SMOOTH but the same array, y_out, is used for input and output

use parameter_inc, only: NRD

implicit none

integer, intent(in) :: n_in, n_out
double precision, intent(in) :: ALFA, x_in(*), x_out(*)
double precision, intent(inout) :: y_out(*)

integer :: J
double precision :: P(NRD)

call SMOOTH(ALFA, n_in, y_out, x_in, n_out, P, x_out)
do j=1, n_out
    y_out(j) = P(j)
enddo

return
end subroutine SMAP

!---------------------------------------------------------------------
subroutine TRANSF(n_in, y_in, x_in, n_out, y_out, x_out)
!---------------------------------------------------------------------
!  The subroutine transmits a function y_in(1:n_in)
!  from an arbitrary grid x_in(1:n_in) to a functon y_out(1:n_out)
!  on an another arbitrary grid x_out(1:n_out) by the method
!  of quadratic interpolation.
! Input: n_in, y_in(1:n_in), x_in(1:n_in), n_out, x_out(1:n_out)
! Output: y_out(1:n_out)
!---------------------------------------------------------------------
!       Note: [x_in(n_in)-x_in(1)]*[x_out(n_out)-x_out(1)] must be > 0
!---------------------------------------------------------------------

implicit none

integer, intent(in) :: n_in, n_out
double precision, intent(in) , dimension(*) :: x_in, y_in, x_out
double precision, intent(out), dimension(*) :: y_out

integer :: I, J
double precision :: YF1, YF2, YF3

I = 1
YF1 = y_in(1)/((x_in(1) - x_in(2))*(x_in(1) - x_in(3)))
YF2 = y_in(2)/((x_in(2) - x_in(1))*(x_in(2) - x_in(3)))
YF3 = y_in(3)/((x_in(3) - x_in(1))*(x_in(3) - x_in(2)))
do j=1, n_out
    if(2.*x_out(j) > x_in(I+1) + x_in(I+2)) then
        do
            I = I + 1
            I = min(I, n_in - 2)
            if (2.*x_out(j) <= x_in(I+1) + x_in(I+2) .or. I >= n_in-2) EXIT
        enddo
        YF1 = y_in(I  )/((x_in(I  ) - x_in(I+1))*(x_in(I  ) - x_in(I+2)))
        YF2 = y_in(I+1)/((x_in(I+1) - x_in(I  ))*(x_in(I+1) - x_in(I+2)))
        YF3 = y_in(I+2)/((x_in(I+2) - x_in(I  ))*(x_in(I+2) - x_in(I+1)))
    endif
    y_out(j) = YF1*(x_out(j) - x_in(I+1))*(x_out(j) - x_in(I+2)) + &
            YF2*(x_out(j) - x_in(I  ))*(x_out(j) - x_in(I+2)) + &
            YF3*(x_out(j) - x_in(I  ))*(x_out(j) - x_in(I+1))
enddo

return
end subroutine TRANSF

!---------------------------------------------------------------------
subroutine QMAP(n_in, x_in, n_out, x_out, y_out)
! Similar to TRANSF but uses the same array, y_out, for input and output.

use parameter_inc, only: NRD

implicit none

integer, intent(in) :: n_in, n_out
double precision, intent(in), dimension(*) :: x_in, x_out
double precision, intent(inout), dimension(*) :: y_out

integer :: j
double precision :: FN(NRD)

call TRANSF(n_in, y_out, x_in, n_out, FN, x_out)

do j=1, n_out
    y_out(j) = FN(j)
enddo

return
end subroutine QMAP

!---------------------------------------------------------------------
subroutine CHEBFT(n_in, y_in, x_in, NCHCF, CHEBCF)
!---------------------------------------------------------------------
! The subroutine returns NCHCF coefficients of a Chebyshev
! polynomial fit to the function y_in(1:n_in) given as a
! function of any "radial" variable on the grid x_in(1:n_in)
! Example:
! call CHEBFT(NA1,TE,FP,5,CHOUT)
! out = PFITN(FP,CHOUT,5)
!-----------------------------------------------------------------------

use debugger, only: markloc, astra_stop

implicit none

integer, parameter :: NMAX=10
double precision, parameter :: PI=3.141592654

integer, intent(in) :: n_in, NCHCF
double precision, intent(in) , dimension(n_in)  :: y_in, x_in
double precision, intent(out), dimension(NCHCF) :: CHEBCF

integer :: K, J
double precision :: YA, YB, FAC, YD, SUM, PIOVN, BMA, BPA
double precision, dimension(NMAX) :: YC, YF

call markloc('CHEBFT')

if (NCHCF > NMAX) call astra_stop('>>> Chebyshev fit error: too high power')

YA = max(x_in(1), x_in(n_in))
YB = min(x_in(1), x_in(n_in))
PIOVN = PI/NCHCF

! The following two lines require that at the boundary
!     the fit concides with the original function
YD = COS(0.5*PIOVN)
YA = (2.*YA - YB*(1. - YD))/(1. + YD)

BMA = 0.5*(YB - YA)
BPA = 0.5*(YB + YA)
do K=1,NCHCF
    YC(K) = BPA + BMA*COS(PIOVN*(K - 0.5))
enddo

call TRANSF(n_in, y_in, x_in, NCHCF, YF, YC)

FAC = 2./NCHCF
do J=1, NCHCF
    SUM = 0.
    do K=1, NCHCF
        SUM = SUM + YF(K)*COS(PIOVN*(K - 0.5)*(J - 1.))
    enddo
    YC(J) = FAC*SUM
enddo

call CHEBPC(YC, CHEBCF, YF, NCHCF)

FAC = 1./BMA
do J=2, NCHCF
    CHEBCF(J) = CHEBCF(J)*FAC
    FAC = FAC/BMA
enddo
do J=1, NCHCF-1
    do K=NCHCF-1, J, -1
        CHEBCF(K) = CHEBCF(K) - BPA*CHEBCF(K+1)
    enddo
enddo

return
end subroutine CHEBFT

!---------------------------------------------------------------------
subroutine CHEBPC(C, D, F, N)
!---------------------------------------------------------------------
! F(1:N) array for internal use
!
!   Input: N number of coefficients
!  C(1:N) 
!   Output: D(1:N) array of Chebyshev coefficients
!---------------------------------------------------------------------

implicit none

integer, intent(in) :: N
double precision, intent(in) , dimension(N) :: C
double precision, intent(out), dimension(N) :: D, F

integer :: J, K
double precision :: SV

D = 0.
F = 0.

D(1) = C(N)
do J=N-1, 2, -1
    do K=N-J+1, 2, -1
        SV = D(K)
        D(K) = 2.*D(K-1) - F(K)
        F(K) = SV
    enddo
    SV = D(1)
    D(1) = -F(1) + C(J)
    F(1) = SV
enddo
do J=N, 2, -1
    D(J) = D(J-1) - F(J)
enddo
D(1) = -F(1) + 0.5*C(1)

return
end subroutine CHEBPC

!---------------------------------------------------------------------
double precision function PFITN(x, cp, n_order)
!---------------------------------------------------------------------
! Nth order Polynomial FIT to a function f(a)
! CP are the N polynomial coefficients given
! PFITN = f(a) = \Sum {CP_j * a^(j-1)} ;   1<=j<=N
!---------------------------------------------------------------------

implicit none

integer, intent(in) :: n_order
double precision, intent(in) :: x, cp(n_order)

integer :: j

PFITN = CP(n_order)
do j=n_order-1, 1, -1
    PFITN = PFITN*x + CP(j)
enddo

return
end function PFITN
