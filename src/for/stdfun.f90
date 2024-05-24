!---------------------------------------------------------------------
integer function N_95_POS(i)
! Returns the radial integer index of the psi_95 position, at the left
! i in input is a dummy integer, so call it as N_95_POS(0)

use const_inc, only: NA1
use status_inc, only: FP_NORM

implicit none

integer, intent(in) :: i

integer :: j

N_95_POS = 1
do j=1, NA1
    if (FP_NORM(j) <= 0.95) N_95_POS = j
enddo

return
end function N_95_POS

!---------------------------------------------------------------------
double precision function V_95_POS(Yin)
! Returns the value of array Y(na1) at the 0.95 psi position

use const_inc, only: NA1, PSIAX
use status_inc, only: FP

implicit none

double precision, intent(in), dimension(NA1) :: Yin

integer :: j
integer, external :: N_95_POS
double precision :: r1, r2, y1, y2
double precision, dimension(NA1) :: rhop

j = N_95_POS(0)

rhop = (FP(1:NA1) - PSIAX)/(FP(NA1) - PSIAX)

y1 = Yin(j)
y2 = Yin(j+1)
r1 = rhop(j)
r2 = rhop(j+1)

V_95_POS = (y1*(r2 - 0.95) + y2*(0.95 - r1))/(r2 -r1)

return
end function V_95_POS

!---------------------------------------------------------------------
double precision function RFA(YA)
! rho=f(a) "rho" at a given radius "a"

use const_inc, only: NA1
use status_inc, only: AMETR, RHO
use numerical_tools, only: QUADIN

implicit none

double precision, intent(in) :: YA

RFA = QUADIN(NA1, AMETR, RHO, YA)

return
end function RFA

!---------------------------------------------------------------------
double precision function RFAN(YAN)
! rho=f(r/a) "rho" at a given radius "r/a"

use const_inc, only: ABC

implicit none

double precision, intent(in) :: YAN
double precision :: RFA

RFAN = RFA(YAN*ABC)

return
end function RFAN

!---------------------------------------------------------------------
double precision function XFA(YA)
! x=f(a) "x"=rho/roc at a given radius "a"

use const_inc, only: ROC

implicit none

double precision, intent(in) :: YA
double precision :: RFA

XFA = RFA(YA)/ROC

return
end function XFA

!---------------------------------------------------------------------
double precision function XFAN(YAN)
! rho=f(a) "rho" at a given radius "a"

use const_inc, only: ABC

implicit none

double precision, intent(in) :: YAN
double precision :: XFA

XFAN = XFA(YAN*ABC)

return
end function XFAN

!---------------------------------------------------------------------
double precision function AFR(YR)
! a=f(rho) "a" at a given radius "rho"

use const_inc, only: NA1
use status_inc, only: RHO, AMETR
use numerical_tools, only: QUADIN

implicit none

double precision, intent(in) :: YR

AFR = QUADIN(NA1, RHO, AMETR, YR)

return
end function AFR

!---------------------------------------------------------------------
double precision function AFX(YX)
! a=f(rho/roc) "a" at a given radius "x"=rho/roc

use const_inc, only: ROC

implicit none

double precision, intent(in) :: YX
double precision :: AFR

AFX = AFR(YX*ROC)

return
end function AFX

!---------------------------------------------------------------------
double precision function FRMAX(YA)
! Usage: CF1=FRMAX(TE); qmax_1/FRMAX(MU);

use const_inc, only: NA1

implicit none

double precision, intent(in) :: YA(*)

FRMAX = MAXVAL(YA(1:NA1))

return
end function FRMAX

!---------------------------------------------------------------------
double precision function FRMIN(YA)
! Usage: CF1=FRMIN(TE); qmin_1/FRMIN(MU);

use const_inc, only: NA1

implicit none

double precision, intent(in) :: YA(*)

FRMIN = MINVAL(YA(1:NA1))

return
end function FRMIN

!---------------------------------------------------------------------
double precision function RFMAX(YA)
! RFMAX [m]: Position of a minimal value of array on 0<=rho/ROC<=1
!   (Pereverzev 17-JUL-97)
! Usage: CAR1=NUES; ...=RFMAX(CAR1);
!  rmnH_RFMAX(MU);  returns position on "rho" in [m]
!  rmnH_RFMAX(MU)/ROC; returns position on "rho" normalized
!  rmnH_AMETR(RFMAX(MU)); returns position on "a" in [m]
!
!  Two expressions: MU(RFMAX(MU))
!    and FRMAX(MU)
!    return the same result.
!
! Using in FORTRAN: RFMAX(TE) or RFMAX(TE(1))

use const_inc, only: NA1
use status_inc, only: RHO

implicit none

double precision, intent(in) :: YA(*)
integer :: jmax

jmax = MAXLOC(YA(1:NA1), 1)
RFMAX = RHO(jmax)

return
end function RFMAX

!---------------------------------------------------------------------
double precision function RFMIN(YA)
! RFMIN [m]: Position of a minimal value of array on 0<=rho/ROC<=1
!   (Pereverzev 17-JUL-97)
! Usage: CF1=RFMIN(CAR1); rmnH_RFMIN(HE);

use const_inc, only: NA1
use status_inc, only: RHO

implicit none

double precision, intent(in) :: YA(*)
integer :: jmin

jmin = MINLOC(YA(1:NA1), 1)
RFMIN = RHO(jmin)

return
end function RFMIN

!---------------------------------------------------------------------
double precision function RFVAL(YA, YVAL)
! RFVAL [m]: Outermost radial (0<=rho<=ROC) position where
!   array YA(*) takes a given value (YVAL) 
!   (Pereverzev 17-MAR-98)
! Usage: (1) CAR1=NUES;        ! Coordinate of a point where NUES=0.1
!      ...=RFVAL(CAR1,.1);
!      (2) rq2_RFVAL(MU,.5)/ROC; ! ------------------------ MU=0.5
!      (3) RFVAL(AMETR,a0)  ! Recalculates "a0" in "rho"

use const_inc, only: NA1
use status_inc, only: RHO

implicit none

double precision, intent(in) :: YVAL, YA(*)
integer :: j
double precision :: YA1, YA2

RFVAL = 0.
YA1 = YA(1) - YVAL
do j=2, NA1
    YA2 = YA(j) - YVAL
    if (YA1*YA2 <= 0.) then
        YA1 = abs(YA(j-1)) + abs(YA(j))
        if (YA1 /= 0.) then
            if (1000.*abs(YA(j-1) - YA(j)) < YA1) then
                RFVAL = RHO(j)
            else
                RFVAL = (RHO(j-1)*YA(j) - RHO(j)*YA(j-1) + &
                    YVAL*(RHO(j) - RHO(j-1)))/(YA(j) - YA(j-1))
            endif
        endif
    endif
    YA1 = YA2
enddo

return
end function RFVAL

!---------------------------------------------------------------------
double precision function AFVAL(YA, YVAL)
! AFVAL [m]: Outermost radial (0<=a<=ABC) position where
!   array YA(*) takes a given value (YVAL) 
!   (Pereverzev 17-MAR-98)
! Usage: CAR1=NUES; ...=AFVAL(CAR1,.1);
!  aq2_AFVAL(MU,.5)/ROC;

use const_inc, only: NA1
use status_inc, only: AMETR

implicit none

double precision, intent(in) :: YVAL, YA(*)
integer :: j
double precision YA1, YA2

AFVAL = 0.
YA1 = YA(1) - YVAL
do j=2, NA1
    YA2 = YA(j) - YVAL
    if (YA1*YA2 <= 0.) then
        YA1 = abs(YA(j-1)) + abs(YA(j))
        if (YA1 /= 0.) then
            if (1000.*abs(YA(j-1) - YA(j)) < YA1) then
                AFVAL = AMETR(j)
            else
                AFVAL = (AMETR(j-1)*YA(j) - AMETR(j)*YA(j-1) + &
                    YVAL*(AMETR(j) - AMETR(j-1)))/(YA(j) - YA(j-1))
            endif
        endif
    endif
    YA1 = YA2
enddo

return
end function AFVAL

!---------------------------------------------------------------------
double precision function RFVEX(YA, YVAL)
! RFVEX [m]: Outermost radial (0<=rho<=ROC) position where
!   array YA(*) takes a given value (YVAL) 
!   (Pereverzev 17-MAR-98)
! Usage: CAR1=NUES; ...=RFVEX(CAR1,.1);
!  rq2_RFVEX(MU,.5)/ROC;

implicit none

double precision, intent(in) :: YVAL, YA(*)
double precision :: RFVAL

RFVEX = RFVAL(YA, YVAL)

return
end function RFVEX

!---------------------------------------------------------------------
double precision function AFVEX(YA,YVAL)
! AFVEX [m]: Outermost radial (0<=rho<=ROC) position where
!   array YA(*) takes a given value (YVAL) 
!   (Pereverzev 17-MAR-98)
! Usage: CAR1=NUES; ...=AFVEX(CAR1,.1);
!  rq2_AFVEX(MU,.5)/ROC;

implicit none

double precision, intent(in) :: YVAL, YA(*)
double precision :: AFVAL

AFVEX = AFVAL(YA, YVAL)

return
end function AFVEX

!---------------------------------------------------------------------
double precision function RFVIN(YA, YVAL)
! RFVIN [m]: Innermost radial (0<=rho<=ROC) position where
!   array YA(*) takes a given value (YVAL) 
!   (Pereverzev 17-MAR-98)
! Usage: CAR1=NUES; ...=RFVIN(CAR1,.1);
!  rq2_RFVIN(MU,.5)/ROC;

use const_inc, only: NA1
use status_inc, only: RHO

implicit none

double precision, intent(in) :: YVAL, YA(*)
integer j
double precision YA1, YA2

RFVIN = 0.
YA1 = YA(NA1) - YVAL
do j=NA1, 2, -1
    YA2 = YA(j-1) - YVAL
    if (YA1*YA2 <= 0.) then
        YA1 = abs(YA(j-1)) + abs(YA(j))
        if (YA1 /= 0.) then
            if (1000.*abs(YA(j-1)-YA(j)) < YA1) then
                RFVIN = RHO(j)
            else
                RFVIN = (RHO(j-1)*YA(j) - RHO(j)*YA(j-1) + &
                    YVAL*(RHO(j) - RHO(j-1)))/(YA(j) - YA(j-1))
            endif
        endif
    endif
enddo
YA1 = YA2

return
end function RFVIN

!---------------------------------------------------------------------
double precision function AFVIN(YA, YVAL)
! AFVIN [m]: Innermost radial (0<=a<=ABC) position where
!   array YA(*) takes a given value (YVAL) 
!   (Pereverzev 17-MAR-98)
! Usage: CAR1=NUES; ...=AFVIN(CAR1,.1);
!  rq2_AFVIN(MU,.5)/ROC;

use const_inc, only: NA1
use status_inc, only: AMETR

implicit none

double precision, intent(in) :: YVAL, YA(*)

integer :: j
double precision :: YA1, YA2

AFVIN = 0.
YA1 = YA(NA1) - YVAL
do j=NA1, 2, -1
    YA2 = YA(j-1) - YVAL
    if (YA1*YA2 <= 0.) then
        YA1 = abs(YA(j-1)) + abs(YA(j))
        if (YA1 /= 0.) then
            if (1000.*abs(YA(j-1) - YA(j)) < YA1) then
                AFVIN = AMETR(j)
            else
                AFVIN = (AMETR(j-1)*YA(j) - AMETR(j)*YA(j-1) + &
                    YVAL*(AMETR(j) - AMETR(j-1)))/(YA(j) - YA(j-1))
            endif
        endif
    endif
enddo
YA1 = YA2

return
end function AFVIN

!---------------------------------------------------------------------
double precision function RECR(YZ, N)
! The function translates ECR frequency into a radius of resonance surface
!   for N-th ECR harmonic at z-plane z=YZ
! Input:
! YZ    [m] vertical pozition with respect to the chamber centre
! FECR  [GHz] ECR frequency 
! N     [ ] ECR harmonic number
! BTOR  [T] toridal field at the chamber centre
! RTOR  [m] major radius of the chamber centre
! Output:
! RECR  [m] rho_res
!    Optional output: a_res, shift, elongation, triangularity
! Example: CF3=RECR(0.,2)/ROC; PEECR=QECR*GAUSS(CF3,0.1);

use const_inc, only: NA1, RTOR, BTOR, FECR
use status_inc, only: AMETR, RHO
use numerical_tools, only: QUADIN

implicit none

integer, intent(in) :: N
double precision, intent(in) :: YZ

double precision :: YA, YR, RZ2A

FECR = 140.
YR = 28.*N*BTOR*RTOR/FECR
YA = RZ2A(YR, YZ, NA1)
RECR = QUADIN(NA1, AMETR, RHO, YA)

return
end function RECR

!---------------------------------------------------------------------
double precision function GAUSS(YX, YW, j)
! Gauss distribution function 
! Note: Both arguments are dimensionless and the distribution are
!       given with respect to the 
!       normalized effective minor radius "rho".
!       The 3rd argument is added by the ASTRA compiler
! Examples:
! PE=GAUSS(.1,0.2); Gaus\GAUSS(0.,.1);
!      (Pereverzev 01-AUG-96)
!    Does not work for multiple calls, due to save?

use const_inc, only: NA1, ROC, HRO
use status_inc, only: RHO, VR

implicit none

integer, intent(in) :: j
double precision, intent(in) :: YX, YW

integer :: jj
double precision :: YPOW, YR

save YPOW
data YPOW/0./

if (j == 1) then
    YPOW = 0.
    do jj = 1, NA1-1
        YR = (RHO(jj)/ROC - YX)/YW 
        YPOW = YPOW + exp(-YR*YR)*VR(jj)
    enddo
    YPOW = YPOW*HRO
endif

if (YPOW <= 0.) then
! Can happen when the function is called from "tmp/detvar.inc" before METRIC
    GAUSS = 0.
else
    GAUSS = exp(-((RHO(j)/ROC - YX)/YW)**2)/YPOW
endif

return
end function GAUSS

!---------------------------------------------------------------------
double precision function ASTEP(YA, j)
! Step function with dimensional [m] radial argument "a"
! Returns 0. for a < YA
!     and 1. otherwise
! Examples:
! PE=ASTEP(.1)-ASTEP(0.2); Step\ASTEP(.1);
! CV1=XQMINB; HE=HE*ASTEP(CV1);
!      (Pereverzev 01-AUG-96)

use status_inc, only: AMETR

implicit none

integer, intent(in) :: j
double precision, intent(in) :: YA

if (AMETR(j) < YA) then
    ASTEP = 0.
else
    ASTEP = 1.
endif

return
end function ASTEP

!---------------------------------------------------------------------
double precision function RSTEP(YR, j)
! RSTEP: Step function with dimensional [m] radial argument "rho"
! Returns 0. for rho < YR
!     and 1. otherwise
! Examples:
! PE=RSTEP(.1)-RSTEP(0.2); Step\RSTEP(.1);
!      (Pereverzev 01-AUG-96)

use status_inc, only: RHO

implicit none

integer, intent(in) :: j
double precision, intent(in) :: YR

if (RHO(j) < YR) then
    RSTEP = 0.
else
    RSTEP = 1.
endif

return
end function RSTEP

!---------------------------------------------------------------------
double precision function XSTEP(YX, j)
! XSTEP: Step function with a dimensionless radial argument
! Returns 0. for rho < YX*ROC
!     and 1. otherwise
! Examples:
! PE=XSTEP(.1)-XSTEP(0.2); Step\XSTEP(.1);
! CV1=XQMINB; HE=HE*XSTEP(CV1);
!      (Pereverzev 01-AUG-96)

use const_inc, only: ROC
use status_inc, only: RHO

implicit none

integer, intent(in) :: j
double precision, intent(in) :: YX

if (RHO(j)/ROC < YX) then
    XSTEP = 0.
else
    XSTEP = 1.
endif

return
end function XSTEP

!---------------------------------------------------------------------
double precision function STEP(YX)
! STEP: Step function with a dimensionless radial argument
! Returns 0. for YX < 0
!     and 1. otherwise
! Examples:
! PE=STEP(CAR1); Step\STEP(CAR1);
!      (Pereverzev 17-FEB-99)

implicit none

double precision, intent(in) :: YX

if (YX < 0.) then
    STEP = 0.
else
    STEP = 1.
endif

return
end function STEP

!---------------------------------------------------------------------
double precision function GRAD(Y, j)
! Gradient
! Only a radially dependent array may be the 1st parameter of the function
! Examples:
!    out\GRAD(CAR3) !Radial profile of gradient CAR3 dCAR3/dRo
!    out_GRAD(CAR3B) !Gradient CAR3 dCAR3/dRo at the boundary
!    out_GRAD(CAR3C) !Gradient CAR3 dCAR3/dRo at the center

use const_inc, only: HRO, NA, NA1

implicit none

integer, intent(in) :: j
double precision, intent(in) :: Y(*)

if (j < NA) then
    GRAD = (Y(j+1) - Y(j))/HRO
else
    GRAD = (Y(NA1) - Y(NA))/HRO
endif

return
end function GRAD

!---------------------------------------------------------------------
double precision function GRADS(Y, j)
! Gradient
! Only a radially dependent array may be the 1st parameter of the function
! Examples:
!    out\GRADS(CAR3) !Radial profile of gradient CAR3 dCAR3/dRo
!    out_GRADS(CAR3B) !Gradient CAR3 dCAR3/dRo at the boundary
!    out_GRADS(CAR3C) !Gradient CAR3 dCAR3/dRo at the center

use const_inc, only: HRO, ROC, NA, NA1

implicit none

integer, intent(in) :: j
double precision, intent(in) :: Y(*)

if (j < NA) then
    GRADS = (Y(j+1) - Y(j))/HRO
else
    GRADS = (Y(NA1) - Y(NA))/(ROC - NA*HRO)
endif

return
end function GRADS

!---------------------------------------------------------------------
double precision function LININT(ARR, YR)
! Line integral {0,R} of any array
! Only a radially dependent array may be the 1st parameter of the function
! Examples:
!    out\Linint(CAR3)     !Radial profile of CAR3 line integral
!    out_Linint(CAR3,Ro); !Line integral {0,Ro} of CAR3
!    out_Linint(CAR3B)    !Total line integral of CAR3 (0,ROC)
!   (Yushmanov 26-DEC-90)

use const_inc, only: HRO
use status_inc, only: AMETR

implicit none

double precision, intent(in) :: YR, ARR(*)
integer :: J, JK
double precision :: YDR

call yrjkdr(YR, JK, YDR)
LININT = 0.d0
do J=2, JK
    LININT = LININT + (AMETR(j) - AMETR(j-1))*(ARR(j) + ARR(j-1))
enddo
LININT = LININT - (arr(JK) + arr(JK-1))*(AMETR(JK) - AMETR(JK-1))*(float(JK) - YR/HRO)

! Line average: divide by 2*ABC

return
end function LININT

!---------------------------------------------------------------------
double precision function VINT(ARR, YR)
! Volume integral {0, R} of any array
! Only a radially dependent array may be the 1st parameter of the function
! Examples:
!    out\Vint(CAR3) !Radial profile of CAR3 volume integral
!    out_Vint(CAR3,Ro); !Volume integral {0,Ro} of CAR3
!    out_Vint(CAR3B)    !Total volume integral of CAR3 (0,ROC)
!   (Yushmanov 26-DEC-90)

use const_inc, only: HRO
use status_inc, only: VR

implicit none

double precision, intent(in) :: YR, ARR(*)
integer :: J, JK
double precision :: YDR

call yrjkdr(YR, JK, YDR)
VINT = 0.
do J=1, JK
    VINT = VINT + ARR(J)*VR(J)
enddo
VINT = HRO*(VINT - ARR(JK)*YDR)

return
end function VINT

!---------------------------------------------------------------------
double precision function VINTO(ARR, YR)
! Volume integral {0,R} of any array
! Only a radially dependent array may be the 1st parameter of the function
! Examples:
!    out\Vint(CAR3) !Radial profile of CAR3 volume integral
!    out_Vint(CAR3,Ro); !Volume integral {0,Ro} of CAR3
!    out_Vint(CAR3B)    !Total volume integral of CAR3 (0,ROC)
!   (Yushmanov 26-DEC-90)

use const_inc, only: HRO
use status_inc, only: VRO

implicit none

double precision, intent(in) :: YR, ARR(*)
integer J, JK
double precision :: YDR

call yrjkdr(YR, JK, YDR)
VINTO = 0.
do J=1, JK
    VINTO = VINTO + ARR(J)*VRO(J)
enddo
VINTO = HRO*(VINTO - ARR(JK)*YDR)

return
end function VINTO

!---------------------------------------------------------------------
double precision function IINT(ARR, YR)
! Integral {0,R} of current density
!  Iint=integral {0,R} (ARR/IPOL**2)dV*IPOL/(GP2*Ro)
!  Only radial dependent array may be a parameter of the function
!  Examples:
!    out\Iint(CU) !Radial profile of toroidal current
!    out_Iint(CD,Ro)    !Toroidal driven current inside {0,Ro}
!    out_Iint(CUB)      !Total toroidal current =Iint(CU,ROC); (=IPL)
!   (Pereverzev 23-OCT-99)

use const_inc, only: HRO, NA, GP2
use status_inc, only: RHO, IPOL, G33

implicit none

double precision, intent(in) :: YR, ARR(*)

integer :: J, JK
double precision :: YA, YDR

IINT = 0.
if (YR <= 0.) return

JK = YR/HRO + 1. - 1.E-4
if (JK > NA) JK = NA

YA = 0.
do J=1, JK
    IINT = IINT + YA
    YA   = ARR(J)*RHO(J)/(G33(J)*IPOL(J)**3)
    YDR = YR - JK*HRO + HRO
enddo
if (JK >= NA) then
    YDR = min(YDR, HRO)
endif
IINT = GP2*IPOL(JK)*(HRO*IINT + YDR*YA)

return
end function IINT

!---------------------------------------------------------------------
integer function NODE(YR)
! Radial node number nearest to YR (radius larger than YR)

use const_inc, only: ROC, HRO, NA1

implicit none

double precision, intent(in) :: YR

if (YR <= 0.) then
    NODE = 1
else if (YR > ROC) then
    NODE = NA1
else
    NODE = int(YR/HRO) + 1
endif
NODE = min(NA1, NODE)

return
end function NODE

!---------------------------------------------------------------------
double precision function RADIAL(ARR, YR)
! Interpolation of ARR to the radial position R
!   It is assumed that the array ARR(j) is given on the grid RHO(j)

use const_inc, only: NA1
use status_inc, only: RHO
use numerical_tools, only: qinterp

implicit none

double precision, intent(in) :: YR, ARR(*)
double precision, dimension(1) :: rad_in, rad_out

rad_in(1) = YR
call qinterp(RHO(1:NA1), ARR(1:NA1), NA1, rad_in, rad_out, 1)
RADIAL = rad_out(1)

return
end function RADIAL

!---------------------------------------------------------------------
double precision function RADINT(ARR, YR)
! Interpolation of ARR to the radial position R
!   It is assumed that the array ARR(j) is given at points j*HRO

use const_inc, only: NA1, HRO

implicit none

double precision, intent(in) :: YR, ARR(*)
integer :: JK
double precision :: YDR

JK = YR/HRO
if (JK >= NA1) then
    RADINT = ARR(NA1)
else if (JK == 0) then
    RADINT = (4.*ARR(1) - ARR(2))/3.*(1. - YR*YR)
else if (JK < 0) then
    RADINT = (4.*ARR(1) - ARR(2))/3.
else
    YDR = YR - JK*HRO
    RADINT = (YDR*ARR(JK+1) + (HRO - YDR)*ARR(JK))/HRO
endif

return
end function RADINT

!---------------------------------------------------------------------
double precision function ATR(ARR, YR)
! ATR synonym for RADIAL

implicit none

double precision, intent(in) :: YR, ARR(*)
double precision :: RADIAL

ATR = RADIAL(ARR, YR)

return
end function ATR

!---------------------------------------------------------------------
double precision function ATX(ARR, YR)
! Interpolation of ARR to the radial position R
!   It is assumed that the array ARR(j) is given on the grid RHO(j)

use const_inc, only: ROC

implicit none

double precision, intent(in) :: YR, ARR(*)
double precision :: RADIAL

ATX = RADIAL(ARR, YR*ROC)

return
end function ATX

!---------------------------------------------------------------------
double precision function CUT(X, Y)
! CUT The function cuts off values of Y out of the range:
!    -X < Y < X
! Usage in a model: Name\cut(.005,CAR1)\.005;
!  or  nues\cut(100.,HEXP)\100

implicit none

double precision, intent(in) :: X, Y

CUT = max(-X, min(X, Y))

return
end function CUT

!---------------------------------------------------------------------
double precision function FRAMP(T1, T2)
! FRAMP: linear ramp from 0 at Time<T1 to 1 at Time>T2
! Example:
! IPL=.1+.2*framp(0.21,0.27)
!   (Yushmanov 26-DEC-90)

use const_inc, only: TIME

implicit none

double precision, intent(in) :: T1, T2

if (T1 >= T2) then
    write(*, *) ' Function FRAMP(t1,t2) error: t1 > t2'
    FRAMP = 0.
    return
endif
if (TIME <= T1) then
    FRAMP = 0.
else if (TIME >= T2) then
    FRAMP = 1.
else
    FRAMP = (TIME - T1)/(T2 - T1)
endif

return
end function FRAMP

!---------------------------------------------------------------------
double precision function FJUMP(T1)
! FJUMP: jump from 0 to 1 at Time=T1
! Example:
! IPL=.1+.2*fjump(.21d0)
!   (Yushmanov 26-DEC-90)

use const_inc, only: TIME

implicit none

double precision, intent(in) :: T1

if (TIME <= T1) then
    FJUMP = 0.
else
    FJUMP = 1.
endif

return
end function FJUMP

!---------------------------------------------------------------------
double precision function FTBOX(T1, T2)
! Box-like time dependence:
! FTBOX = | 0   if  time <= T1  or time >= T2
!         | 1   otherwise
! Example:
! CV3 = 0.2*ftbox(0.21d0,5.d-1)
!
!       CAR3=FBOX(
!   (Pereverzev 26-JAN-06)

use const_inc, only: TIME

implicit none

double precision, intent(in) :: T1, T2

if (TIME <= T1 .or. TIME >= T2) then
    FTBOX = 0.
else
    FTBOX = 1.
endif

return
end function FTBOX

!---------------------------------------------------------------------
double precision function FXBOX(X1, X2, j)
! FXBOX: Box-like radial (argument x=rho_norm)  dependence:
!       FXBOX = | 0   if  x <= X1  or x >= X2
!  | 1   otherwise
! Example:
!        CAR3=FXBOX(5.d-1,7.d-1)
!     (Pereverzev 26-JAN-06)

use const_inc, only: ROC
use status_inc, only: RHO

implicit none

integer, intent(in) :: j
double precision, intent(in) :: X1,X2

if (RHO(j) <= X1*ROC .or. RHO(j) >= X2*ROC) then
    FXBOX = 0.
else
    FXBOX = 1.
endif

return
end function FXBOX

!---------------------------------------------------------------------
double precision function TIMDER(Y)
! Time derivative
!  Examples:
!     cv1=TimDer(IPL);  cv2=TIMDER(WTOTB);
!     dIdt_cv1;   dWdt_cv2;
! Note! the function cannot be used in the time output directly !!!!!!!!
!   (Yushmanov 13-FEB-91)
!hanged by Pereverzev 15.10.98

use const_inc, only: TIME

implicit none

integer, parameter :: NLOC=2200
double precision, intent(in) :: Y

integer :: IY, ITSME, ICALL
double precision, dimension(NLOC) :: YO, YT

save ITSME, ICALL, YO, YT
data ITSME/0/ ICALL/0/ YO/NLOC*0./ IY/0/

call getid(Y, itsme, IY)

TIMDER = 0.
if (IY < 0 .or. IY > NLOC) then
    write(*, *) "            Calling from TIMDER"
    write(*, *) ' too many time derivatives >', NLOC
    return
endif

if (ICALL == 0.) then
    ICALL  = 1
    TIMDER = 0.
    YO(IY) = Y
    YT(IY) = TIME
else if (TIME > YT(IY)) then
    TIMDER = (Y - YO(IY))/(TIME - YT(IY))
    YO(IY) = Y
    YT(IY) = TIME
endif

return
end function TIMDER

!---------------------------------------------------------------------
double precision function TIMINT(Y)
! Time integral
!  Examples:
!  cv3=TimInt(QITOTB)
!      A_cv3 ! Total ion energy input
!
! Note: the function cannot be used in the time output directly
!   (Yushmanov 13-FEB-91)
! Changed by Pereverzev 15.10.98

use const_inc, only: TIME

implicit none

integer, parameter :: NLOC=2200
double precision, intent(in) :: Y

integer :: IY, ITSME, ICALL
double precision, dimension(NLOC) :: YO, YT

save ITSME, YO, YT, ICALL
data ITSME/0/ ICALL/0/ YO/NLOC*0./ IY/0/

call getid(Y, itsme, IY)

TIMINT = 0.
if (IY < 0 .or. IY > NLOC) then
    write(*, *) "            Calling from TIMINT"
    write(*, *) ' too many time integrals >', NLOC
    return
endif

if (ICALL == 0.) then
    ICALL  = 1
    TIMINT = 0.
    YT(IY) = TIME
else
    TIMINT = YO(IY) + Y*(TIME - YT(IY))
    YO(IY) = TIMINT
    YT(IY) = TIME
endif

return
end function TIMINT

!---------------------------------------------------------------------
double precision function TIMAVG(Y, YTINT)
! Timavg: time average
!  Examples:
!  cv3=Timavg(NNCL,1.)
!      N_cv3 ! Neutral gas density averaged over 1 sec
!     Pereverzev 01.06.02
!---------------------------------------------------------------------
!   Y      - quantity to be averaged over YTINT sec
!   YTINT  - sampling time 
!---------------------------------------------------------------------
!   IY     - ID of the input quantity Y
!   YI(IY) - integral of Y
!   YT(IY) - time of the previous calling

use const_inc, only: TIME

implicit none

integer, parameter :: NLOC=2200
double precision, intent(in) :: Y, YTINT

integer :: IY, ITSME, j
double precision :: YST
double precision, dimension(NLOC) :: YI, YT

save ITSME, YI, YT
data ITSME/0/ YI/NLOC*-1.E9/  IY/0/

call getid(Y, itsme, IY)

TIMAVG = 0.
if (IY < 0 .or. IY > NLOC) then
    write(*, *) "            Calling from TIMAVG"
    write(*, *) ' too many time averages >', NLOC
    return
endif

if (YI(IY) < -0.9E9) then  ! 1st call for IY
    YI(IY) = 0.
    YST = TIME
    TIMAVG = Y
    YT(IY) = TIME
else
    j = TIME/YTINT
    YST = j*YTINT
    if (YST < YT(IY)) then  ! intermediate call
        YI(IY) = YI(IY) + Y*(TIME - YT(IY))
        TIMAVG = YI(IY)/(TIME - YST)
    else     ! YTINT starting call
        YI(IY) = YI(IY) + Y*(YST - YT(IY))
        TIMAVG = YI(IY)/YTINT
        YI(IY) = Y*(TIME - YST)
    endif
    YT(IY) = TIME
endif

return
end function TIMAVG

!---------------------------------------------------------------------
double precision function FIXVAL(Y, YTIME)
! Returns the value Yvar(time<=Ytime)
!  Example:
!     dli_LINTB-FIXVAL(LINTB,1.)
! This output produces an increment of l_i(t) with respect to t=1sec
!   (Pereverzev 15-OCT-98)

use const_inc, only: TIME

implicit none

integer, parameter :: NLOC=2200
double precision, intent(in) :: Y, YTIME

integer IY, ITSME
double precision :: YO(NLOC)

save ITSME, YO
data ITSME/0/ IY/0/
! IY is the ID (ordinal number) of "Y"

call getid(Y, itsme, IY)

FIXVAL = 0.
if (IY < 0 .or. IY > NLOC) then
    write(*, *) ' >>> FIXVAL >>> too many calls: > ', NLOC
    return
endif

if (TIME <= YTIME) YO(IY) = Y
FIXVAL = YO(IY)

return
end function FIXVAL

!---------------------------------------------------------------------
double precision function FTAV(Y, YTAV)
! Returns g_i=g_{i-1}exp(-tau/YTAV)+f_i*[1-exp(-tau/YTAV)]
!  i.e.
!  g_i=f_i  if  tau << YTAV
!  g_i=g_{i-1} if  tau >> YTAV
!---------------------------------------------------------------------
! Floating Time AVerage, exponential decay with YTAV
! 
!  Examples:
!  CV2=FTAV(cv1,CF1)
!  CAR2=FTAV(UPL,3*tau);
!  CAR3=HARL
!  HE=FTAV(CAR3,.01);
! Note:
!      Do not use ...=FTAV(UPLB,.1)    (UPLB - ASTRA abbreviation)
!      Use CV1=UPLB; CV2=FTAV(CV1,.1)
!
!      Do not use ...=FTAV(HARL,.1) (HARL - ASTRA formula)
!      Use CAR3=HARL; HE=...+FTAV(CAR3,.01)
!
! G.W. Pacher (18/01/1994)
! Changed by Pereverzev 15.10.98

use const_inc, only: TAU

implicit none

integer, parameter :: NLOC=2200
double precision, intent(in) :: Y, YTAV

integer :: IY, ITSME, ICALL
double precision :: YO(NLOC)

save ITSME, ICALL, YO
data ITSME/0/ ICALL/0/ YO/NLOC*0./ IY/0/

call getid(Y, itsme, IY)
! IY is the ID (ordinal number) of "Y"
! itsme is the ID (ordinal number) of the calling function (FTAV)

FTAV = 0.
if (IY < 0 .or. IY > NLOC) then
    write(*, *) ' >>> FTAV >>> vuffer overflow: >', NLOC
    return
endif

if (ICALL == 0) then
    FTAV = Y
    YO(IY) = FTAV
    ICALL  = 1
else
    if (YTAV <= .1*TAU) then
        FTAV = Y    ! Return the input value
    else
        FTAV = Y + (YO(IY) - Y)*EXP(-TAU/YTAV) ! Return a weighted value
    endif
    YO(IY) = FTAV    ! Save the previous value
endif

return
end function FTAV

!---------------------------------------------------------------------
double precision function FTAV3(Y)
! Returns the average of the quantity Y over JBASE time steps
! The argument to FTAV3 can be either an ASTRA variable or ASTRA array
! 
!  Examples:
!  CV2=FTAV3(cv1)
! Note:
!      Do not use ...=FTAV3(UPLB)      (UPLB - ASTRA abbreviation)
!      Use CV1=UPLB; CV2=FTAV3(CV1)
!reated by Pereverzev 8.08.2007

use const_inc, only: TIME

implicit none

integer, parameter :: JBASE=1000
double precision, intent(in) :: Y

integer j, JY
double precision YY(JBASE), YTIME, YAV

save YY, YTIME, YAV, JY
data JY/0/ YTIME/-1.d37/ YAV/0.d0/

if (JY /= 0)  then
    FTAV3 = YAV
    if (YTIME == TIME) return
endif

if (jy < JBASE) then
    JY = JY + 1
    YY(jy) = Y
    YAV = 0.
    do j=1, JY
        YAV = YAV + YY(j)
    enddo
    YAV = YAV/JY
else
    YAV = YAV + (Y - YY(1))/JY
    do j=2, JY
        YY(j-1) = YY(j)
    enddo
    YY(jy) = Y
endif
FTAV3 = YAV
YTIME = TIME

return
end function FTAV3

!---------------------------------------------------------------------
double precision function FTAV2(Y)
! Returns an average of the quantity Y over JBASE time steps
! The argument to FTAV2 can be either an ASTRA variable or ASTRA array
! 
!  Examples:
!  CV2=FTAV2(cv1)
! Note:
!      Using <dt>_FTAV2(TAU) is allowed although can give a wrong
!  result: It averages over time output times instead of
!  all calculation times.
!      Do not use ...=FTAV2(UPLB)      (UPLB - ASTRA abbreviation)
!      Use CV1=UPLB; CV2=FTAV2(CV1)
!reated by Pereverzev 8.08.2007

use const_inc, only: TIME

implicit none

integer, parameter :: JBASE=100
double precision, intent(in) :: Y

integer j, JY
double precision :: YY(JBASE), YTIME, YAV

save YY, YTIME, YAV, JY
data JY/0/ YTIME/-1.d37/ YAV/0.d0/

if (JY /= 0) then
    FTAV2 = YAV
    if (YTIME == TIME) return
endif

if (jy < JBASE) then
    JY = JY + 1
    YY(jy) = Y
    YAV = 0.
    do j=1, JY
        YAV = YAV + YY(j)
    enddo
    YAV = YAV/JY
else
    YAV = YAV + (Y - YY(1))/JY
    do j=2, JY
        YY(j-1) = YY(j)
    enddo
    YY(jy) = Y
endif
FTAV2 = YAV
YTIME = TIME

return
end function FTAV2

!---------------------------------------------------------------------
double precision function FTMIN(Y)
! The function returns min_t[Y]
! The argument to FTMIN can be either an ASTRA variable or ASTRA array
! 
!  Examples:
!  CV2=FTMIN(cv1)
!  CAR2=FTMIN(UPL);
!  CAR3=HARL
!  HE=FTMIN(CAR3);
! Note:
!      Do not use ...=FTMIN(UPLB)      (UPLB - ASTRA abbreviation)
!      Use CV1=UPLB; CV2=FTMIN(CV1)
!
!      Do not use ...=FTMIN(HARL)     (HARL - ASTRA formula)
!      Use CAR3=HARL; HE=...+FTMIN(CAR3)
!
! Changed by Pereverzev 15.10.98

implicit none

integer, parameter :: NLOC=2200
double precision, intent(in) :: Y

integer :: IY, ITSME
double precision :: YO(NLOC)

save ITSME,YO
data ITSME/0/ YO/NLOC*1.E37/  IY/0/

! IY is the ID (ordinal number) of "Y"
call getid(Y, itsme, IY)

if (IY < 0 .or. IY > NLOC) then
    FTMIN = 0.
    write(*, *)' >>> FTMIN >>> too many calls: >', NLOC
else
    FTMIN = 1.E37
    if (YO(IY) <= Y) then
        FTMIN = YO(IY)
    else
        FTMIN = Y
        YO(IY) = Y
    endif
endif

return
end function FTMIN

!---------------------------------------------------------------------
double precision function FTMAX(Y)
! The function returns max_t[Y]
! The argument to FTMAX can be either an ASTRA variable or ASTRA array
! 
!  Examples:
!  CV2=FTMAX(cv1)
!  CAR2=FTMAX(UPL);
!  CAR3=HARL
!  HE=FTMAX(CAR3);
! Note:
!      Do not use ...=FTMAX(UPLB)      (UPLB - ASTRA abbreviation)
!      Use CV1=UPLB; CV2=FTMAX(CV1)
!
!      Do not use ...=FTMAX(HARL)     (HARL - ASTRA formula)
!      Use CAR3=HARL; HE=...+FTMAX(CAR3)
!
! G.V.Pereverzev 15.10.98

implicit none

integer, parameter :: NLOC=2200
double precision, intent(in) :: Y

integer :: IY, ITSME
double precision :: YO(NLOC)

save ITSME, YO
data ITSME/0/ YO/NLOC*-1.E37/  IY/0/
! IY is the ID (ordinal number) of "Y"

call getid(Y, itsme, IY)

if (IY < 0 .or. IY > NLOC) then
    FTMAX = 0.
    write(*, *) ' >>> FTMAX >>> too many calls: >', NLOC
else
    FTMAX = -1.E37
    if (YO(IY) >= Y) then
        FTMAX = YO(IY)
    else
        FTMAX = Y
        YO(IY) = Y
    endif
endif

return
end function FTMAX

!---------------------------------------------------------------------
double precision function ADTRMC(x, r, y, na1)
! Adiabatic compression term

implicit none

integer, intent(in) :: na1
double precision, intent(in) :: x(na1), r
double precision, intent(out) :: y(na1)

integer :: i
double precision :: dum1, dum2, z1

do i=2, na1-1
    dum1 = (y(i+1) - y(i))/(x(i+1) - x(i))
    dum2 = (y(i) - y(i-1))/(x(i) - x(i-1))
    z1 = (dum1 + dum2)/2.0
    y(i) = x(i)*z1*r
enddo
i = 1
dum1 = (y(i+1) - y(i))/(x(i+1) - x(i))
dum2 = (y(i) - y(i+1))/(2*x(i))
z1 = (dum1 + dum2)/2.0
y(i) = x(i)*z1*r

return
end function ADTRMC

!---------------------------------------------------------------------
double precision function RZ2A(R_in, Z_in, nx_in)
! Returns A(r, z) where
!
! r = R0 + del(A) + A*[cos(theta)-tri(A)*sin^2(theta)]
! z = UPD + A*elo(A)*sin(theta)
!
! and del(j), elo(j), tri(j) are given as arrays[1:N] 
!      on the grid A=AMETR(j)

use const_inc, only: RTOR, AB
use status_inc, only: AMETR, SHIF, SHIV, ELON, TRIA
use numerical_tools, only: QUADIN

implicit none

integer, intent(in) :: nx_in
double precision, intent(in) :: R_in, Z_in

integer :: j
double precision :: Y1, YAS, YAO, YA
double precision :: YHOR, YVER, YELO, YTRI

YAS = ((Z_in - SHIV(nx_in))/ELON(nx_in))**2
YA = sqrt(YAS + (R_in - RTOR - SHIF(nx_in))**2)

do j=1, 20
    YAO = YA
    YHOR = QUADIN(nx_in, AMETR(1: nx_in), SHIF(1: nx_in), YA)
    YELO = QUADIN(nx_in, AMETR(1: nx_in), ELON(1: nx_in), YA)
    YTRI = QUADIN(nx_in, AMETR(1: nx_in), TRIA(1: nx_in), YA)
    YVER = QUADIN(nx_in, AMETR(1: nx_in), SHIV(1: nx_in), YA)
    YAS = ((Z_in - YVER)/YELO)**2
    if (YA > 1.E-4) then
        Y1 = YAS/YA
    else
        Y1 = YAS
    endif
    YA = sqrt(YAS + (R_in - RTOR - YHOR + YTRI*Y1)**2)
    if (abs(YA - YAO) < 1.E-10) EXIT
enddo
RZ2A = min(YA, AB)

return
end function RZ2A
