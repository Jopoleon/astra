!======================================================================|
! BETAJR []: Beta poloidal (r)         (Pereverzev 12-JAN-90)
!
! In this subroutine, the dimensionless \beta_J is defined as
!
!   2*c*c                        2*c*c
!  ------ * {Int(p*dS) - pS} = - ------ * Int[S*dp]
!    I*I                          I*I
!
! Equivalent definition:
!
!    c*c           p*dV              c*c         V*dp
!  -------- * {Int[----] - pV} = - ------- * Int[----]
!  \pi*I*I          r              \pi*I*I        r
!
!----------------------------------------------------------------------|
! The same in SI units:
!
!     8\pi                           8\pi
!  ---------*{Int(p*dS) - pS} = - --------- * Int[S*dp]
!  \mu_0*I*I                      \mu_0*I*I
!
! or
!
!      4            p*dV                4           V*dp
!  --------- * {Int[----] - pV} = - --------- * Int[----]
!  \mu_0*I*I         r              \mu_0*I*I        r
!
!----------------------------------------------------------------------|
! Note!
!   The simplified representation S ~= \pi\lambda*a^2 is employed here.
!   It is valid for 3-moment equilibrium only.
!======================================================================|
double precision function BETAJR(YR)

use const_inc, only: ROC, HRO, GP2, NA1, BTOR, RTOR, NA
use status_inc, only: VR, ELON, TE, TI, NE, NI, PBLON, PBPER, MU, &
   G22, IPOL, AMETR

implicit  none

integer J, JR
double precision YR, YB, YP, YP1
!----------------------------------------------------------------------|
JR = YR/HRO+0.5
if (JR .lt. 2) JR=2
if (JR .gt. NA) JR=NA
YB = 0.
YP = NE(1)*TE(1)+NI(1)*TI(1)
do J=1,JR
   YP1 = NE(J+1)*TE(J+1)+NI(J+1)*TI(J+1)
   YB = YB+(YP-YP1)*ELON(J)*AMETR(J)**2
   YP = YP1
enddo
BETAJR = 6.4E-4*GP2*YB*(RTOR/(G22(JR)*IPOL(JR)*BTOR*JR*HRO*MU(JR)))**2

return
end function betajr
