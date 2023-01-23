!=======================================================================
double precision   function   WBPOLR(YR1)
!-----------------------------------------------------------------------
! WBPOL [MJ]: Energy contents in the poloidal magnetic field
!        /                               /
!         1    |                          1    | d\psi    IPOL
! Wi = -------*|(Bpol**2)dV (Gauss)  = -------*|(-----)^2 ---- G22 d\rho (SI)
!       8*pi   |                        2\mu_0 | d\rho    RTOR
!        /                               /
!              V                               0
!
!   Wi = 0.5*L_i*(I_pl)^2;  SI units: J,  Hn,  A
!
!   (G.Pereverzev 14-NOV-94)
!      (corrected  6-JUL-99)
!-----------------------------------------------------------------------

use const_inc, only: HRO, GP, BTOR, RTOR
use status_inc, only: MU, IPOL, G22

implicit none

double precision, intent(in) :: YR1
integer :: J, JK, J1
double precision :: YR

YR = YR1/HRO
JK = YR
WBPOLR = (YR*YR - JK*JK)*(YR*YR + JK*JK)*MU(JK+1)**2 * IPOL(JK+1)*G22(JK+1)/(JK+1)
DO J=1, JK
   J1 = J + J - 1
   WBPOLR = WBPOLR + (J1*MU(J)**2*(2.*J*J - J1))*IPOL(J)*G22(J)/J
enddo
WBPOLR = 1.25*GP*WBPOLR*HRO*HRO*HRO*BTOR*BTOR/RTOR

return
end function WBPOLR
