! WC [MJ]: Exp. thermal enegry stored in the plasma core (WTOT - Pedestal)
! 0.0016*3/2*Integral {0, R} (Ni*Ti+Ne*Te-Nib*Tib-Nexb*Teb)dV
!   (Pereverzev 26-MAR-00)
double precision function WCR(YR)

use parameter_inc, only: NRD
use const_inc, only: NA1
use status_inc, only: NI, TI, NE, TE

implicit none

double precision, intent(in) :: YR
integer j
double precision :: VINT, YY(NRD)

do j = 1, NA1
   YY(j) = NI(J)*TI(J) - NI(NA1)*TI(NA1) + NE(J)*TE(J) - NE(NA1)*TE(NA1)
enddo
WCR = 0.0024*VINT(YY, YR)

end function WCR
