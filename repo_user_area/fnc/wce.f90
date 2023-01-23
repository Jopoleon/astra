! WCE [MJ]: Ion thermal enegry stored in the plasma core (WE - Pedestal)
!  0.0016 * 3/2* Integral {0, R} (Ne*Te-Neb*Teb) dV
!   (Pereverzev 26-MAR-00)
double precision function WCER(YR)

use parameter_inc, only: NRD
use const_inc, only: NA1
use status_inc, only: NE, TE

implicit none

double precision, intent(in) :: YR
integer  J
double precision :: VINT, YY(NRD)

do j=1, NA1
   YY(j) = NE(J)*TE(J) - NE(NA1)*TE(NA1)
enddo
WCER = 0.0024*VINT(YY, YR)

end function WCER
