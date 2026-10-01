! WCE [MJ]: Ion thermal enegry stored in the plasma core (WE - Pedestal)
!  0.0016 * 3/2* Integral {0, R} (Ne*Te-Neb*Teb) dV
!   (Pereverzev 26-MAR-00)
double precision function WCER(YR)

use scalars, only: NA1
use status, only: NE, TE
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR
integer :: j
double precision, dimension(NA1) :: YY

do j=1, NA1
   YY(j) = NE(J)*TE(J) - NE(NA1)*TE(NA1)
enddo
WCER = 0.0024*VINT(YY, YR)

end function WCER
