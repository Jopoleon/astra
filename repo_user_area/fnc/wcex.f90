! WCEX [MJ]: Ion thermal enegry stored in the plasma core (WEX - Pedestal)
!  0.0016 * 3/2* Integral {0, R} (Nex*Tex-Nexb*Texb) dV
!   (Pereverzev 26-MAR-00)
double precision function WCEXR(YR)

use parameter_inc, only: NRD
use const_inc, only: NA1
use status_inc, only: NEX, TEX

implicit none

double precision, intent(in) :: YR
integer  J
double precision :: VINT, YY(NRD)

do j=1, NA1
   YY(j) = NEX(J)*TEX(J) - NEX(NA1)*TEX(NA1)
enddo
WCEXR = 0.0024*VINT(YY, YR)

end function WCEXR
