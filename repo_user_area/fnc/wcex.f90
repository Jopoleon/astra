! WCEX [MJ]: Ion thermal energy stored in the plasma core (WEX - Pedestal)
!  0.0016 * 3/2* Integral {0, R} (Nex*Tex-Nexb*Texb) dV
!   (Pereverzev 26-MAR-00)
double precision function WCEXR(YR)

use scalars, only: NA1
use status, only: NEX, TEX
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR
integer :: j
double precision, dimension(NA1) :: YY

do j=1, NA1
   YY(j) = NEX(J)*TEX(J) - NEX(NA1)*TEX(NA1)
enddo
WCEXR = 0.0024*VINT(YY, YR)

end function WCEXR
