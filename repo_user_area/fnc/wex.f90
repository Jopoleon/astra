! WEX [MJ]: Integral {0:R} ( 3/2*NEX*TEX ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION WEXR(YR)

use status, only: NEX, TEX
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

WEXR = VINT(NEX*TEX, YR)
WEXR = 0.0024*WEXR

end function WEXR
