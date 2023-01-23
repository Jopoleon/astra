! WEX [MJ]: Integral {0:R} ( 3/2*NEX*TEX ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION WEXR(YR)

use status_inc, only: NEX, TEX

implicit none

double precision, intent(in) :: YR
double precision VINT

WEXR = VINT(NEX*TEX, YR)
WEXR = 0.0024*WEXR

end function WEXR
