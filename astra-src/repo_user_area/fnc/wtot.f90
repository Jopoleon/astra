! WTOT [MJ]: Integral {0:R} ( 3/2*(NE*TE+NI*TI) ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION WTOTR(YR)

use status, only: NE, TE, NI, TI
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

WTOTR = VINT(NE*TE, YR) + VINT(NI*TI, YR)
WTOTR = 0.0024*WTOTR

end function WTOTR
