! WIX [MJ]: Integral {0:R} ( 3/2*NIX*TIX ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION WIXR(YR)

use status_inc, only: NI, TIX

implicit none

double precision, intent(in) :: YR
double precision VINT

WIXR = VINT(NI*TIX, YR)
WIXR = 0.0024*WIXR

end function WIXR
