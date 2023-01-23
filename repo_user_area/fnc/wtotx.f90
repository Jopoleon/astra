! WTOTX [MJ]: Integral {0:R} ( 3/2*NEX*(TEX+TIX) ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION WTOTXR(YR)

use status_inc, only: NEX, TEX, NI, TIX

implicit none

double precision, intent(in) :: YR
double precision VINT

WTOTXR = VINT(NEX*TEX, YR) + VINT(NI*TIX, YR)
WTOTXR = 0.0024*WTOTXR

end function WTOTXR
