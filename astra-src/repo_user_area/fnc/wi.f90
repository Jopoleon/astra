! WTOTX [MJ]: Integral {0:R} ( 3/2*NI*TI ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION WIR(YR)

use status, only: NI, TI
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

WIR = VINT(NI*TI, YR)
WIR = 0.0024*WIR

end function WIR
