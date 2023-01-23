! WTOTX [MJ]: Integral {0:R} ( 3/2*NI*TI ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION WIR(YR)

use status_inc, only: NI, TI

implicit none

double precision, intent(in) :: YR
double precision VINT

WIR = VINT(NI*TI, YR)
WIR = 0.0024*WIR

end function WIR
