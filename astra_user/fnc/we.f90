! WE [MJ]: Integral {0:R} ( 3/2*NE*TE ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION WER(YR)

use status_inc, only: NE, TE

implicit none

double precision, intent(in) :: YR
double precision VINT

WER = VINT(NE*TE, YR)
WER = 0.0024*WER

end function WER
