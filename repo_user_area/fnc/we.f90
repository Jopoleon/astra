! WE [MJ]: Integral {0:R} ( 3/2*NE*TE ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION WER(YR)

use status, only: NE, TE
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

WER = VINT(NE*TE, YR)
WER = 0.0024*WER

end function WER
