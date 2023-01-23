! QITOT [MW]:  Integral {0, R} (PITOT) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QITOTR(YR)

use status_inc, only: PITOT

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QITOTR = VINT(PITOT, YR)

end function QITOTR
