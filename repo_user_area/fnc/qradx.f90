! QRADX [MW]:  Integral {0, R} (PRADX) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QRADXR(YR)

use status_inc, only: PRADX

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QRADXR = VINT(PRADX, YR)

end function QRADXR
