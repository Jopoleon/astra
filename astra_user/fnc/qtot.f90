! QTOT [MW]:  Integral {0, R} ( PETOT+PITOT ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QTOTR(YR)

use status_inc, only: PE, PI

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QTOTR = VINT(PE + PI, YR)

end function QTOTR
