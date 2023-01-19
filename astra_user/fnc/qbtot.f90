! QBTOT [MW]:  Integral {0, R} (PBEAM) dV
!   (Polevoy 28.09.89)
double precision FUNCTION QBTOTR(YR)

use status_inc, only: PBEAM

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QBTOTR = VINT(PBEAM, YR)

end function QBTOTR
