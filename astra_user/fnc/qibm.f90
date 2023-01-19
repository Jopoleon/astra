! QIBM [MW]:  Integral {0, R} (PIBM) dV
!   (Pereverzev 20-MAY-08)
double precision FUNCTION QIBMR(YR)

use status_inc, only: PIBM

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QIBMR = VINT(PIBM, YR)

end function QIBMR
