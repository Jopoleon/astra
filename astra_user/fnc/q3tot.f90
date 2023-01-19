! Q3TOT [[F9]/s]:  Integral {0, R} (SF3TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q3TOTR(YR)

use status_inc, only: SF3TOT
implicit none

double precision, intent(in) :: YR
double precision :: VINT

Q3TOTR = VINT(SF3TOT, YR)

end function Q3TOTR
