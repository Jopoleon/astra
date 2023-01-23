! Q4TOT [[F9]/s]:  Integral {0, R} (SF4TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q4TOTR(YR)

use status_inc, only: SF4TOT
implicit none

double precision, intent(in) :: YR
double precision :: VINT

Q4TOTR = VINT(SF4TOT, YR)

end function Q4TOTR
