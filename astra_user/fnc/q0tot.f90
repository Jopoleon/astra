! Q0TOT [[F9]/s]:  Integral {0, R} (SF0TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q0TOTR(YR)

use status_inc, only: SF0TOT
implicit none

double precision, intent(in) :: YR
double precision :: VINT

Q0TOTR = VINT(SF0TOT, YR)

end function Q0TOTR
