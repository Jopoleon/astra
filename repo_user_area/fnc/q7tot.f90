! Q7TOT [[F9]/s]:  Integral {0, R} (SF7TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q7TOTR(YR)

use status_inc, only: SF7TOT
implicit none

double precision, intent(in) :: YR
double precision :: VINT

Q7TOTR = VINT(SF7TOT, YR)

end function Q7TOTR
