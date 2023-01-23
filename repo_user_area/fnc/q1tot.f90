! Q1TOT [[F9]/s]:  Integral {0, R} (SF1TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q1TOTR(YR)

use status_inc, only: SF1TOT
implicit none

double precision, intent(in) :: YR
double precision :: VINT

Q1TOTR = VINT(SF1TOT, YR)

end function Q1TOTR
