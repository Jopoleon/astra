! Q1TOT [[F9]/s]:  Integral {0, R} (SF1TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q1TOTR(YR)

use status, only: SF1TOT
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

Q1TOTR = VINT(SF1TOT, YR)

end function Q1TOTR
