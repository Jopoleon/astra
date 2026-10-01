! Q7TOT [[F9]/s]:  Integral {0, R} (SF7TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q7TOTR(YR)

use status, only: SF7TOT
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

Q7TOTR = VINT(SF7TOT, YR)

end function Q7TOTR
