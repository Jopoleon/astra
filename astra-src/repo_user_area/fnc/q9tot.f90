! Q9TOT [[F9]/s]:  Integral {0, R} (SF9TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q9TOTR(YR)

use status, only: SF9TOT
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

Q9TOTR = VINT(SF9TOT, YR)

end function Q9TOTR
