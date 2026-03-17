! Q5TOT [[F9]/s]:  Integral {0, R} (SF5TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q5TOTR(YR)

use status, only: SF5TOT
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

Q5TOTR = VINT(SF5TOT, YR)

end function Q5TOTR
