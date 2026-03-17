! Q8TOT [[F9]/s]:  Integral {0, R} (SF8TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q8TOTR(YR)

use status, only: SF8TOT
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

Q8TOTR = VINT(SF8TOT, YR)

end function Q8TOTR
