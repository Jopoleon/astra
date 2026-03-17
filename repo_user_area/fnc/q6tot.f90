! Q6TOT [[F9]/s]:  Integral {0, R} (SF6TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q6TOTR(YR)

use status, only: SF6TOT
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

Q6TOTR = VINT(SF6TOT, YR)

end function Q6TOTR
