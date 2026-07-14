! Q2TOT [[F9]/s]:  Integral {0, R} (SF2TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q2TOTR(YR)

use status, only: SF2TOT
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

Q2TOTR = VINT(SF2TOT, YR)

end function Q2TOTR
