! Q2TOT [[F9]/s]:  Integral {0, R} (SF2TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q2TOTR(YR)

use status_inc, only: SF2TOT
implicit none

double precision, intent(in) :: YR
double precision :: VINT

Q2TOTR = VINT(SF2TOT, YR)

end function Q2TOTR
