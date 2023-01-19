! Q5TOT [[F9]/s]:  Integral {0, R} (SF5TOT) dV
!   (Pereverzev 09-May-2008)
double precision FUNCTION Q5TOTR(YR)

use status_inc, only: SF5TOT
implicit none

double precision, intent(in) :: YR
double precision :: VINT

Q5TOTR = VINT(SF5TOT, YR)

end function Q5TOTR
