! NET [10#19/m#3]: Total number of electrons (r)
! Integral{0:R}(NE)dV
!   (Pereverzev 16-OCT-09)
double precision FUNCTION NETR(YR)

use status_inc, only: NE

implicit none

double precision, intent(in) :: YR
double precision :: VINT

NETR = VINT(NE, YR)

end function NETR
