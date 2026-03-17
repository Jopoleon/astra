! NET [10#19/m#3]: Total number of electrons (r)
! Integral{0:R}(NE)dV
!   (Pereverzev 16-OCT-09)
double precision FUNCTION NETR(YR)

use status, only: NE
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

NETR = VINT(NE, YR)

end function NETR
