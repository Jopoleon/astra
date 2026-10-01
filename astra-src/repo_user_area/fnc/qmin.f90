! QMIN []: min(q). Minimal value of safety factor "q" on [0, a]
!   (Pereverzev 13-JAN-97)
! Usage: qmin_QMINB;
!
double precision FUNCTION QMINR(YR)

use scalars, only: NA1
use status, only: MU

implicit none

double precision, intent(in) :: YR
QMINR = MINVAL(1./MU(1:NA1))

end function QMINR
