! QNDNT [10#19/s]: d/dt(Volume integral {0, R}  NE )
!   (Yushmanov 15-FEB-89)
double precision FUNCTION QNDNTR(YR)

use const_inc, only: HRO, TAU
use status_inc, only: NE, NEO

implicit none

double precision, intent(in) :: YR
double precision YQ, YQO, VINT, VINTO

YQ  = VINT (NE , YR)
YQO = VINTO(NEO, YR)

QNDNTR = (YQ - YQO)/TAU*HRO

end function QNDNTR
