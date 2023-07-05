! QUDUT []: d/dt(Volume integral {0, R}  Upar*Ups )
!   (E Fable July 2023)
double precision FUNCTION QUDUTR(YR)

use const_inc, only: HRO, TAU
use status_inc, only: UPAR, UPARO, UPS0, UPS0O

implicit none

double precision, intent(in) :: YR
double precision YQ, YQO, VINT, VINTO

YQ  = VINT (UPAR*UPS0 , YR)
YQO = VINTO(UPARO*UPS0O, YR)

QUDUTR = (YQ - YQO)/TAU

end function QUDUTR
