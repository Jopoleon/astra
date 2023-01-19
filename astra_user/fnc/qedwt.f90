! QEDWT [MW]: dWi/dt
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QEDWTR(YR)

use const_inc, only: TAU, HRO
use status_inc, only: NE, NEO, TE, TEO

implicit none

double precision, intent(in) :: YR

double precision :: YQ, YQO, VINT, VINTO

YQO = VINTO(NEO*TEO, YR)
YQ  = VINT (NE *TE , YR)

QEDWTR = (YQ - YQO)/TAU*HRO*.0024

end function QEDWTR
