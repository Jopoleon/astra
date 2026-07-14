! QEDWT [MW]: dWi/dt
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QEDWTR(YR)

use scalars, only: TAU, HRO
use status, only: NE, NEO, TE, TEO
use standard_functions, only: VINT, VINTO

implicit none

double precision, intent(in) :: YR
double precision :: YQ, YQO

YQO = VINTO(NEO*TEO, YR)
YQ  = VINT (NE *TE , YR)
QEDWTR = (YQ - YQO)/TAU*.0024

end function QEDWTR
