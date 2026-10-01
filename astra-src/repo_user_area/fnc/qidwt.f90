! QIDWT [MW]: dWi/dt
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QIDWTR(YR)

use scalars, only: TAU, HRO
use status, only: NI, NIO, TI, TIO
use standard_functions, only: VINT, VINTO

implicit none

double precision, intent(in) :: YR

double precision :: YQ, YQO

YQO = VINTO(NIO*TIO, YR)
YQ  = VINT (NI *TI , YR)
QIDWTR = (YQ - YQO)/TAU*.0024

end function QIDWTR
