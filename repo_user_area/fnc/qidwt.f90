! QIDWT [MW]: dWi/dt
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QIDWTR(YR)

use const_inc, only: TAU, HRO
use status_inc, only: NI, NIO, TI, TIO

implicit none

double precision, intent(in) :: YR

double precision :: YQ, YQO, VINT, VINTO

YQO = VINTO(NIO*TIO, YR)
YQ  = VINT (NI *TI , YR)

QIDWTR = (YQ - YQO)/TAU*.0024

end function QIDWTR
