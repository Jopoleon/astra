! TAUEI [s]: Ion energy confinement time at radial position R [m]
!       TAUEI=Wi/(Qitot-dWi/dt)
!   (Yushmanov 11-MAY-87)
double precision function TAUEIR(YR)

use const_inc, only: TAU
use status_inc, only: PITOT, NI, NIO, TI, TIO

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VINTO, YQ, YW, YWO

YQ  = VINT(PITOT, YR)
YW  = VINT(NI*TI, YR)
YWO = VINTO(NIO*TIO, YR)

if(YQ == 0.) then
   TAUEIR = 0.
else
   TAUEIR = TAU*YW/((YQ*TAU*417. + YWO) - YW)
endif

end function TAUEIR
