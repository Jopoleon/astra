! TAUEE [s]: Electron energy confinement time at radial position R [m]
!       TAUEE=We/(Qetot-dWe/dt)
!   (Yushmanov 11-MAY-87)
double precision function TAUEER(YR)

use const_inc, only: TAU
use status_inc, only: PETOT, NE, NEO, TE, TEO

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VINTO, YQ, YW, YWO

YQ  = VINT(PETOT, YR)
YW  = VINT(NE*TE, YR)
YWO = VINTO(NEO*TEO, YR)

if(YQ == 0.) then
   TAUEER = 0.
else
   TAUEER = TAU*YW/((YQ*TAU*417. + YWO) - YW)
endif

end function TAUEER
