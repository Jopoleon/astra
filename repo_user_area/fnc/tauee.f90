! TAUEE [s]: Electron energy confinement time at radial position R [m]
!       TAUEE=We/(Qetot-dWe/dt)
!   (Yushmanov 11-MAY-87)
double precision function TAUEER(YR)

use scalars, only: TAU
use status, only: PETOT, NE, NEO, TE, TEO
use standard_functions, only: VINT, VINTO

implicit none

double precision, intent(in) :: YR
double precision :: YQ, YW, YWO

YQ  = VINT(PETOT, YR)
YW  = VINT(NE*TE, YR)
YWO = VINTO(NEO*TEO, YR)

if (YQ == 0.) then
   TAUEER = 0.
else
   TAUEER = TAU*YW/((YQ*TAU*417. + YWO) - YW)
endif

end function TAUEER
