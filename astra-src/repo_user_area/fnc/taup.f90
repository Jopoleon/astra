! TAUP [s]: Particle confinement time at radial position R [m]
!       TAUP=Vint(Ne)/(Vint(SNTOT)-dVint(Ne)/dt)
!   (Yushmanov 11-MAY-87)
double precision function TAUPR(YR)

use scalars, only: TAU
use status, only: SNTOT, NE, NEO
use standard_functions, only: VINT, VINTO

implicit none

double precision, intent(in) :: YR
double precision :: YQ, YW, YWO

YQ  = VINT(SNTOT, YR)
YW  = VINT(NE, YR)
YWO = VINTO(NEO, YR)

if (YQ == 0.) then
   TAUPR = 0.
else
   TAUPR = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUPR
