! TAUP [s]: Particle confinement time at radial position R [m]
!       TAUP=Vint(Ne)/(Vint(SNTOT)-dVint(Ne)/dt)
!   (Yushmanov 11-MAY-87)
double precision function TAUPR(YR)

use const_inc, only: TAU
use status_inc, only: SNTOT, NE, NEO

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VINTO, YQ, YW, YWO

YQ  = VINT(SNTOT, YR)
YW  = VINT(NE, YR)
YWO = VINTO(NEO, YR)

if (YQ == 0.) then
   TAUPR = 0.
else
   TAUPR = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUPR
