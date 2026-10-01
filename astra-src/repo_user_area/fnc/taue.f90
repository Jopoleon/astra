! TAUE [s]: Full energy confinement time at radial position R [m]
!       TAUE=Wtot/(Qtot-dWtot/dt)
!   (Yushmanov 29-DEC-90)
double precision function TAUER(YR)

use scalars, only: TAU
use status, only: PITOT, PETOT, NI, NE, NIO, NEO, TI, TE, TIO, TEO
use standard_functions, only: VINT, VINTO

implicit none

double precision, intent(in) :: YR
double precision :: YQ, YW, YWO

YQ  = VINT(PITOT + PETOT, YR)
YW  = VINT(NI*TI + NE*TE, YR)
YWO = VINTO(NIO*TIO + NEO*TEO, YR)

if (YQ == 0.) then
   TAUER = 0.
else
   TAUER = TAU*YW/((YQ*TAU*417. + YWO) - YW)
endif

end function TAUER
