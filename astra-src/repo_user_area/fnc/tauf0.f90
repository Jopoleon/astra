! TAUF0 [s]: Confinement_time@(rho) for the quantity F0
!       TAUF0=Vint(F0)/(Vint(SF0TOT)-dVint(F0)/dt)
!   (Pereverzev 09-OCT-08)
double precision function TAUF0R(YR)

use scalars, only: TAU
use status, only: SF0TOT, F0, F0O
use standard_functions, only: VINT, VINTO

implicit none

double precision, intent(in) :: YR
double precision :: YQ, YW, YWO

YQ  = VINT(SF0TOT, YR)
YW  = VINT(F0, YR)
YWO = VINTO(F0O, YR)

if (YQ == 0.) then
   TAUF0R = 0.
else
   TAUF0R = TAU*YW/((YQ*TAU + YWO) - YW)
endif

end function TAUF0R
